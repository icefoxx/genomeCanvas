"""Model adapters: one internal conversation format, translated per vendor API.

Internal messages (stored in the database, append-only):

* ``{"role": "user", "text": str}``
* ``{"role": "assistant", "text": str, "calls": [{"id", "name", "input"}], "raw": {...}}`` where
  ``raw`` keeps the vendor's own content blocks so they can be sent back unchanged
* ``{"role": "tool", "results": [{"id", "name", "text", "is_error", "image"?}]}`` where
  ``image`` is ``{"media_type", "data"}`` (base64) filled in by the loop just before a request

Each adapter's ``complete()`` streams text through ``on_text`` and returns a :class:`Reply`.
"""

from __future__ import annotations

import json
import ssl
import urllib.error
import urllib.request
from dataclasses import dataclass, field

from .settings import Provider


class ProviderError(Exception):
    """A model call failed; the message is safe to show to the user."""


@dataclass
class Reply:
    text: str
    calls: list[dict]
    stop: str  # end | tool_use | max_tokens | refusal
    raw: dict | None = None
    usage: dict = field(default_factory=lambda: dict(input=0, output=0))
    model: str = ""


def adapter(provider: Provider):
    kinds = {"anthropic": AnthropicAdapter, "openai": OpenAIAdapter, "opencode": OpenAIAdapter}
    return kinds[provider.kind](provider)


USER_AGENT = "GenomeCanvas-Studio"  # urllib's default agent is blocked by some gateways' bot rules
TLS_HINT = (
    "the server does not trust the TLS certificate; if this network inspects HTTPS "
    '(e.g. Cloudflare Gateway), set ca_bundle = "system" for this provider'
)


def ssl_context(provider: Provider):
    """TLS context for a provider's ca_bundle setting; None means Python's defaults."""
    if not provider.ca_bundle:
        return None
    if provider.ca_bundle == "system":
        try:
            import truststore
        except ImportError as exc:
            raise ProviderError('ca_bundle = "system" needs the truststore package') from exc
        return truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    return ssl.create_default_context(cafile=provider.ca_bundle)


def request_headers(provider: Provider, session=None) -> dict:
    """A gateway's extra headers, plus the conversation id where the gateway asks for one."""
    headers = {str(k): str(v) for k, v in provider.headers.items()}
    if provider.session_header and session:
        headers[provider.session_header] = session
    return headers


# ----------------------------------------------------------------------------- Anthropic


class AnthropicAdapter:
    """Claude through the official SDK; a custom base_url (local gateway) gets a plain request."""

    def __init__(self, provider: Provider):
        self.p = provider

    def client(self):
        try:
            import anthropic
        except ImportError as exc:
            raise ProviderError("The anthropic package is not installed on the server") from exc
        kwargs = dict(max_retries=2)
        if self.p.api_key():
            kwargs["auth_token" if self.p.auth == "bearer" else "api_key"] = self.p.api_key()
        if self.p.base_url:
            kwargs["base_url"] = self.p.base_url
        context = ssl_context(self.p)
        if context is not None:
            kwargs["http_client"] = anthropic.DefaultHttpxClient(verify=context)
        return anthropic.Anthropic(**kwargs)

    def wire(self, messages):
        out = []
        for m in messages:
            if m["role"] == "user":
                out.append({"role": "user", "content": [{"type": "text", "text": m["text"]}]})
            elif m["role"] == "assistant":
                raw = m.get("raw") or {}
                if raw.get("provider") == self.p.name and raw.get("content"):
                    content = _echo_after_fallback(raw["content"])
                else:
                    content = [{"type": "text", "text": m["text"]}] if m["text"] else []
                    content += [
                        {"type": "tool_use", "id": c["id"], "name": c["name"], "input": c["input"]}
                        for c in m.get("calls", [])
                    ]
                out.append({"role": "assistant", "content": content})
            else:
                blocks = []
                for r in m["results"]:
                    content = [{"type": "text", "text": r["text"]}]
                    if r.get("image") and self.p.vision:
                        content.append(
                            {"type": "image", "source": {"type": "base64", **r["image"]}}
                        )
                    blocks.append(
                        {
                            "type": "tool_result",
                            "tool_use_id": r["id"],
                            "content": content,
                            "is_error": bool(r.get("is_error")),
                        }
                    )
                out.append({"role": "user", "content": blocks})
        return out

    def complete(self, model, system, messages, tools, on_text, session=None) -> Reply:
        import anthropic

        client = self.client()
        kwargs = dict(
            model=model,
            max_tokens=self.p.max_tokens,
            system=[{"type": "text", "text": system}],
            messages=self.wire(messages),
            tools=[
                {"name": t["name"], "description": t["description"], "input_schema": t["schema"]}
                for t in tools
            ],
        )
        if self.p.effort:
            kwargs["output_config"] = {"effort": self.p.effort}
        extra = request_headers(self.p, session)
        if extra:
            kwargs["extra_headers"] = extra
        if self.p.official:
            # Cache the stable prefix (system, tools, history) across steps and turns, and let the
            # API re-run a policy decline on a suitable model instead of stopping the chat.
            kwargs["cache_control"] = {"type": "ephemeral"}
            messages_api = client.beta.messages
            kwargs.update(betas=["server-side-fallback-2026-07-01"], fallbacks="default")
        else:
            messages_api = client.messages
        streamed = []
        try:
            try:
                with messages_api.stream(**kwargs) as events:
                    for event in events:
                        if event.type == "text":
                            streamed.append(event.text)
                            on_text(event.text)
                    message = events.get_final_message()
            except (anthropic.APIStatusError, anthropic.APIConnectionError):
                raise
            except Exception:  # noqa: BLE001 - a stream that ended without a final message
                # Ask once more without streaming; only show text the user has not seen yet.
                message = messages_api.create(**kwargs)
                full = "".join(getattr(b, "text", "") for b in message.content if b.type == "text")
                seen = "".join(streamed)
                if full.startswith(seen) and len(full) > len(seen):
                    on_text(full[len(seen) :])
        except anthropic.APIStatusError as exc:
            raise ProviderError(f"{self.p.name}: {exc.message}") from exc
        except anthropic.APIConnectionError as exc:
            reason = (
                TLS_HINT
                if "CERTIFICATE_VERIFY_FAILED" in repr(exc.__cause__)
                else "cannot reach the model service"
            )
            raise ProviderError(f"{self.p.name}: {reason}") from exc
        content = [b.model_dump(mode="json", exclude_none=True) for b in message.content]
        usage = message.usage
        stop = {"tool_use": "tool_use", "max_tokens": "max_tokens", "refusal": "refusal"}.get(
            message.stop_reason, "end"
        )
        return Reply(
            text="".join(b.get("text", "") for b in content if b.get("type") == "text"),
            calls=[
                dict(id=b["id"], name=b["name"], input=b.get("input") or {})
                for b in content
                if b.get("type") == "tool_use"
            ],
            stop=stop,
            raw=dict(provider=self.p.name, model=message.model, content=content),
            usage=dict(
                input=(usage.input_tokens or 0)
                + (getattr(usage, "cache_read_input_tokens", 0) or 0)
                + (getattr(usage, "cache_creation_input_tokens", 0) or 0),
                output=usage.output_tokens or 0,
            ),
            model=message.model,
        )


def _echo_after_fallback(content):
    """After a mid-output fallback, drop model-internal blocks before the last fallback marker."""
    marks = [i for i, b in enumerate(content) if b.get("type") == "fallback"]
    if not marks:
        return content
    internal = {"thinking", "redacted_thinking", "tool_use"}
    return [b for i, b in enumerate(content) if i > marks[-1] or b.get("type") not in internal]


# ----------------------------------------------------------------------------- OpenAI-compatible

OPENAI_URL = "https://api.openai.com/v1"


class OpenAIAdapter:
    """Chat Completions over HTTP: OpenAI and compatible servers (vLLM, Ollama, LM Studio)."""

    def __init__(self, provider: Provider):
        self.p = provider

    def wire(self, system, messages):
        out = [{"role": "system", "content": system}]
        for m in messages:
            if m["role"] == "user":
                out.append({"role": "user", "content": m["text"]})
            elif m["role"] == "assistant":
                item = {"role": "assistant", "content": m["text"] or None}
                if m.get("calls"):
                    item["tool_calls"] = [
                        {
                            "id": c["id"],
                            "type": "function",
                            "function": {"name": c["name"], "arguments": json.dumps(c["input"])},
                        }
                        for c in m["calls"]
                    ]
                out.append(item)
            else:
                images = []
                for r in m["results"]:
                    text = ("ERROR: " if r.get("is_error") else "") + r["text"]
                    out.append({"role": "tool", "tool_call_id": r["id"], "content": text})
                    if r.get("image") and self.p.vision:
                        images.append(r["image"])
                if images:  # tool messages carry text only; figures follow as a user message
                    out.append(
                        {
                            "role": "user",
                            "content": [{"type": "text", "text": "Figure(s) from the tools above:"}]
                            + [
                                {
                                    "type": "image_url",
                                    "image_url": {
                                        "url": f"data:{i['media_type']};base64,{i['data']}"
                                    },
                                }
                                for i in images
                            ],
                        }
                    )
        return out

    def complete(self, model, system, messages, tools, on_text, session=None) -> Reply:
        body = dict(
            model=model,
            messages=self.wire(system, messages),
            tools=[
                {
                    "type": "function",
                    "function": {
                        "name": t["name"],
                        "description": t["description"],
                        "parameters": t["schema"],
                    },
                }
                for t in tools
            ],
            stream=True,
        )
        # OpenAI's own API takes max_completion_tokens; compatible servers take max_tokens.
        body["max_completion_tokens" if self.p.official else "max_tokens"] = self.p.max_tokens
        if self.p.stream_usage:
            body["stream_options"] = {"include_usage": True}
        headers = {
            "Content-Type": "application/json",
            "Accept": "text/event-stream",
            "User-Agent": USER_AGENT,
            **request_headers(self.p, session),
        }
        if self.p.api_key():
            headers["Authorization"] = f"Bearer {self.p.api_key()}"
        url = (self.p.base_url or OPENAI_URL).rstrip("/") + "/chat/completions"
        request = urllib.request.Request(url, json.dumps(body).encode(), headers, method="POST")
        text, calls, finish, usage, served = [], {}, None, dict(input=0, output=0), model
        try:
            with urllib.request.urlopen(
                request, timeout=600, context=ssl_context(self.p)
            ) as response:
                for raw_line in response:
                    line = raw_line.decode("utf-8", "replace").strip()
                    if not line.startswith("data:"):
                        continue
                    data = line[5:].strip()
                    if data == "[DONE]":
                        break
                    chunk = json.loads(data)
                    served = chunk.get("model") or served
                    if chunk.get("usage"):
                        usage = dict(
                            input=chunk["usage"].get("prompt_tokens") or 0,
                            output=chunk["usage"].get("completion_tokens") or 0,
                        )
                    for choice in chunk.get("choices") or []:
                        delta = choice.get("delta") or {}
                        if delta.get("content"):
                            text.append(delta["content"])
                            on_text(delta["content"])
                        for part in delta.get("tool_calls") or []:
                            slot = calls.setdefault(
                                part.get("index", 0), dict(id="", name="", args="")
                            )
                            slot["id"] = part.get("id") or slot["id"]
                            function = part.get("function") or {}
                            slot["name"] += function.get("name") or ""
                            slot["args"] += function.get("arguments") or ""
                        finish = choice.get("finish_reason") or finish
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode("utf-8", "replace")[:500]
            try:
                detail = json.loads(detail).get("error", {}).get("message") or detail
            except (ValueError, AttributeError):
                pass
            raise ProviderError(f"{self.p.name}: {detail}") from exc
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            if isinstance(getattr(exc, "reason", exc), ssl.SSLCertVerificationError):
                raise ProviderError(f"{self.p.name}: {TLS_HINT}") from exc
            raise ProviderError(f"{self.p.name}: cannot reach the model service ({exc})") from exc
        parsed = []
        for index in sorted(calls):
            slot = calls[index]
            try:
                arguments = json.loads(slot["args"] or "{}")
            except json.JSONDecodeError:
                arguments = {"_invalid_json": slot["args"][:2000]}
            parsed.append(
                dict(id=slot["id"] or f"call_{index}", name=slot["name"], input=arguments)
            )
        stop = {"length": "max_tokens", "content_filter": "refusal"}.get(finish)
        stop = stop or ("tool_use" if parsed else "end")
        return Reply(text="".join(text), calls=parsed, stop=stop, usage=usage, model=served)
