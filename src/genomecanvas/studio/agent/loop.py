"""The chat agent: one user message in, tool calls and figures out, the same loop for every model.

Events sent to the browser (``emit``): ``text`` (streamed reply), ``tool`` (a tool started),
``job`` (a figure was queued), ``figure`` (a figure is ready), ``error``, ``done``.
"""

from __future__ import annotations

import time
import traceback
from pathlib import Path

from ..db import new_id
from .providers import ProviderError, adapter
from .settings import agent_settings
from .tools import TOOLS, Toolbox, load_image

PREFACE = """You are the GenomeCanvas Studio assistant for the project "{name}" (genome: {genome}).
You help biologists turn this project's sequencing data into figures. Reply in the user's language.

You act only through tools. The guide below was written for a command line; here:
* `gcplot.py survey` is the survey_data tool (it covers this project's files only);
* `gcplot.py run` is run_plot, with the same request JSON; refer to files by the names or labels
  survey_data returned, never by paths, and never invent file names;
* the family skills are available through read_skill.
Each figure appears in the chat as soon as it is rendered, so do not paste links or file paths;
refer to figures by what they show. You cannot change, upload or delete files.
Keep messages short. Ask only the questions whose answers change the figure.
Any text you write, including before a tool call, is shown to the user as it streams: write it for
them (at most one short sentence on what you are about to do) or write nothing.
Keep what the user explicitly asked for (number of clusters, samples to drop, scale, region set).
If a figure suggests a different setting would be better, say so and ask; do not apply it yourself.

"""


class ClientGone(Exception):
    """The browser closed the stream (the user pressed stop or left)."""


class ChatAgent:
    def __init__(self, app):
        self.app = app
        self.settings = agent_settings(app.config.agent)
        self.guide = self._guide()

    def _guide(self):
        path = Path(self.app.config.static_dir).parent / "skills" / "genomecanvas" / "SKILL.md"
        if not path.is_file():
            return ""
        text = path.read_text(encoding="utf-8")
        return text.split("---", 2)[2].strip() if text.startswith("---") else text

    def system_prompt(self, project):
        preface = PREFACE.format(name=project["name"], genome=project.get("genome") or "unknown")
        return preface + self.guide

    # ------------------------------------------------------------------ conversations
    def conversations(self, user, project_id):
        return self.app.db.all(
            "SELECT id, title, provider, model, created, updated FROM agent_conversations "
            "WHERE user_id = ? AND project_id = ? ORDER BY updated DESC LIMIT 50",
            user["id"],
            project_id,
        )

    def conversation(self, user, project_id, conversation_id):
        row = self.app.db.one(
            "SELECT * FROM agent_conversations WHERE id = ? AND user_id = ? AND project_id = ?",
            conversation_id,
            user["id"],
            project_id,
        )
        if not row:
            raise ValueError("Conversation not found")
        return row

    def used_today(self, user):
        row = self.app.db.one(
            "SELECT tokens FROM agent_usage WHERE user_id = ? AND day = ?",
            user["id"],
            time.strftime("%Y-%m-%d"),
        )
        return row["tokens"] if row else 0

    def _count(self, user, tokens):
        day = time.strftime("%Y-%m-%d")
        self.app.db.run(
            "INSERT INTO agent_usage (user_id, day, tokens) VALUES (?, ?, ?) "
            "ON CONFLICT(user_id, day) DO UPDATE SET tokens = tokens + excluded.tokens",
            user["id"],
            day,
            tokens,
        )

    # ------------------------------------------------------------------ one turn
    def run(self, user, project, text, provider_name, model, emit, conversation_id=None, ip=None):
        provider = self.settings.provider(project["id"], provider_name, model)
        limit = self.settings.daily_tokens_per_user
        if limit and self.used_today(user) >= limit:
            raise ValueError("You have used today's assistant allowance; try again tomorrow")
        text = str(text or "").strip()
        if not text:
            raise ValueError("Type a message first")
        now = time.time()
        if conversation_id:
            row = self.conversation(user, project["id"], conversation_id)
            messages = row["messages"] or []
        else:
            conversation_id = new_id()
            messages = []
            self.app.db.insert(
                "agent_conversations",
                id=conversation_id,
                project_id=project["id"],
                user_id=user["id"],
                title=text[:80],
                provider=provider.name,
                model=model,
                messages=[],
                usage={},
                created=now,
                updated=now,
            )
        emit(dict(type="conversation", id=conversation_id))
        messages.append({"role": "user", "text": text})
        tools = Toolbox(self.app, user, project, ip, emit)
        llm = adapter(provider)
        system = self.system_prompt(project)
        totals = dict(input=0, output=0)
        used_tools = []
        stopped, failure = False, None

        def save():
            self.app.db.update(
                "agent_conversations",
                "id",
                conversation_id,
                messages=messages,
                provider=provider.name,
                model=model,
                updated=time.time(),
            )

        try:
            for _ in range(self.settings.max_steps):
                reply = llm.complete(
                    model,
                    system,
                    self._with_images(messages),
                    TOOLS,
                    lambda chunk: emit(dict(type="text", text=chunk)),
                    session=conversation_id,
                )
                totals["input"] += reply.usage["input"]
                totals["output"] += reply.usage["output"]
                self._count(user, reply.usage["input"] + reply.usage["output"])
                messages.append(
                    dict(role="assistant", text=reply.text, calls=reply.calls, raw=reply.raw)
                )
                save()
                if reply.stop == "refusal":
                    emit(dict(type="error", message="The model declined this request."))
                    break
                if not reply.calls:
                    break
                if reply.stop == "max_tokens":
                    self._close_calls(messages, reply.calls, "The reply was cut off; not run.")
                    emit(dict(type="error", message="The reply was too long and was cut off."))
                    break
                results = []
                for call in reply.calls:
                    emit(dict(type="tool", name=call["name"], input=call["input"]))
                    used_tools.append(call["name"])
                    output, image, failed = tools.call(call["name"], call["input"])
                    result = dict(id=call["id"], name=call["name"], text=output, is_error=failed)
                    if image:
                        result["image_ref"] = image
                    results.append(result)
                messages.append(dict(role="tool", results=results))
                save()
            else:
                emit(dict(type="error", message="Stopped after the maximum number of steps."))
        except ClientGone:
            stopped = True
            self._close_open(messages, "The user stopped the reply.")
            save()
        except ProviderError as exc:
            self._close_open(messages, "The model call failed.")
            save()
            failure = str(exc)
        except Exception:  # noqa: BLE001 - keep the conversation usable whatever went wrong
            traceback.print_exc()
            self._close_open(messages, "The assistant hit an internal error.")
            save()
            failure = "The assistant hit an internal error. The conversation was saved; try again."
        self.app.db.audit(
            user,
            "agent.turn",
            project["id"],
            dict(
                provider=provider.name,
                model=model,
                tokens=totals,
                tools=used_tools,
                stopped=stopped,
                error=failure,
            ),
            ip,
        )
        if stopped:
            raise ClientGone
        if failure:
            emit(dict(type="error", message=failure))
        emit(dict(type="done", conversation=conversation_id, usage=totals))
        return conversation_id

    def _with_images(self, messages):
        """Copy of the history with figure images loaded (they are stored as job references)."""
        out = []
        for m in messages:
            if m["role"] == "tool" and any(r.get("image_ref") for r in m["results"]):
                results = []
                for r in m["results"]:
                    r = dict(r)
                    ref = r.pop("image_ref", None)
                    if ref:
                        image = load_image(self.app, ref)
                        if image:
                            r["image"] = image
                        else:
                            r["text"] += " (the figure has since expired)"
                    results.append(r)
                m = dict(m, results=results)
            out.append(m)
        return out

    @staticmethod
    def _close_calls(messages, calls, reason):
        messages.append(
            dict(
                role="tool",
                results=[
                    dict(id=c["id"], name=c["name"], text=reason, is_error=True) for c in calls
                ],
            )
        )

    def _close_open(self, messages, reason):
        """Every tool call needs a result before the next request; fill any that are missing."""
        if messages and messages[-1]["role"] == "assistant" and messages[-1].get("calls"):
            self._close_calls(messages, messages[-1]["calls"], reason)

    # ------------------------------------------------------------------ display
    @staticmethod
    def transcript(row):
        """The conversation as the chat window shows it: text, tool steps and figures."""
        items = []
        for m in row["messages"] or []:
            if m["role"] == "user":
                items.append(dict(type="user", text=m["text"]))
            elif m["role"] == "assistant":
                if m["text"]:
                    items.append(dict(type="assistant", text=m["text"]))
                items += [dict(type="tool", name=c["name"]) for c in m.get("calls", [])]
            else:
                for r in m["results"]:
                    ref = r.get("image_ref")
                    if ref:
                        items.append(
                            dict(
                                type="figure",
                                job=ref["job"],
                                url=f"/results/{ref['job']}/figure.png",
                            )
                        )
        return items
