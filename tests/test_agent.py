"""Studio chat agent: settings, the tool loop over HTTP (scripted model), both model adapters."""

from __future__ import annotations

import json
import threading
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

pytest.importorskip("pyBigWig")

from test_studio import Client, data_root, make_config  # noqa: E402,F401 - fixture import

from genomecanvas.studio.agent import loop as agent_loop  # noqa: E402
from genomecanvas.studio.agent.providers import (  # noqa: E402
    AnthropicAdapter,
    OpenAIAdapter,
    Reply,
    _echo_after_fallback,
    ssl_context,
)
from genomecanvas.studio.agent.settings import Provider, agent_settings  # noqa: E402
from genomecanvas.studio.server import Studio, make_handler  # noqa: E402

AGENT = dict(
    max_steps=6,
    daily_tokens_per_user=0,
    local_only_projects=["P2"],
    providers=[
        dict(name="cloud", kind="anthropic", models=["claude-opus-5-5"]),
        dict(name="lab", kind="openai", base_url="http://gpu:8000/v1", models=["m"], local=True),
    ],
)


def test_settings_validation_and_local_only():
    settings = agent_settings(AGENT)
    assert [p.name for p in settings.available("P1")] == ["cloud", "lab"]
    assert [p.name for p in settings.available("P2")] == ["lab"]
    with pytest.raises(ValueError, match="not available"):
        settings.provider("P2", "cloud", "claude-opus-5-5")
    with pytest.raises(ValueError, match="does not offer"):
        settings.provider("P1", "cloud", "gpt-x")
    with pytest.raises(ValueError, match="kind"):
        agent_settings(dict(providers=[dict(name="x", kind="cohere", models=["m"])]))
    with pytest.raises(ValueError, match="unknown keys"):
        agent_settings(dict(providers=[dict(name="x", kind="openai", models=["m"], key="sk")]))
    assert not agent_settings({}).enabled
    keyed = agent_settings(
        dict(
            providers=[
                dict(name="k", kind="openai", models=["m"], api_key_env="GC_TEST_NO_SUCH_KEY")
            ]
        )
    )
    assert keyed.available("P1") == [] and keyed.missing_keys() == ["k (GC_TEST_NO_SUCH_KEY)"]


class Scripted:
    """A model that follows a script and records what it was sent."""

    def __init__(self, steps):
        self.steps, self.seen = list(steps), []

    def __call__(self, provider):
        return self

    def complete(self, model, system, messages, tools, on_text, session=None):
        self.seen.append(
            dict(
                session=session, system=system, messages=messages, tools=[t["name"] for t in tools]
            )
        )
        step = self.steps.pop(0)(messages)
        if step.text:
            on_text(step.text)
        return step


@pytest.fixture()
def agent_studio(tmp_path, data_root):  # noqa: F811 - pytest fixture
    app = Studio(make_config(tmp_path, data_root, agent=AGENT))
    app.start()
    app.db.create_user("alice@lab.org", password="alice-password1")
    app.db.create_user("bob@lab.org", password="bob-password12")
    for project, folder in (("P1", "smith/P1"), ("P2", "jones/P2")):
        app.projects.create(project, data_dir=folder)
        app.projects.scan(project)
    alice, bob = (app.db.user_by_email(e) for e in ("alice@lab.org", "bob@lab.org"))
    app.projects.set_member("P1", alice["id"], "editor")
    app.projects.set_member("P2", alice["id"], "editor")
    app.projects.set_member("P1", bob["id"], "viewer")
    server = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(app))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    yield app, f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()
    app.jobs.stop()


def chat(client, project, body):
    """POST a chat message and collect the Server-Sent Events."""
    request = urllib.request.Request(
        f"{client.base}/api/projects/{project}/agent",
        data=json.dumps(body).encode(),
        method="POST",
        headers={"X-GenomeCanvas": "studio", "Content-Type": "application/json"},
    )
    with client.opener.open(request, timeout=180) as response:
        assert response.headers.get_content_type() == "text/event-stream"
        return [
            json.loads(line[5:])
            for line in response.read().decode().splitlines()
            if line.startswith("data:")
        ]


def test_chat_turn_surveys_plots_and_saves(agent_studio, monkeypatch):
    app, base = agent_studio
    survey_seen = {}

    def ask_survey(messages):
        return Reply(
            "Let me look at your data.", [dict(id="t1", name="survey_data", input={})], "tool_use"
        )

    def plot_first_recommendation(messages):
        report = json.loads(messages[-1]["results"][0]["text"])
        survey_seen.update(report)
        request = next(r["request"] for r in report["recommendations"] if r["request"])
        return Reply("", [dict(id="t2", name="run_plot", input=dict(request=request))], "tool_use")

    def describe(messages):
        return Reply("Here is the profile heatmap.", [], "end", usage=dict(input=1000, output=50))

    model = Scripted([ask_survey, plot_first_recommendation, describe])
    monkeypatch.setattr(agent_loop, "adapter", model)
    alice = Client(base).login("alice@lab.org", "alice-password1")
    events = chat(
        alice, "P1", dict(message="Plot my data", provider="cloud", model="claude-opus-5-5")
    )
    kinds = [e["type"] for e in events]
    assert kinds[0] == "conversation" and kinds[-1] == "done", kinds
    assert "error" not in kinds, events
    assert [e["name"] for e in events if e["type"] == "tool"] == ["survey_data", "run_plot"]
    figure = next(e for e in events if e["type"] == "figure")
    assert alice.call("GET", figure["url"])[0] == 200
    # The survey shows project file names, never server paths.
    assert str(app.config.data_root) not in json.dumps(survey_seen)
    assert survey_seen["summary"]["signals"] == 4
    # The model got the rendered figure back as an image with the tool result.
    tool_turn = model.seen[2]["messages"][-1]
    assert (
        tool_turn["role"] == "tool"
        and tool_turn["results"][0]["image"]["media_type"] == "image/png"
    )
    assert "Smith lab ChIP-seq" in model.seen[0]["system"]
    # Saved, listed and replayable; usage counted; audited.
    conversation = events[0]["id"]
    assert {seen["session"] for seen in model.seen} == {conversation}  # gateways get the chat id
    status, listing = alice.call("GET", "/api/projects/P1/agent")
    assert status == 200 and listing["conversations"][0]["id"] == conversation
    assert [p["name"] for p in listing["providers"]] == ["cloud", "lab"]
    assert listing["used_today"] == 1050
    status, replay = alice.call("GET", f"/api/projects/P1/agent/{conversation}")
    assert [i["type"] for i in replay["items"]] == [
        "user",
        "assistant",
        "tool",
        "tool",
        "figure",
        "assistant",
    ]
    assert app.db.one("SELECT * FROM audit WHERE action = 'agent.turn'")
    # A follow-up continues the same conversation.
    monkeypatch.setattr(agent_loop, "adapter", Scripted([lambda m: Reply("Sure.", [], "end")]))
    events = chat(
        alice,
        "P1",
        dict(
            message="Thanks", provider="cloud", model="claude-opus-5-5", conversation=conversation
        ),
    )
    assert events[-1]["conversation"] == conversation
    assert alice.call("DELETE", f"/api/projects/P1/agent/{conversation}")[0] == 200


def test_chat_access_rules(agent_studio, monkeypatch):
    app, base = agent_studio
    monkeypatch.setattr(agent_loop, "adapter", Scripted([lambda m: Reply("hi", [], "end")] * 3))
    body = dict(message="hello", provider="cloud", model="claude-opus-5-5")
    alice = Client(base).login("alice@lab.org", "alice-password1")
    bob = Client(base).login("bob@lab.org", "bob-password12")
    assert bob.call("POST", "/api/projects/P1/agent", body)[0] == 403  # viewers cannot run jobs
    status, error = alice.call("POST", "/api/projects/P2/agent", body)  # local-only project
    assert status == 400 and "not available" in error["error"]
    status, listing = alice.call("GET", "/api/projects/P2/agent")
    assert [p["name"] for p in listing["providers"]] == ["lab"]
    conversation = chat(alice, "P1", body)[0]["id"]
    assert bob.call("GET", f"/api/projects/P1/agent/{conversation}")[0] == 400  # not his
    app.agent.settings.daily_tokens_per_user = 1
    app.agent._count(app.db.user_by_email("alice@lab.org"), 5)
    assert alice.call("POST", "/api/projects/P1/agent", body)[0] == 429


def test_tool_errors_go_back_to_the_model(agent_studio, monkeypatch):
    app, base = agent_studio
    bad = dict(plot="enrichment/profile", files=[{"name": "ghost.bw"}], settings={})
    model = Scripted(
        [
            lambda m: Reply(
                "", [dict(id="t1", name="run_plot", input=dict(request=bad))], "tool_use"
            ),
            lambda m: Reply("That file does not exist.", [], "end"),
        ]
    )
    monkeypatch.setattr(agent_loop, "adapter", model)
    alice = Client(base).login("alice@lab.org", "alice-password1")
    chat(alice, "P1", dict(message="plot ghost", provider="cloud", model="claude-opus-5-5"))
    result = model.seen[1]["messages"][-1]["results"][0]
    assert result["is_error"] and "not a file in this project" in result["text"]


HISTORY = [
    {"role": "user", "text": "plot"},
    {
        "role": "assistant",
        "text": "Surveying.",
        "calls": [{"id": "c1", "name": "survey_data", "input": {}}],
        "raw": {
            "provider": "cloud",
            "content": [
                {"type": "thinking", "thinking": "", "signature": "sig"},
                {"type": "text", "text": "Surveying."},
                {"type": "tool_use", "id": "c1", "name": "survey_data", "input": {}},
            ],
        },
    },
    {
        "role": "tool",
        "results": [
            {
                "id": "c1",
                "name": "survey_data",
                "text": "{}",
                "is_error": False,
                "image": {"media_type": "image/png", "data": "AAAA"},
            }
        ],
    },
]


def test_anthropic_wire_format():
    same = AnthropicAdapter(Provider(name="cloud", kind="anthropic", models=["m"])).wire(HISTORY)
    assert same[1]["content"][0]["type"] == "thinking"  # own blocks echoed unchanged
    result = same[2]["content"][0]
    assert result["type"] == "tool_result" and result["content"][1]["type"] == "image"
    other = AnthropicAdapter(
        Provider(name="gateway", kind="anthropic", models=["m"], vision=False)
    ).wire(HISTORY)
    assert [b["type"] for b in other[1]["content"]] == ["text", "tool_use"]  # rebuilt
    assert [b["type"] for b in other[2]["content"][0]["content"]] == ["text"]  # no vision
    marked = _echo_after_fallback(
        [{"type": "thinking"}, {"type": "text"}, {"type": "fallback"}, {"type": "tool_use"}]
    )
    assert [b["type"] for b in marked] == ["text", "fallback", "tool_use"]


class FakeOpenAI(BaseHTTPRequestHandler):
    requests: list = []

    def log_message(self, *args):
        pass

    def do_POST(self):  # noqa: N802
        FakeOpenAI.requests.append(
            dict(
                path=self.path,
                auth=self.headers.get("Authorization"),
                session=self.headers.get("x-opencode-session"),
                extra=self.headers.get("X-Team"),
                agent=self.headers.get("User-Agent"),
                body=json.loads(self.rfile.read(int(self.headers["Content-Length"]))),
            )
        )
        chunks = [
            {"choices": [{"delta": {"content": "Plotting "}}]},
            {
                "choices": [
                    {
                        "delta": {
                            "tool_calls": [
                                {
                                    "index": 0,
                                    "id": "call_9",
                                    "function": {
                                        "name": "run_plot",
                                        "arguments": '{"request": {"plot"',
                                    },
                                }
                            ]
                        }
                    }
                ]
            },
            {
                "choices": [
                    {
                        "delta": {
                            "tool_calls": [
                                {"index": 0, "function": {"arguments": ': "tracks/stacked"}}'}}
                            ]
                        },
                        "finish_reason": "tool_calls",
                    }
                ]
            },
            {"choices": [], "usage": {"prompt_tokens": 120, "completion_tokens": 30}, "model": "m"},
        ]
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.end_headers()
        for chunk in chunks:
            self.wfile.write(f"data: {json.dumps(chunk)}\n\n".encode())
        self.wfile.write(b"data: [DONE]\n\n")


def test_openai_compatible_adapter(monkeypatch):
    server = ThreadingHTTPServer(("127.0.0.1", 0), FakeOpenAI)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    monkeypatch.setenv("LAB_KEY", "secret")
    provider = Provider(
        name="lab",
        kind="openai",
        models=["m"],
        api_key_env="LAB_KEY",
        base_url=f"http://127.0.0.1:{server.server_address[1]}/v1",
    )
    streamed = []
    tools = [dict(name="run_plot", description="d", schema={"type": "object"})]
    reply = OpenAIAdapter(provider).complete("m", "sys", HISTORY, tools, streamed.append)
    server.shutdown()
    assert streamed == ["Plotting "] and reply.text == "Plotting "
    assert reply.calls == [
        dict(id="call_9", name="run_plot", input={"request": {"plot": "tracks/stacked"}})
    ]
    assert reply.stop == "tool_use" and reply.usage == dict(input=120, output=30)
    sent = FakeOpenAI.requests[-1]
    assert sent["path"] == "/v1/chat/completions" and sent["auth"] == "Bearer secret"
    body = sent["body"]
    assert body["max_tokens"] == 16000 and "max_completion_tokens" not in body
    roles = [m["role"] for m in body["messages"]]
    assert roles == ["system", "user", "assistant", "tool", "user"]  # figure follows as user image
    assert body["messages"][2]["tool_calls"][0]["function"]["name"] == "survey_data"
    assert body["messages"][4]["content"][1]["image_url"]["url"].startswith(
        "data:image/png;base64,"
    )


def test_unexpected_failure_keeps_conversation_usable(agent_studio, monkeypatch):
    app, base = agent_studio

    def crash(messages):
        raise AssertionError("stream ended without a final message")

    survey = Reply("", [dict(id="t1", name="survey_data", input={})], "tool_use")
    monkeypatch.setattr(agent_loop, "adapter", Scripted([lambda m: survey, crash]))
    alice = Client(base).login("alice@lab.org", "alice-password1")
    body = dict(message="plot", provider="cloud", model="claude-opus-5-5")
    events = chat(alice, "P1", body)
    assert [e["type"] for e in events][-2:] == ["error", "done"]
    assert "internal error" in events[-2]["message"]
    audit = app.db.one("SELECT detail FROM audit WHERE action = 'agent.turn'")["detail"]
    assert "internal error" in audit
    # The saved history is valid, so the next message works.
    monkeypatch.setattr(agent_loop, "adapter", Scripted([lambda m: Reply("ok", [], "end")]))
    events = chat(alice, "P1", dict(body, conversation=events[0]["id"]))
    assert events[-1]["type"] == "done" and "error" not in [e["type"] for e in events]


def test_opencode_provider_presets():
    settings = agent_settings(
        dict(
            providers=[
                dict(name="zen", kind="opencode", models=["kimi-k3"]),
                dict(name="go", kind="opencode", plan="go", models=["glm-5.3"]),
            ]
        )
    )
    zen, go = settings.providers
    assert zen.base_url == "https://opencode.ai/inference/openai/v1" and not zen.session_header
    assert zen.api_key_env == "OPENCODE_API_KEY" and not zen.official
    assert go.base_url == "https://opencode.ai/inference/go/openai/v1"
    assert go.session_header == "x-opencode-session"
    with pytest.raises(ValueError, match="only work inside OpenCode"):
        agent_settings(dict(providers=[dict(name="z", kind="opencode", models=["beam-free"])]))
    with pytest.raises(ValueError, match="plan"):
        agent_settings(dict(providers=[dict(name="z", kind="opencode", plan="pro", models=["m"])]))


def test_gateway_session_and_extra_headers(monkeypatch):
    server = ThreadingHTTPServer(("127.0.0.1", 0), FakeOpenAI)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    monkeypatch.setenv("OPENCODE_API_KEY", "oc-key")
    (provider,) = agent_settings(
        dict(
            providers=[
                dict(
                    name="go",
                    kind="opencode",
                    plan="go",
                    models=["glm-5.3"],
                    base_url=f"http://127.0.0.1:{server.server_address[1]}/v1",
                    headers={"X-Team": "core"},
                )
            ]
        )
    ).providers
    tools = [dict(name="run_plot", description="d", schema={"type": "object"})]
    OpenAIAdapter(provider).complete("glm-5.3", "sys", HISTORY[:1], tools, lambda t: None, "conv42")
    server.shutdown()
    sent = FakeOpenAI.requests[-1]
    assert sent["auth"] == "Bearer oc-key"
    assert sent["session"] == "conv42" and sent["extra"] == "core"
    assert sent["agent"] == "GenomeCanvas-Studio"  # not urllib's default, which gateways block


def test_anthropic_bearer_auth_for_gateways(monkeypatch):
    monkeypatch.setenv("GW_KEY", "gw-secret")
    gateway = Provider(
        name="oc-claude",
        kind="anthropic",
        models=["claude-opus-5-5"],
        base_url="https://opencode.ai/inference/anthropic",
        api_key_env="GW_KEY",
        auth="bearer",
    )
    client = AnthropicAdapter(gateway).client()
    assert client.auth_token == "gw-secret" and client.api_key is None
    assert str(client.base_url).startswith("https://opencode.ai/inference/anthropic")


def test_ca_bundle_choices(tmp_path):
    import ssl

    import certifi

    plain = Provider(name="a", kind="openai", models=["m"])
    assert ssl_context(plain) is None
    pem = Provider(name="b", kind="openai", models=["m"], ca_bundle=certifi.where())
    assert isinstance(ssl_context(pem), ssl.SSLContext)
    pytest.importorskip("truststore")
    system = Provider(name="c", kind="openai", models=["m"], ca_bundle="system")
    assert isinstance(ssl_context(system), ssl.SSLContext)
