"""Screenshots of the application pages: the fallback if the live demo fails (``genomeasm docs screenshots``).

A server is started on its own port and a headless Chrome is driven through the DevTools protocol:
for every page the picture is taken after a real waiting time, because Streamlit draws a page while
it computes it. Two passes: the first warms every page up, the second is kept.
"""

from __future__ import annotations

import asyncio
import base64
import json
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

from websockets.asyncio.client import connect

from genomeasm.appkit.paths import ROOT
from genomeasm.docgen import md

PAGES = [
    ("01_guide", ""),
    ("02_learn", "learn"),
    ("03_data", "data"),
    ("04_simulator", "simulator"),
    ("05_overlap_graph", "overlap_graph"),
    ("06_solver_lab", "solver_lab"),
    ("07_experiments", "experiments"),
    ("08_assembly_viewer", "assembly_viewer"),
    ("09_reproducibility", "reproducibility"),
]
WIDTH, HEIGHT = 1440, 900


def _wait_until_up(url: str, expect: bytes | None, seconds: float = 60.0) -> bytes:
    deadline = time.time() + seconds
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(url, timeout=2) as reply:
                body: bytes = reply.read().strip()
                if expect is None or body == expect:
                    return body
        except OSError:
            time.sleep(1)
    raise RuntimeError(f"{url} did not answer")


async def _capture(ws_url: str, delay: float) -> bytes:
    async with connect(ws_url, max_size=None) as socket:
        counter = 0

        async def call(method: str, params: dict[str, object] | None = None) -> dict[str, object]:
            nonlocal counter
            counter += 1
            await socket.send(json.dumps({"id": counter, "method": method, "params": params or {}}))
            while True:
                message = json.loads(await socket.recv())
                if message.get("id") == counter:
                    return dict(message)

        size: dict[str, object] = {
            "width": WIDTH,
            "height": HEIGHT,
            "deviceScaleFactor": 1,
            "mobile": False,
        }
        await call("Emulation.setDeviceMetricsOverride", size)
        await asyncio.sleep(delay)
        reply = await call("Page.captureScreenshot", {"format": "png"})
    result = dict(reply["result"])  # type: ignore[call-overload]
    return base64.b64decode(str(result["data"]))


def _open_page(port: int, url: str) -> dict[str, str]:
    request = urllib.request.Request(f"http://127.0.0.1:{port}/json/new?{url}", method="PUT")
    with urllib.request.urlopen(request, timeout=20) as reply:
        return dict(json.loads(reply.read()))


def take_screenshots(out_dir: Path, app_port: int = 8599, debug_port: int = 9333) -> list[Path]:
    chrome = md.chrome_path()
    if chrome is None:
        raise RuntimeError("Chrome is needed for screenshots")
    out_dir.mkdir(parents=True, exist_ok=True)
    server = subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", str(ROOT / "app" / "streamlit_app.py"),
         "--server.port", str(app_port), "--server.headless", "true", "--browser.gatherUsageStats", "false"],
        cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )  # fmt: skip
    written: list[Path] = []
    with tempfile.TemporaryDirectory() as profile:
        browser = subprocess.Popen(
            [str(chrome), "--headless=new", "--disable-gpu", "--no-sandbox", f"--user-data-dir={profile}",
             f"--remote-debugging-port={debug_port}", "about:blank"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )  # fmt: skip
        try:
            _wait_until_up(f"http://localhost:{app_port}/_stcore/health", b"ok")
            _wait_until_up(f"http://127.0.0.1:{debug_port}/json/version", None)
            for keep in (False, True):
                for name, path in PAGES:
                    target = _open_page(debug_port, f"http://localhost:{app_port}/{path}")
                    try:
                        data = asyncio.run(
                            _capture(target["webSocketDebuggerUrl"], 14.0 if keep else 10.0)
                        )
                    finally:
                        close = f"http://127.0.0.1:{debug_port}/json/close/{target['id']}"
                        urllib.request.urlopen(close, timeout=20).read()
                    if keep and len(data) > 20_000:
                        file = out_dir / f"{name}.png"
                        file.write_bytes(data)
                        written.append(file)
        finally:
            browser.terminate()
            server.terminate()
            for process in (browser, server):
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    process.kill()
    return written
