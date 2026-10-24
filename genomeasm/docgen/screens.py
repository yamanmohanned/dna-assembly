"""Screenshots of the application pages: the fallback if the live demo fails (``genomeasm docs screenshots``).

A server is started on its own port, headless Chrome takes one picture per page (the default state of
the page, with the settings the demo starts from), and the server is stopped again.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

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


def _wait_until_up(port: int, seconds: float = 60.0) -> None:
    deadline = time.time() + seconds
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(
                f"http://localhost:{port}/_stcore/health", timeout=2
            ) as reply:
                if reply.read().strip() == b"ok":
                    return
        except OSError:
            time.sleep(1)
    raise RuntimeError("the application did not start")


def take_screenshots(out_dir: Path, port: int = 8599) -> list[Path]:
    chrome = md.chrome_path()
    if chrome is None:
        raise RuntimeError("Chrome is needed for screenshots")
    out_dir.mkdir(parents=True, exist_ok=True)
    server = subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", str(ROOT / "app" / "streamlit_app.py"),
         "--server.port", str(port), "--server.headless", "true", "--browser.gatherUsageStats", "false"],
        cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )  # fmt: skip
    written: list[Path] = []
    try:
        _wait_until_up(port)
        for name, path in PAGES:
            target = (out_dir / f"{name}.png").resolve()
            with tempfile.TemporaryDirectory() as profile:
                subprocess.run(
                    [str(chrome), "--headless=new", "--disable-gpu", "--no-sandbox", "--hide-scrollbars",
                     f"--user-data-dir={profile}", "--window-size=1440,900", "--virtual-time-budget=25000",
                     f"--screenshot={target}", f"http://localhost:{port}/{path}"],
                    check=False, capture_output=True, timeout=180,
                )  # fmt: skip
            if target.is_file() and target.stat().st_size > 5000:
                written.append(target)
    finally:
        server.terminate()
        try:
            server.wait(timeout=15)
        except subprocess.TimeoutExpired:
            server.kill()
    return written
