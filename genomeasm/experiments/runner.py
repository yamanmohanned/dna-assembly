"""Running experiment tasks in parallel, resumably.

A *task* is one read set (instance spec + seed index) and the solver runs to make on it. Tasks
run in separate processes (``spawn``), each with one numba thread, five at a time. Results are
appended to a JSONL work file as they arrive, so an interrupted experiment resumes where it stopped;
at the end the work file is compressed into ``results/``.
"""

from __future__ import annotations

import gzip
import json
import multiprocessing as mp
import os
import shutil
import time
import traceback
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np


def _json_default(value: object) -> object:
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def init_worker() -> None:
    """Compile (or load from cache) all kernels once per worker, before anything is timed."""
    from genomeasm.solvers.warmup import warmup

    warmup()


def run_task(task: dict[str, Any]) -> list[dict[str, Any]]:
    """Run every solver of a task on its read set; errors are recorded, not raised."""
    from genomeasm.experiments import design

    cfg, spec = task["config"], design.InstanceSpec.from_config(task["config"], task["spec"])
    head = {"stage": task["stage"], "instance": spec.id, "seed_index": task["seed_index"]}
    records: list[dict[str, Any]] = []
    try:
        built = design.build(spec, task["read_seed"])
        inst = built.instance
        order = inst.true_order()
        if task.get("emit_instance", True):
            records.append(
                {
                    **head,
                    "kind": "instance",
                    "n": inst.n,
                    "n_input": inst.n_input,
                    "edges": int((inst.ov > 0).sum()),
                    "build_seconds": round(built.build_seconds, 3),
                    "F_true": int(inst.fitness(order)),
                    "bio_true": design.bio_metrics(built, order),
                }
            )
    except Exception:  # noqa: BLE001
        return [{**head, "kind": "error", "error": traceback.format_exc()[-600:]}]
    for run in task["runs"]:
        try:
            result = design.run_algorithm(
                run["algorithm"], built, cfg, run["params"], task["seconds"], task["solver_seed"],
                task.get("lkh_search_seconds"),
            )  # fmt: skip
            record = design.evaluate(built, result, with_bio=task.get("bio", True))
            record.update(head)
            record.update(kind="run", requested=run["algorithm"], trial=run.get("trial"))
            record["n"] = inst.n
            records.append(record)
        except Exception:  # noqa: BLE001
            records.append(
                {**head, "kind": "error", "algorithm": run["algorithm"], "trial": run.get("trial"),
                 "error": traceback.format_exc()[-600:]}
            )  # fmt: skip
    return records


def task_key(task: dict[str, Any]) -> str:
    return f"{task['stage']}|{task['spec']['id']}|{task['seed_index']}|{task['group']}"


def execute(
    tasks: list[dict[str, Any]],
    work_file: Path,
    workers: int,
    say: Callable[[str], None] = print,
) -> None:
    """Run the tasks that are not yet in ``work_file`` and append their records to it."""
    work_file.parent.mkdir(parents=True, exist_ok=True)
    done: set[str] = set()
    if work_file.is_file():
        for line in work_file.read_text(encoding="utf-8").splitlines():
            rec = json.loads(line)
            if rec.get("kind") == "task_done":
                done.add(rec["key"])
    pending = [t for t in tasks if task_key(t) not in done]
    say(f"{len(tasks)} tasks, {len(done)} already done, {len(pending)} to run")
    if not pending:
        return
    os.environ["NUMBA_NUM_THREADS"] = "1"
    started = time.time()
    context = mp.get_context("spawn")
    with (
        context.Pool(workers, initializer=init_worker) as pool,
        open(work_file, "a", encoding="utf-8", newline="\n") as out,
    ):
        for count, (task, records) in enumerate(zip_unordered(pool, pending), start=1):
            for record in records:
                out.write(json.dumps(record, default=_json_default) + "\n")
            out.write(json.dumps({"kind": "task_done", "key": task_key(task)}) + "\n")
            out.flush()
            if count % max(1, len(pending) // 20) == 0 or count == len(pending):
                elapsed = time.time() - started
                left = elapsed / count * (len(pending) - count) / 60
                say(f"  {count}/{len(pending)} tasks, {elapsed / 60:.1f} min, ~{left:.0f} min left")


def zip_unordered(pool: Any, tasks: list[dict[str, Any]]) -> Any:
    """Yield (task, records) as tasks finish."""
    for index, records in pool.imap_unordered(_indexed, list(enumerate(_payloads(tasks)))):
        yield tasks[index], records


def _payloads(tasks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return tasks


def _indexed(item: tuple[int, dict[str, Any]]) -> tuple[int, list[dict[str, Any]]]:
    index, task = item
    return index, run_task(task)


def read_records(path: Path) -> list[dict[str, Any]]:
    """Read a work file or its compressed copy (without the bookkeeping lines)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        rows = [json.loads(line) for line in handle if line.strip()]
    return [r for r in rows if r.get("kind") != "task_done"]


def compress(work_file: Path, target: Path) -> None:
    """Keep the raw records of a finished experiment as a .jsonl.gz under results/."""
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(work_file, "rb") as src, gzip.GzipFile(target, "wb", mtime=0) as dst:
        shutil.copyfileobj(src, dst)
