"""Read-only progress and integrity inspection."""

from pathlib import Path

from .storage import read_json, receipt_valid


def status(path: str | Path) -> dict:
    """Report completed/missing units without reading molecular data into memory."""
    path = Path(path).resolve()
    if path.is_file():
        root, plan_path = path.parent, path
    else:
        root = path
        plan_path = (
            root / "plan.json"
            if (root / "plan.json").exists()
            else root / "build-plan.json"
        )
    if not plan_path.exists():
        return {
            "path": str(root),
            "state": "staging" if (root / "build.json").exists() else "missing",
            "completed_tasks": 0,
            "total_tasks": None,
        }
    plan = read_json(plan_path)
    build = plan["kind"] == "build"
    units = [
        task["shard"] if build else f"task-{task['id']:08d}" for task in plan["tasks"]
    ]
    completed = sum((root / "receipts" / f"{unit}.json").exists() for unit in units)
    published = (
        (root / "manifest.json").exists()
        if build
        else (root / "receipts/merge.json").exists()
    )
    return {
        "path": str(root),
        "kind": plan["kind"],
        "state": "complete"
        if published and completed == len(units)
        else "ready"
        if completed == len(units)
        else "incomplete",
        "completed_tasks": completed,
        "total_tasks": len(units),
    }


def verify(path: str | Path) -> dict:
    """Checksum all staged/catalog/shard/task outputs and reject changed plans.

    Returns ``ok`` plus individual failed receipt names. Verification is read-only.
    A complete status is not an integrity guarantee; verification reads every
    output byte and therefore costs sequential I/O proportional to index size.
    """
    path = Path(path).resolve()
    root = path.parent if path.is_file() else path
    failures = []
    count = 0
    for receipt in sorted((root / "receipts").glob("*.json")):
        count += 1
        try:
            expected = read_json(receipt)["signature"]
            if not receipt_valid(root, receipt.stem, expected):
                failures.append(receipt.stem)
        except (ValueError, KeyError, OSError):
            failures.append(receipt.stem)
    try:
        if (root / "plan.json").exists():
            from .operations import _load_plan, _task_signature

            _, plan, _ = _load_plan(root / "plan.json")
            if not receipt_valid(root, "queries", plan["signature"]):
                failures.append("queries")
            for task in plan["tasks"]:
                if task.get("query_chunk") and not receipt_valid(
                    root, task["query_chunk"], task["query_chunk_signature"]
                ):
                    failures.append(task["query_chunk"])
                unit = f"task-{task['id']:08d}"
                if not receipt_valid(root, unit, _task_signature(plan, task)):
                    failures.append(unit)
            if not receipt_valid(root, "merge", plan["signature"]):
                failures.append("merge")
        elif (root / "build-plan.json").exists():
            from .index import _build_plan

            _, plan = _build_plan(root / "build-plan.json")
            for staged in plan["manifest"]["stage_units"]:
                if not receipt_valid(root, staged["id"], staged["signature"]):
                    failures.append(staged["id"])
            for shard in plan["manifest"]["shards"]:
                if not receipt_valid(root, shard["id"], shard["signature"]):
                    failures.append(shard["id"])
            if (
                not (root / "manifest.json").exists()
                or read_json(root / "manifest.json") != plan["manifest"]
            ):
                failures.append("manifest")
        else:
            failures.append("plan")
    except (ValueError, KeyError, OSError) as error:
        failures.append(f"plan: {error}")
    return {
        "ok": not failures,
        "checked_receipts": count,
        "failures": sorted(set(failures)),
        **status(path),
    }
