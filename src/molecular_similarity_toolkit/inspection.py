"""Progress and the shared integrity contract for indexes and scientific runs."""

from pathlib import Path

from .storage import check_units, read_json


def run_dependencies(root, plan, reference, *, include_merge=True):
    """Validate every required scientific unit, independent of task cardinality."""
    check_units(
        root,
        plan["required_units"],
        skip=() if include_merge else ("merge",),
        exact=True,
    )
    from .index import MoleculeIndex

    seen = set()
    for index in (reference, MoleculeIndex.open(plan["queries"])):
        if index.path not in seen:
            check_units(index.path, index.manifest["required_units"])
            seen.add(index.path)


def status(path: str | Path) -> dict:
    """Report staged chunks/selected rows and planned task progress without hashing."""
    path = Path(path).resolve()
    root = path.parent if path.is_file() else path
    try:
        if (root / "manifest.json").exists() and read_json(root / "manifest.json").get(
            "compact"
        ):
            from .index import MoleculeIndex

            index = MoleculeIndex.open(root)
            return {
                "path": str(root),
                "kind": "index",
                "state": "complete",
                "compact": True,
                "completed_tasks": len(index.manifest["shards"]),
                "total_tasks": len(index.manifest["shards"]),
            }
        plan_path = root / (
            "plan.json" if (root / "plan.json").exists() else "build-plan.json"
        )
        if not plan_path.exists():
            progress = (
                read_json(root / "progress.json")
                if (root / "progress.json").exists()
                else {}
            )
            return {
                "path": str(root),
                "state": "staging" if (root / "build.json").exists() else "missing",
                "completed_tasks": 0,
                "total_tasks": None,
                **progress,
            }
        plan = read_json(plan_path)
        build = plan["kind"] == "build"
        names = [t["shard"] if build else f"task-{t['id']:08d}" for t in plan["tasks"]]
        completed = sum((root / "receipts" / f"{name}.json").exists() for name in names)
        published = (
            (root / "manifest.json").exists()
            if build
            else (root / "receipts/merge.json").exists()
        )
        progress = (
            read_json(root / "progress.json")
            if build and (root / "progress.json").exists()
            else {}
        )
        return {
            **progress,
            "path": str(root),
            "kind": plan["kind"],
            "state": "complete"
            if published and completed == len(names)
            else "ready"
            if completed == len(names)
            else "incomplete",
            "completed_tasks": completed,
            "total_tasks": len(names),
        }
    except (ValueError, KeyError, TypeError, OSError) as error:
        return {
            "path": str(root),
            "state": "corrupt",
            "error": str(error),
            "completed_tasks": None,
            "total_tasks": None,
        }


def verify(path: str | Path) -> dict:
    """Read/checksum required artifacts and their scientific inputs; never mutate.

    Includes prepared queries even with zero pool tasks. Full indexes include
    build dependencies; compact indexes need only runtime catalog/shard units.
    Returns ``ok`` and failures; the CLI exits nonzero on any failure.
    """
    path = Path(path).resolve()
    root = path.parent if path.is_file() else path
    failures = []
    try:
        if (root / "plan.json").exists():
            from .operations import _load_plan

            _, plan, reference = _load_plan(root / "plan.json")
            run_dependencies(root, plan, reference)
        elif (root / "manifest.json").exists() and read_json(
            root / "manifest.json"
        ).get("compact"):
            from .index import MoleculeIndex

            index = MoleculeIndex.open(root)
            check_units(root, index.manifest["required_units"], exact=True)
        elif (root / "build-plan.json").exists():
            from .index import _build_plan

            _, plan = _build_plan(root / "build-plan.json")
            check_units(root, plan["required_units"], exact=True)
            if (
                not (root / "manifest.json").exists()
                or read_json(root / "manifest.json") != plan["manifest"]
            ):
                raise ValueError("manifest missing or changed")
        else:
            raise ValueError("complete index or operation plan is missing")
    except (ValueError, KeyError, TypeError, OSError) as error:
        failures.append(str(error))
    return {
        "ok": not failures,
        "failures": failures,
        "checked_receipts": len(list((root / "receipts").glob("*.json"))),
        **status(path),
    }
