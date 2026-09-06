"""Freeze or verify explicit study inputs without assigning evaluation independence.

Example:
  python scripts/freeze_study.py freeze --root . --plan study.json --output frozen.json
  python scripts/freeze_study.py verify --root . --manifest frozen.json
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gembench.study import freeze_study, verify_study  # noqa: E402


def _read_json(path):
    def unique_object(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key!r}")
            result[key] = value
        return result

    def reject_constant(value):
        raise ValueError(f"Nonfinite JSON value: {value}")

    with Path(path).open(encoding="utf-8") as fh:
        return json.load(fh, object_pairs_hook=unique_object, parse_constant=reject_constant)


def _write_new_json(path, value):
    """Publish a complete manifest without overwriting an existing freeze."""
    path = Path(path)
    if not path.parent.is_dir():
        raise ValueError("The output directory must already exist")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix=".study-freeze-", suffix=".tmp", delete=False) as fh:
            temporary = Path(fh.name)
            json.dump(value, fh, indent=2, ensure_ascii=False, allow_nan=False)
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        # An atomic hard link fails if another writer already created the path.
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    freeze = modes.add_parser("freeze", help="Create a new manifest; existing manifests are never overwritten")
    freeze.add_argument("--root", required=True, help="Root directory against which declared relative paths are resolved")
    freeze.add_argument("--plan", required=True, help="JSON study declaration")
    freeze.add_argument("--output", required=True, help="New manifest JSON file")
    verify = modes.add_parser("verify", help="Check manifest integrity and frozen file bytes")
    verify.add_argument("--root", required=True)
    verify.add_argument("--manifest", required=True)
    args = parser.parse_args(argv)
    try:
        if args.mode == "freeze":
            manifest = freeze_study(args.root, _read_json(args.plan))
            _write_new_json(args.output, manifest)
            print(json.dumps({"manifest": str(Path(args.output)), "study_id": manifest["plan"]["study_id"],
                              "n_files": len(manifest["files"]), "content_fingerprint": manifest["content_fingerprint"],
                              "evaluation_caveat": manifest["evaluation_caveat"]}, indent=2))
            return 0
        report = verify_study(args.root, _read_json(args.manifest))
        print(json.dumps(report, indent=2, allow_nan=False))
        return 0 if report["valid"] else 1
    except (OSError, ValueError, TypeError) as error:
        print(json.dumps({"valid": False, "errors": [str(error)]}, indent=2), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
