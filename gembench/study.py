"""Freeze explicitly named study inputs and verify their exact bytes.

This records a declaration and its files. It does not establish that a dataset
was uninspected, prevent outcome-driven curation, or authenticate a manifest.
Only named files are frozen; repository provenance is a separate snapshot.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import re
import stat
import subprocess


SCHEMA_VERSION = 1
MANIFEST_TYPE = "gembench-study-freeze"
DEVELOPMENT_ORGANISMS = ("Btheta", "Putida", "MR1", "Smeli", "Keio")
EVALUATION_CAVEAT = (
    "File hashes record declared inputs, not evaluation independence. The five named "
    "development organisms are exposed. Other downloaded organisms have unknown prior "
    "exposure and remain quarantined pending independent evaluation. Documentary claims "
    "about prospective external evaluation require independent review."
)
_RESERVED_PLAN_KEYS = {
    "version", "manifest_type", "files", "content_fingerprint", "repository", "created_at",
    "evaluation_caveat", "held_out_organisms", "heldout_organisms", "holdout_organisms",
}


def _json_bytes(value):
    def check(item):
        if isinstance(item, dict):
            if any(not isinstance(key, str) for key in item):
                raise ValueError("JSON object keys must be strings")
            for child in item.values():
                check(child)
        elif isinstance(item, list):
            for child in item:
                check(child)
        elif item is None or isinstance(item, (str, bool, int)):
            pass
        elif isinstance(item, float) and math.isfinite(item):
            pass
        else:
            raise ValueError("Study metadata must contain only finite JSON values")

    check(value)
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                          allow_nan=False).encode("utf-8")
    except (UnicodeError, ValueError) as error:
        raise ValueError("Study metadata must be valid UTF-8 JSON") from error


def _relative_path(value):
    if not isinstance(value, str) or not value or "\0" in value or "\\" in value:
        raise ValueError("Input paths must be nonempty relative POSIX paths")
    path = PurePosixPath(value)
    if path.is_absolute() or PureWindowsPath(value).drive or ".." in path.parts:
        raise ValueError(f"Unsafe path outside the study root: {value!r}")
    if not path.parts or path.parts[0] in {".git", ".venv"}:
        raise ValueError(f"Reserved or invalid input path: {value!r}")
    return path.as_posix()


def _names(value, field):
    if not isinstance(value, list) or any(not isinstance(v, str) or not v.strip() for v in value):
        raise ValueError(f"{field} must be a list of nonempty organism identifiers")
    if len(set(value)) != len(value):
        raise ValueError(f"Duplicate organism identifiers in {field}")
    return value


def _normalize_plan(plan):
    if not isinstance(plan, dict):
        raise ValueError("The study plan must be a JSON object")
    normalized = json.loads(_json_bytes(plan))
    if type(normalized.get("schema_version")) is not int or normalized["schema_version"] != SCHEMA_VERSION:
        raise ValueError(f"Unsupported study plan schema_version; expected {SCHEMA_VERSION}")
    reserved = set(normalized) & _RESERVED_PLAN_KEYS
    if reserved:
        raise ValueError(f"Reserved study plan keys: {sorted(reserved)}")
    study_id = normalized.get("study_id")
    if not isinstance(study_id, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", study_id):
        raise ValueError("study_id must be a 1-128 character identifier beginning with a letter or digit")
    role = normalized.get("evaluation_role")
    if role not in {"development", "prospective_external"}:
        raise ValueError("evaluation_role must be development or documented prospective_external")
    development = _names(normalized.get("development_organisms"), "development_organisms")
    if set(development) != set(DEVELOPMENT_ORGANISMS):
        raise ValueError(f"development_organisms must declare the known exposed set: {list(DEVELOPMENT_ORGANISMS)}")
    quarantined = _names(normalized.get("quarantined_organisms"), "quarantined_organisms")
    if set(development) & set(quarantined):
        raise ValueError("Known development organisms cannot also be quarantined")
    raw_paths = normalized.get("paths")
    if not isinstance(raw_paths, list) or not raw_paths:
        raise ValueError("paths must be a nonempty list of explicitly named files")
    paths = [_relative_path(value) for value in raw_paths]
    if len(set(paths)) != len(paths):
        raise ValueError("Duplicate input paths after normalization")
    normalized["paths"] = sorted(paths)
    if role == "prospective_external":
        exposure = normalized.get("candidate_exposure")
        if not isinstance(exposure, dict):
            raise ValueError("prospective_external requires documented candidate_exposure")
        for field in ("external_custodian", "uninspected_evidence"):
            if not isinstance(exposure.get(field), str) or not exposure[field].strip():
                raise ValueError(f"candidate_exposure.{field} must document the prospective external evaluation")
        evidence = exposure.get("evidence_paths")
        if not isinstance(evidence, list) or not evidence:
            raise ValueError("candidate_exposure.evidence_paths must name frozen documentation files")
        evidence = [_relative_path(path) for path in evidence]
        if len(set(evidence)) != len(evidence) or not set(evidence) <= set(paths):
            raise ValueError("Exposure evidence paths must be unique and included in paths")
        exposure["evidence_paths"] = sorted(evidence)
    return normalized


def _root_path(root):
    resolved = Path(root).resolve(strict=True)
    if not resolved.is_dir():
        raise ValueError("The study root must be a directory")
    return resolved


def _resolve_input(root, relative):
    candidate = root / _relative_path(relative)
    try:
        resolved = candidate.resolve(strict=False)
    except (OSError, RuntimeError) as error:
        raise ValueError(f"Cannot safely resolve input path {relative!r}") from error
    if not resolved.is_relative_to(root):
        raise ValueError(f"Input symlink escapes the study root: {relative!r}")
    _relative_path(resolved.relative_to(root).as_posix())
    info = resolved.stat()  # preserves FileNotFoundError for verification reports
    if not stat.S_ISREG(info.st_mode):
        raise ValueError(f"Study inputs must be regular files: {relative!r}")
    return resolved


def _hash_file(root, relative):
    resolved = _resolve_input(root, relative)
    digest = hashlib.sha256()
    size = 0
    with resolved.open("rb") as fh:
        before = os.fstat(fh.fileno())
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
        after = os.fstat(fh.fileno())
    identity = lambda info: (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns)
    if (identity(before) != identity(after) or identity(after) != identity(resolved.stat()) or
            size != after.st_size or _resolve_input(root, relative) != resolved):
        raise ValueError(f"Input changed while being frozen or verified: {relative!r}")
    return {"path": relative, "size_bytes": size, "sha256": digest.hexdigest()}


def _git(root, *args):
    try:
        return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True,
                              encoding="utf-8", timeout=30, check=False)
    except (OSError, subprocess.TimeoutExpired, UnicodeError) as error:
        raise ValueError("Could not read repository provenance") from error


def _repository_snapshot(root):
    top = _git(root, "rev-parse", "--show-toplevel")
    empty = {"available": False, "head": None, "dirty": None, "tracked_changes": [],
             "untracked_paths": [], "excluded_untracked_roots": [".venv"]}
    if top.returncode:
        return empty
    head = _git(root, "rev-parse", "--verify", "HEAD")
    status = _git(root, "status", "--porcelain=v1", "-z", "--untracked-files=all", "--", ".",
                  ":(exclude).venv", ":(exclude).venv/**")
    if status.returncode:
        raise ValueError("Could not read Git working-tree status")
    entries = status.stdout.split("\0")
    changes, untracked = [], []
    i = 0
    while i < len(entries):
        entry = entries[i]
        i += 1
        if not entry:
            continue
        code, path = entry[:2], entry[3:]
        if code == "??":
            untracked.append(path)
        else:
            change = {"path": path, "status": code}
            if "R" in code or "C" in code:
                if i >= len(entries) or not entries[i]:
                    raise ValueError("Incomplete Git rename record")
                change["original_path"] = entries[i]
                i += 1
            changes.append(change)
    return {"available": True, "head": head.stdout.strip() if head.returncode == 0 else None,
            "dirty": bool(changes or untracked), "tracked_changes": sorted(changes, key=lambda item: item["path"]),
            "untracked_paths": sorted(untracked), "excluded_untracked_roots": [".venv"]}


def _fingerprint(plan, files):
    content = {"schema_version": SCHEMA_VERSION, "manifest_type": MANIFEST_TYPE,
               "plan": plan, "files": sorted(files, key=lambda item: item["path"])}
    return hashlib.sha256(_json_bytes(content)).hexdigest()


def _timestamp(created_at):
    if created_at is None:
        return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    if not isinstance(created_at, str) or not created_at.endswith("Z"):
        raise ValueError("created_at must be a UTC ISO timestamp ending in Z")
    try:
        datetime.fromisoformat(created_at[:-1] + "+00:00")
    except ValueError as error:
        raise ValueError("Invalid created_at timestamp") from error
    return created_at


def freeze_study(root, plan, *, created_at=None):
    """Return a manifest for named files; never discover or select evaluation data.

    The content fingerprint binds the normalized plan and exact file records.
    Timestamp and Git provenance are separate, so copied checkouts with identical
    named inputs have the same content fingerprint. They are not signatures.
    """
    root = _root_path(root)
    plan = _normalize_plan(plan)
    resolved = [_resolve_input(root, path) for path in plan["paths"]]
    if len(set(resolved)) != len(resolved):
        raise ValueError("Duplicate input paths resolve to the same file")
    files = [_hash_file(root, path) for path in plan["paths"]]
    return {"schema_version": SCHEMA_VERSION, "manifest_type": MANIFEST_TYPE, "created_at": _timestamp(created_at),
            "plan": plan, "files": files, "repository": _repository_snapshot(root),
            "content_fingerprint": _fingerprint(plan, files), "evaluation_caveat": EVALUATION_CAVEAT}


def _validate_manifest(manifest):
    required = {"schema_version", "manifest_type", "created_at", "plan", "files", "repository",
                "content_fingerprint", "evaluation_caveat"}
    if not isinstance(manifest, dict) or set(manifest) != required:
        raise ValueError("Manifest fields are missing, reserved, or unsupported")
    _json_bytes(manifest)
    if type(manifest["schema_version"]) is not int or manifest["schema_version"] != SCHEMA_VERSION:
        raise ValueError("Unsupported manifest schema_version")
    if manifest["manifest_type"] != MANIFEST_TYPE or manifest["evaluation_caveat"] != EVALUATION_CAVEAT:
        raise ValueError("Unsupported manifest type or evaluation caveat")
    if not isinstance(manifest["created_at"], str):
        raise ValueError("Manifest created_at must be a UTC timestamp")
    _timestamp(manifest["created_at"])
    plan = _normalize_plan(manifest["plan"])
    files = manifest["files"]
    if not isinstance(files, list) or not files:
        raise ValueError("Manifest files must be a nonempty list")
    paths = []
    for record in files:
        if not isinstance(record, dict) or set(record) != {"path", "size_bytes", "sha256"}:
            raise ValueError("Invalid or incomplete file record")
        path = _relative_path(record["path"])
        if path != record["path"]:
            raise ValueError("Manifest file paths must be normalized")
        paths.append(path)
        if type(record["size_bytes"]) is not int or record["size_bytes"] < 0:
            raise ValueError("File size_bytes must be a nonnegative integer")
        if not isinstance(record["sha256"], str) or not re.fullmatch(r"[0-9a-f]{64}", record["sha256"]):
            raise ValueError("File sha256 must contain exactly 64 lowercase hexadecimal characters")
    if len(set(paths)) != len(paths) or sorted(paths) != plan["paths"]:
        raise ValueError("Manifest must have exactly one file record for every declared path")
    repository = manifest["repository"]
    repo_fields = {"available", "head", "dirty", "tracked_changes", "untracked_paths", "excluded_untracked_roots"}
    if not isinstance(repository, dict) or set(repository) != repo_fields or type(repository["available"]) is not bool:
        raise ValueError("Invalid repository provenance")
    head = repository["head"]
    if head is not None and (not isinstance(head, str) or not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", head)):
        raise ValueError("Invalid repository commit identifier")
    if repository["available"] and type(repository["dirty"]) is not bool:
        raise ValueError("Repository dirty status must be boolean")
    if not repository["available"] and (head is not None or repository["dirty"] is not None):
        raise ValueError("Unavailable repository provenance cannot declare a commit or dirty status")
    if not isinstance(repository["tracked_changes"], list) or not isinstance(repository["untracked_paths"], list):
        raise ValueError("Repository path records must be lists")
    tracked_paths = []
    for change in repository["tracked_changes"]:
        if (not isinstance(change, dict) or set(change) not in
                ({"path", "status"}, {"path", "status", "original_path"})):
            raise ValueError("Invalid tracked-change provenance record")
        tracked_paths.append(_relative_path(change["path"]))
        if "original_path" in change:
            _relative_path(change["original_path"])
        if not isinstance(change["status"], str) or not re.fullmatch(r"[ MADRCUT?!]{2}", change["status"]):
            raise ValueError("Invalid Git change status")
    untracked_paths = [_relative_path(path) for path in repository["untracked_paths"]]
    if len(set(tracked_paths)) != len(tracked_paths) or len(set(untracked_paths)) != len(untracked_paths):
        raise ValueError("Duplicate repository provenance paths")
    if not repository["available"] and (tracked_paths or untracked_paths):
        raise ValueError("Unavailable repository provenance cannot declare changed files")
    if repository["available"] and repository["dirty"] != bool(tracked_paths or untracked_paths):
        raise ValueError("Repository dirty status disagrees with its path records")
    if repository["excluded_untracked_roots"] != [".venv"]:
        raise ValueError("Unsupported repository provenance exclusions")
    claimed = manifest["content_fingerprint"]
    if not isinstance(claimed, str) or not re.fullmatch(r"[0-9a-f]{64}", claimed):
        raise ValueError("Invalid content_fingerprint")
    if claimed != _fingerprint(plan, files):
        raise ValueError("Manifest content fingerprint does not match the declared plan and file records")
    return plan, files


def verify_study(root, manifest):
    """Return valid/changed/missing/errors; a malformed manifest never passes."""
    report = {"valid": False, "checked_files": 0, "missing": [], "changed": [], "errors": []}
    try:
        plan, records = _validate_manifest(manifest)
        root = _root_path(root)
    except (ValueError, OSError, TypeError) as error:
        report["errors"].append(str(error))
        return report
    report["content_fingerprint"] = manifest["content_fingerprint"]
    resolved_paths = set()
    for expected in records:
        path = expected["path"]
        try:
            resolved = _resolve_input(root, path)
            if resolved in resolved_paths:
                raise ValueError(f"Duplicate input paths resolve to the same file: {path!r}")
            resolved_paths.add(resolved)
            actual = _hash_file(root, path)
        except FileNotFoundError:
            report["missing"].append(path)
            continue
        except (ValueError, OSError) as error:
            report["errors"].append(f"{path}: {error}")
            continue
        report["checked_files"] += 1
        if actual != expected:
            report["changed"].append({"path": path, "expected_sha256": expected["sha256"],
                                      "actual_sha256": actual["sha256"], "expected_size_bytes": expected["size_bytes"],
                                      "actual_size_bytes": actual["size_bytes"]})
    report["valid"] = not (report["missing"] or report["changed"] or report["errors"])
    return report
