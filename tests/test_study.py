from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from gembench.study import DEVELOPMENT_ORGANISMS, freeze_study, verify_study
from scripts.freeze_study import main


def plan(**updates):
    declaration = {
        "schema_version": 1,
        "study_id": "development-2026-09",
        "evaluation_role": "development",
        "development_organisms": list(DEVELOPMENT_ORGANISMS),
        "quarantined_organisms": ["uninspected-candidate"],
        "paths": ["model.json", "method.py"],
        "question": "Does a predeclared correction improve the development benchmark?",
        "decision_rule": {"metric": "balanced_accuracy", "threshold": 0.02},
        "caveats": ["Development performance is not independent evaluation."],
    }
    declaration.update(updates)
    return declaration


def populate(root):
    root.mkdir(parents=True, exist_ok=True)
    (root / "model.json").write_bytes(b'{"model":1}\r\n')
    (root / "method.py").write_bytes(b"method = 1\n")
    return root


def test_freeze_exact_bytes_metadata_and_relocated_checkout(tmp_path):
    first_root = populate(tmp_path / "first")
    second_root = populate(tmp_path / "different-name")
    declaration = plan()
    first = freeze_study(first_root, declaration, created_at="2026-09-06T12:00:00Z")
    second = freeze_study(second_root, declaration, created_at="2026-09-07T12:00:00Z")
    assert first["content_fingerprint"] == second["content_fingerprint"]
    assert first["created_at"] != second["created_at"]
    assert first["plan"]["decision_rule"] == declaration["decision_rule"]
    assert declaration["paths"] == ["model.json", "method.py"]  # caller data unchanged
    assert first["plan"]["paths"] == ["method.py", "model.json"]
    assert first["files"][1] == {
        "path": "model.json", "size_bytes": 13,
        "sha256": hashlib.sha256(b'{"model":1}\r\n').hexdigest(),
    }
    assert str(first_root) not in json.dumps(first)
    assert verify_study(second_root, first)["valid"]
    assert "not evaluation independence" in first["evaluation_caveat"]


def test_verification_reports_same_size_changes_and_missing_inputs(tmp_path):
    populate(tmp_path)
    manifest = freeze_study(tmp_path, plan())
    (tmp_path / "method.py").write_bytes(b"method = 2\n")
    (tmp_path / "model.json").unlink()
    report = verify_study(tmp_path, manifest)
    assert not report["valid"]
    assert report["checked_files"] == 1
    assert report["missing"] == ["model.json"]
    assert report["changed"][0]["path"] == "method.py"
    assert report["changed"][0]["expected_size_bytes"] == report["changed"][0]["actual_size_bytes"]
    assert not report["errors"]


@pytest.mark.parametrize("path", ["../outside", "/etc/passwd", "C:/outside", "C:\\outside", ".git/config", ".venv/package.py", ".", ""])
def test_rejects_unsafe_and_reserved_input_paths(tmp_path, path):
    populate(tmp_path)
    with pytest.raises(ValueError):
        freeze_study(tmp_path, plan(paths=[path]))


def test_rejects_duplicate_missing_directory_and_outside_symlink_inputs(tmp_path):
    root = populate(tmp_path / "root")
    outside = tmp_path / "outside.txt"
    outside.write_bytes(b"outside")
    (root / "escape").symlink_to(outside)
    (root / "alias").symlink_to(root / "method.py")
    for paths in (["model.json", "./model.json"], ["method.py", "alias"], ["escape"]):
        with pytest.raises(ValueError):
            freeze_study(root, plan(paths=paths))
    with pytest.raises(FileNotFoundError):
        freeze_study(root, plan(paths=["missing.txt"]))
    (root / "folder").mkdir()
    with pytest.raises(ValueError, match="regular files"):
        freeze_study(root, plan(paths=["folder"]))
    # An in-root symlink to a single regular input is supported.
    manifest = freeze_study(root, plan(paths=["alias"]))
    assert verify_study(root, manifest)["valid"]
    (root / "alias").unlink()
    (root / "alias").symlink_to(outside)
    report = verify_study(root, manifest)
    assert not report["valid"]
    assert "escapes" in report["errors"][0]


@pytest.mark.parametrize("updates", [
    {"schema_version": 2}, {"schema_version": True}, {"version": 1},
    {"files": []}, {"paths": []}, {"study_id": "../escape"},
    {"evaluation_role": "held_out"}, {"held_out_organisms": ["candidate"]},
    {"development_organisms": ["Btheta"]},
    {"quarantined_organisms": ["Keio"]},
    {"question": float("nan")},
])
def test_rejects_invalid_declarations(tmp_path, updates):
    populate(tmp_path)
    with pytest.raises(ValueError):
        freeze_study(tmp_path, plan(**updates))


def test_external_declaration_requires_frozen_custodian_and_exposure_evidence(tmp_path):
    populate(tmp_path)
    external = plan(evaluation_role="prospective_external")
    with pytest.raises(ValueError, match="candidate_exposure"):
        freeze_study(tmp_path, external)
    external["candidate_exposure"] = {"external_custodian": "Independent custodian",
                                      "uninspected_evidence": "Documented access log and custodian attestation",
                                      "evidence_paths": ["attestation.txt"]}
    with pytest.raises(ValueError, match="included in paths"):
        freeze_study(tmp_path, external)
    (tmp_path / "attestation.txt").write_text("Test evidence declaration, not proof of independence.")
    external["paths"].append("attestation.txt")
    frozen = freeze_study(tmp_path, external)
    assert frozen["plan"]["evaluation_role"] == "prospective_external"
    assert frozen["plan"]["quarantined_organisms"] == ["uninspected-candidate"]
    assert "require independent review" in frozen["evaluation_caveat"]
    assert verify_study(tmp_path, frozen)["valid"]


@pytest.mark.parametrize("mutation", [
    lambda m: m["files"].pop(),
    lambda m: m["files"].append(deepcopy(m["files"][0])),
    lambda m: m["files"][0].pop("sha256"),
    lambda m: m["files"][0].update(sha256="not-a-hash"),
    lambda m: m["files"][0].update(sha256="0" * 64),
    lambda m: m["files"][0].update(size_bytes=-1),
    lambda m: m["files"][0].update(size_bytes=True),
    lambda m: m["files"][0].update(size_bytes=999),
    lambda m: m["plan"].update(question="Changed analysis plan"),
    lambda m: m.update(content_fingerprint="f" * 64),
    lambda m: m.update(schema_version=100),
    lambda m: m.update(schema_version=True),
    lambda m: m.update(created_at=None),
    lambda m: m.update(unsupported_field=True),
    lambda m: m["repository"].update(head="invalid"),
    lambda m: m["repository"].update(tracked_changes=[None]),
])
def test_malformed_or_tampered_manifest_never_passes(tmp_path, mutation):
    populate(tmp_path)
    manifest = freeze_study(tmp_path, plan())
    mutation(manifest)
    report = verify_study(tmp_path, manifest)
    assert not report["valid"]
    assert report["errors"]
    assert report["checked_files"] == 0


def test_removing_a_declared_path_and_its_record_still_breaks_fingerprint(tmp_path):
    populate(tmp_path)
    manifest = freeze_study(tmp_path, plan())
    removed = manifest["files"].pop()["path"]
    manifest["plan"]["paths"].remove(removed)
    report = verify_study(tmp_path, manifest)
    assert not report["valid"]
    assert "fingerprint" in report["errors"][0]


def git(root, *args):
    return subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True, text=True).stdout


def test_git_commit_dirty_and_untracked_provenance_excludes_runtime(tmp_path):
    populate(tmp_path)
    git(tmp_path, "init")
    git(tmp_path, "add", "model.json", "method.py")
    git(tmp_path, "-c", "user.name=Study Test", "-c", "user.email=study-test@example.invalid", "commit", "-m", "Fixture")
    clean = freeze_study(tmp_path, plan())
    assert clean["repository"]["head"] == git(tmp_path, "rev-parse", "HEAD").strip()
    assert clean["repository"]["available"]
    assert clean["repository"]["dirty"] is False
    (tmp_path / ".venv").mkdir()
    (tmp_path / ".venv" / "runtime.txt").write_text("Runtime is excluded from provenance status")
    assert freeze_study(tmp_path, plan())["repository"]["dirty"] is False
    (tmp_path / "notes with spaces.txt").write_text("New development notes")
    (tmp_path / "method.py").write_text("method = 2\n")
    dirty = freeze_study(tmp_path, plan())
    assert dirty["repository"]["dirty"]
    assert dirty["repository"]["tracked_changes"] == [{"path": "method.py", "status": " M"}]
    assert dirty["repository"]["untracked_paths"] == ["notes with spaces.txt"]
    assert verify_study(tmp_path, dirty)["valid"]


def test_cli_roundtrip_rejects_overwriting_and_reports_changes(tmp_path, capsys):
    populate(tmp_path)
    declaration = tmp_path / "plan.json"
    declaration.write_text(json.dumps(plan()))
    output = tmp_path / "manifest.json"
    args = ["freeze", "--root", str(tmp_path), "--plan", str(declaration), "--output", str(output)]
    assert main(args) == 0
    before = output.read_bytes()
    assert main(args) == 2
    assert output.read_bytes() == before
    assert not list(tmp_path.glob(".study-freeze-*.tmp"))
    verification = ["verify", "--root", str(tmp_path), "--manifest", str(output)]
    assert main(verification) == 0
    (tmp_path / "method.py").write_text("method = 3\n")
    assert main(verification) == 1
    capsys.readouterr()


def test_cli_rejects_duplicate_json_keys(tmp_path, capsys):
    populate(tmp_path)
    declaration = tmp_path / "plan.json"
    text = json.dumps(plan())
    declaration.write_text(text[:-1] + ', "schema_version": 1}')
    output = tmp_path / "manifest.json"
    assert main(["freeze", "--root", str(tmp_path), "--plan", str(declaration), "--output", str(output)]) == 2
    assert "Duplicate JSON key" in capsys.readouterr().err
    assert not output.exists()
