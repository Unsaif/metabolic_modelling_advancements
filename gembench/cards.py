"""Benchmark cards: the provenance, protocol and leakage record written next to every result."""
from __future__ import annotations

import datetime as _dt
import json
import os
import platform
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional

import cobra


def software_versions() -> Dict[str, str]:
    import numpy, optlang, scipy, sklearn
    return {"python": platform.python_version(), "cobra": cobra.__version__, "optlang": optlang.__version__,
            "numpy": numpy.__version__, "scipy": scipy.__version__, "scikit-learn": sklearn.__version__}


@dataclass
class ModelProvenance:
    model_id: str
    file: str
    source: str
    version_note: str = ""
    n_reactions: int = 0
    n_metabolites: int = 0
    n_genes: int = 0
    sha256: str = ""


@dataclass
class LeakageCard:
    """What the model has 'seen' of the ground truth, and what a frontier model has seen."""
    ground_truth_used_in_model_curation: str      # 'no' | 'yes' | 'partly' | 'unknown', with explanation
    ground_truth_public_since: str                # date/year the data became public
    frontier_model_training_exposure: str         # e.g. 'almost certainly in pre-2026 training corpora'
    held_out_recommendation: str                  # how to obtain a clean held-out estimate
    notes: List[str] = field(default_factory=list)


def carbon_fitness_leakage(org: str, variant: str, *, patched: bool,
                           medium_completion: bool = False) -> LeakageCard:
    """Describe the evaluated arm, including data used after reconstruction."""
    if patched:
        exposure = ("yes: these patch sets were developed and accepted/held/rejected using Fitness Browser "
                    "phenotypes also used for scoring; annotation support does not remove this reuse")
    elif variant == "curated" or org == "Keio":
        exposure = ("partly/unknown: published curated models used phenotype data; exact overlap with "
                    "these measurements has not been audited")
    elif variant == "gapfilled":
        exposure = ("partly: gap filling uses observed wild-type growth in a reference medium; "
                    "the shared media/protocol were refined using these benchmark organisms")
    else:
        exposure = ("no targeted fitness fitting documented for the archived automatic reconstruction; "
                    "however the evaluation media/protocol were refined on these benchmark organisms")
    notes = ["Retrospective development evaluation. No independent held-out or prospective test is established.",
             "Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.",
             "Compare gene-level scores on matched genes, conditions and finite observations; report wild-type growth coverage separately."]
    if medium_completion:
        notes.append("Medium completion adds assumed gene-less transport; recipe presence alone does not establish uptake or its energy cost.")
    return LeakageCard(
        ground_truth_used_in_model_curation=exposure,
        ground_truth_public_since="Fitness Browser releases include Price et al. 2018; exact dates differ by experiment",
        frontier_model_training_exposure="unknown; public accessibility does not establish inclusion in a particular model's training data",
        held_out_recommendation=("Freeze code, media rules and patches before selecting independent data. Use organisms/experiments "
                                 "not inspected during development, or nested evaluation that repeats every data-guided selection step. "
                                 "A retrospective split of already-inspected data does not restore independence."),
        notes=notes)


@dataclass
class BenchmarkCard:
    benchmark: str
    created: str
    dataset_provenance: Dict[str, str]
    model: ModelProvenance
    protocol: Dict[str, Any]
    leakage: LeakageCard
    results: Dict[str, Any]
    software: Dict[str, str] = field(default_factory=software_versions)
    warnings: List[str] = field(default_factory=list)

    def write(self, path_json: str, path_md: Optional[str] = None) -> None:
        os.makedirs(os.path.dirname(path_json), exist_ok=True)
        with open(path_json, "w") as fh:
            json.dump(asdict(self), fh, indent=2, default=_jsonable)
        if path_md:
            with open(path_md, "w") as fh:
                fh.write(self.to_markdown())

    def to_markdown(self) -> str:
        d = asdict(self)
        lines = [f"# Benchmark card — {self.benchmark}", "", f"Created: {self.created}", ""]
        lines += ["## Model", ""] + [f"- {k}: {v}" for k, v in d["model"].items()] + [""]
        lines += ["## Dataset provenance", ""] + [f"- {k}: {v}" for k, v in d["dataset_provenance"].items()] + [""]
        lines += ["## Protocol", ""] + [f"- {k}: {_short(v)}" for k, v in d["protocol"].items()] + [""]
        lines += ["## Leakage", ""] + [f"- {k}: {v}" for k, v in d["leakage"].items()] + [""]
        lines += ["## Results", ""] + [f"- {k}: {_short(v)}" for k, v in d["results"].items()] + [""]
        if self.warnings:
            lines += ["## Warnings", ""] + [f"- {w}" for w in self.warnings] + [""]
        lines += ["## Software", ""] + [f"- {k}: {v}" for k, v in d["software"].items()] + [""]
        return "\n".join(lines)


def now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _short(v: Any) -> str:
    s = json.dumps(v, default=_jsonable) if not isinstance(v, str) else v
    return s if len(s) < 400 else s[:400] + " …"


def _jsonable(o: Any):
    try:
        import numpy as np
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
    except Exception:
        pass
    if hasattr(o, "__dict__"):
        return o.__dict__
    return str(o)


def sha256_of(path: str) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
