"""Cross-check of check_energy_gate.py with the repository's own code: gembench's medium completion and
energy_from_nothing, applied to both BiGG views exactly as gembench.protocols.carbon_fitness_generic.run applies
them (this is the only script here that imports gembench; it is a cross-check, not the check).

Usage: python -I crosscheck_gate_with_repo_code.py   (writes crosscheck_gate_with_repo_code.json)
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import MODELS, ROOT, dump, read_model  # noqa: E402

sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import run_transfer_study as RT  # noqa: E402
from gembench.checks import energy_from_nothing  # noqa: E402
from gembench.protocols.carbon_fitness_generic import complete_medium_transport  # noqa: E402


def main():
    out = {}
    for label, cfg in MODELS.items():
        m = read_model(cfg["view"])
        with RT.reference_tables(RT.load_config(cfg["org"])):
            for ex in m.exchanges:
                ex.bounds = (0.0, 1000.0)
            added = complete_medium_transport(m, [cfg["medium"]], ["pnto__R", "fol", "hco3"])
            egc = energy_from_nothing(m)
        out[label] = {"medium_completion_added": added, "energy_from_nothing": egc,
                      "gate_fails": any(v > 1e-6 for v in egc.values()),
                      "currencies_tested": sorted(egc)}
        print(label, out[label], flush=True)
    dump(out, "crosscheck_gate_with_repo_code.json")


if __name__ == "__main__":
    main()
