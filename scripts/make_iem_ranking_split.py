"""Split the 57 simulated IEMs into a development set and a held-out set for work on the healthy reference state.

Usage:
  python scripts/make_iem_ranking_split.py   # writes data/iem/iem_ranking_split_v1.json

Stratified systematic sampling: the IEMs are ordered by stratum (whether the lab profile contains a known decrease,
then the profile size: 1, 2-3, 4-6, 7 or more biomarkers) and, within a stratum, by a seeded random key; every third
IEM from a seeded random offset goes to development. This gives 19 development and 38 held-out IEMs, each stratum
split close to one in three. The split depends only on the lab profiles and the seed, not on any prediction.
"""
from __future__ import annotations

import hashlib
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROFILES = os.path.join(ROOT, "data", "iem", "iem_ranking_profiles_v0.2.json")
OUT = os.path.join(ROOT, "data", "iem", "iem_ranking_split_v1.json")
SEED = 20261007


def size_bin(n):
    return "1" if n == 1 else "2-3" if n <= 3 else "4-6" if n <= 6 else "7+"


def make_split(lab, seed=SEED):
    rng = np.random.RandomState(seed)
    iems = sorted(lab)
    keys = dict(zip(iems, rng.random_sample(len(iems))))
    stratum = {i: (any(d == "Decreased" for _, d in lab[i]), size_bin(len(lab[i]))) for i in iems}
    order = sorted(iems, key=lambda i: (stratum[i][0], ["1", "2-3", "4-6", "7+"].index(stratum[i][1]), keys[i]))
    offset = int(rng.randint(3))
    dev = sorted(order[offset::3])
    held = sorted(i for i in iems if i not in dev)
    strata = {}
    for i in iems:
        k = f"{'decrease' if stratum[i][0] else 'increases only'}, {stratum[i][1]}"
        s = strata.setdefault(k, {"development": 0, "held_out": 0})
        s["development" if i in dev else "held_out"] += 1
    return dev, held, strata, offset


def main():
    with open(PROFILES) as fh:
        lab = json.load(fh)["profiles"]["lab"]
    dev, held, strata, offset = make_split(lab)
    out = {"seed": SEED, "offset": offset, "method": __doc__.split("\n\n")[2].replace("\n", " "),
           "profiles_file": os.path.relpath(PROFILES, ROOT),
           "profiles_sha256": hashlib.sha256(open(PROFILES, "rb").read()).hexdigest(),
           "development": dev, "held_out": held, "strata": strata}
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"{len(dev)} development, {len(held)} held out -> {os.path.relpath(OUT, ROOT)}")
    for k, v in sorted(strata.items()):
        print(f"  {k}: {v}")
    print("development:", " ".join(dev))


if __name__ == "__main__":
    main()
