"""Shared helpers for the independent verification of IEM study v0.4.

Written by the independent verifier from the definitions in the task statement and the protocol text.
Imports nothing from gembench or from the project's scripts/*.py.
"""
import hashlib
import json
import math
import os
from collections import OrderedDict

import numpy as np

ROOT = "/home/claude/mma"
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, "results/wbm_iem")
PROTOCOL = os.path.join(ROOT, "data/iem/iem_protocol_v0.2.json")

FILES = OrderedDict([
    ("harvey_max_v0.3", "results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json"),
    ("harvey_max_v0.2", "results/wbm_iem/Harvey_1_03d_iem_results_v0.2.json"),
    ("harvey_max_v0.2b", "results/wbm_iem/Harvey_1_03d_iem_results_v0.2b.json"),
    ("harvetta_max_A", "results/wbm_iem/Harvetta_1_03d_iem_results_v0.4.json"),
    ("harvey_min_B", "results/wbm_iem/Harvey_1_03d_iem_results_v0.4_min.json"),
    ("harvetta_min_C", "results/wbm_iem/Harvetta_1_03d_iem_results_v0.4_min.json"),
])
TOL = 1e-6
DIRS = ("Increased", "Decreased")


def load_json(rel):
    with open(os.path.join(ROOT, rel)) as fh:
        return json.load(fh)


def sha256_file(rel):
    h = hashlib.sha256()
    with open(os.path.join(ROOT, rel), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def expected_label(text):
    """runIEM_HH: 'Incre' -> Increased, 'Decre' -> Decreased, otherwise no direction."""
    if "Incre" in text:
        return "Increased"
    if "Decre" in text:
        return "Decreased"
    return "Unchanged"


def finite(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def z(v):
    """Protocol: values with |f| <= 1e-6 are set to 0."""
    return 0.0 if abs(v) <= TOL else v


def call_protocol(h, d, zero=True):
    if zero:
        h, d = z(h), z(d)
    diff = d - h
    if diff > TOL:
        return "Increased"
    if diff < -TOL:
        return "Decreased"
    return "Unchanged"


def rel_change(h, d, zero=False):
    if zero:
        h, d = z(h), z(d)
    s = max(abs(h), abs(d))
    return (d - h) / s if s > 0 else 0.0


def call_rel(h, d, tau, zero=False):
    r = rel_change(h, d, zero)
    if r > tau:
        return "Increased"
    if r < -tau:
        return "Decreased"
    return "Unchanged"


def protocol_rows(protocol):
    """Flat list of (iem, call_index, position, reaction, label, expected) in protocol order."""
    out = []
    for p in protocol:
        for k, (rxn, text) in enumerate(p["biomarkers"]):
            out.append(dict(iem=p["iem"], call_index=p["call_index"], pos=k, reaction=rxn, label=text,
                            expected=expected_label(text)))
    return out


def index_results(records):
    by = {}
    for r in records:
        key = (r["iem"], r["call_index"])
        if key in by:
            raise SystemExit(f"duplicate record {key}")
        by[key] = r
    return by


def biomarker(by, row):
    rec = by[(row["iem"], row["call_index"])]
    b = rec["biomarkers"][row["pos"]]
    if b["reaction"] != row["reaction"]:
        raise SystemExit(f"position mismatch {row}")
    return b


def max_ok(b):
    return (b.get("status_healthy") == "Optimal" and b.get("status_disease") == "Optimal"
            and finite(b.get("healthy")) and finite(b.get("disease")))


def min_ok(b):
    return (b.get("status_healthy_min") == "Optimal" and b.get("status_disease_min") == "Optimal"
            and finite(b.get("healthy_min")) and finite(b.get("disease_min")))


# ---------------------------------------------------------------- MATLAB parsing (runIEM_HH printed values)
def cell_text(c):
    a = np.asarray(c)
    if a.size == 0:
        return ""
    if a.dtype.kind in "US":
        return str(a.ravel()[0])
    raise TypeError(f"unexpected cell content {a!r}")


def str2num(s):
    """MATLAB str2num: a number, NaN, or None for [] (empty, 'NA', 'ND', ...)."""
    s = s.strip()
    if s == "":
        return None
    low = s.lower()
    if low == "nan":
        return float("nan")
    if low in ("inf", "+inf"):
        return float("inf")
    if low == "-inf":
        return float("-inf")
    try:
        return float(s)
    except ValueError:
        return None


def matlab_call(h, d):
    """runIEM_HH: H_D = H - D; Increased if H_D < -1e-6, Decreased if H_D > 1e-6; NaN or [] -> unchanged."""
    if h is None or d is None:
        return "Unchanged"
    hd = h - d
    if hd < -TOL:
        return "Increased"
    if hd > TOL:
        return "Decreased"
    return "Unchanged"


def parse_matlab_iem(path):
    import scipy.io as sio
    d = sio.loadmat(path, squeeze_me=False)
    out = OrderedDict()
    for key in sorted(k for k in d if k.startswith("IEMSol_")):
        cells = np.asarray(d[key], dtype=object)
        header = [[cell_text(cells[i, j]) for j in range(cells.shape[1])] for i in range(4)]
        rows = []
        for j in range(4, cells.shape[0], 2):
            hname, dname = cell_text(cells[j, 0]), cell_text(cells[j + 1, 0])
            if not (hname.startswith("Healthy:") and dname.startswith("Disease:")):
                raise SystemExit(f"unexpected row names {key} {j}: {hname} / {dname}")
            rxn = hname[len("Healthy:"):]
            if dname[len("Disease:"):] != rxn:
                raise SystemExit(f"row pair mismatch {key} {j}")
            lab = cell_text(cells[j, 2])
            prefix = next((p for p in ("Disease - Reported:", "Healthy - Reported:") if lab.startswith(p)), None)
            if prefix is None:
                raise SystemExit(f"unexpected label {key} {j}: {lab}")
            hs, ds = cell_text(cells[j, 1]), cell_text(cells[j + 1, 1])
            rows.append(dict(reaction=rxn, label=lab[len(prefix):], label_prefix=prefix, healthy_str=hs,
                             disease_str=ds, healthy=str2num(hs), disease=str2num(ds)))
        out[key[len("IEMSol_"):]] = dict(header=header, rows=rows)
    scal = {}
    for k in ("Accuracy", "Precision", "FalseDiscoveryRate", "NumBiomarkers", "NumDiseases"):
        if k in d:
            scal[k] = float(np.asarray(d[k]).ravel()[0])
    return out, scal


def mfin(v):
    return v is not None and math.isfinite(v)


# ---------------------------------------------------------------- model / bounds
def cellstr(a):
    out = []
    for x in np.asarray(a, dtype=object).ravel():
        x = np.asarray(x)
        out.append(str(x.ravel()[0]) if x.size else "")
    return out


def load_model_bounds(path):
    import scipy.io as sio
    d = sio.loadmat(path, squeeze_me=False, struct_as_record=True,
                    variable_names=None)
    keys = [k for k in d if not k.startswith("__")]
    if len(keys) != 1:
        raise SystemExit(f"model file has variables {keys}")
    m = d[keys[0]][0, 0]
    return (keys[0], cellstr(m["rxns"]), np.asarray(m["lb"], dtype=np.float64).ravel(),
            np.asarray(m["ub"], dtype=np.float64).ravel(), m)


def load_bounds(path):
    import scipy.io as sio
    d = sio.loadmat(path, squeeze_me=False)
    return cellstr(d["rxns"]), np.asarray(d["lb"], dtype=np.float64).ravel(), np.asarray(d["ub"], dtype=np.float64).ravel()


def sha_bounds(lb, ub):
    return hashlib.sha256(np.ascontiguousarray(lb, dtype="<f8").tobytes() +
                          np.ascontiguousarray(ub, dtype="<f8").tobytes()).hexdigest()


def _strkeys(o):
    if isinstance(o, dict):
        return {(k if isinstance(k, (str, int, float, bool)) or k is None else " | ".join(map(str, k))): _strkeys(v)
                for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_strkeys(v) for v in o]
    return o


def dump(name, obj):
    path = os.path.join(HERE, name)
    with open(path, "w") as fh:
        json.dump(_strkeys(obj), fh, indent=1, default=str)
        fh.write("\n")
    return path
