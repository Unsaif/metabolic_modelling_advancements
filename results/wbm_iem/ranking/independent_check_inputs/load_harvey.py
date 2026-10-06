"""Load Harvey 1.03d once and cache mets/metNames/rxns as JSON (independent check)."""
import json, sys, scipy.io as sio
path = "/home/claude/mma/external/COBRA.models/mat/Harvey_1_03d.mat"
d = sio.loadmat(path, squeeze_me=False, struct_as_record=False)
keys = [k for k in d if not k.startswith("__")]
print("keys", keys)
m = d[keys[0]][0, 0]
print("fields", m._fieldnames[:60])
def cellstr(x):
    out = []
    for e in x.ravel():
        if hasattr(e, "size") and e.size == 0:
            out.append("")
        else:
            v = e
            while hasattr(v, "shape") and v.size == 1 and v.dtype == object:
                v = v.ravel()[0]
            out.append(str(v.ravel()[0]) if hasattr(v, "ravel") and v.dtype != object and v.dtype.kind == "U" else str(v))
    return out
mets = cellstr(m.mets); names = cellstr(m.metNames); rxns = cellstr(m.rxns)
print(len(mets), len(names), len(rxns))
print(mets[:5], names[:5], rxns[:5])
json.dump({"mets": mets, "metNames": names, "rxns": rxns},
          open("/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/harvey_ids.json", "w"))
