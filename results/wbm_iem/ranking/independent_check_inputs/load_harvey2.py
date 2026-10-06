import json, scipy.io as sio
path = "/home/claude/mma/external/COBRA.models/mat/Harvey_1_03d.mat"
d = sio.loadmat(path, squeeze_me=False, struct_as_record=False)
m = d["male"][0, 0]
def cellstr(x):
    out = []
    for e in x.ravel():
        if hasattr(e, "size") and e.size == 0:
            out.append(""); continue
        v = e
        while hasattr(v, "dtype") and v.dtype == object and v.size >= 1:
            v = v.ravel()[0]
        if hasattr(v, "dtype"):
            out.append(str(v.ravel()[0]) if v.size else "")
        else:
            out.append(str(v))
    return out
res = {}
for f in ["mets", "metNames", "rxns", "metHMDBID", "metKEGGID", "metChEBIID", "metFormulas", "rxnNames"]:
    res[f] = cellstr(getattr(m, f))
    print(f, len(res[f]), res[f][:3])
json.dump(res, open("/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/harvey_ids.json", "w"))
