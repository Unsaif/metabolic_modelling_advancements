import json, re, scipy.io as sio, numpy as np
d = sio.loadmat("/home/claude/mma/external/COBRA.models/mat/Harvey_1_03d.mat", squeeze_me=False, struct_as_record=False)
m = d["male"][0, 0]
H = json.load(open("/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/harvey_ids.json"))
mets, rxns = H["mets"], H["rxns"]
import scipy.sparse as sp
S = sp.csr_matrix(m.S)
for base in ("3hivac", "CE2028"):
    idx = [i for i, x in enumerate(mets) if re.match(rf"^(?:[A-Za-z0-9]+_)?{base}\[", x)]
    rx = set()
    for i in idx:
        rx |= set(S[i].indices.tolist())
    names = sorted(rxns[j] for j in rx)
    nonlumen = [n for n in names if not re.search(r"(luSI|luLI|luC|\[fe\]|_fe\]|\[d\]|EX_|Diet_|Excretion|\[lu\])", n)]
    print(base, "mets:", len(idx), "reactions:", len(names))
    print("   e.g.", names[:25])
