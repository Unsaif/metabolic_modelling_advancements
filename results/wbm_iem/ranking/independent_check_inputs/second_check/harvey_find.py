import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_v02")
import re
from collections import defaultdict
from common import get_all

d = get_all()
H = d["harvey"]
mets, names, rxns = H["mets"], H["metNames"], H["rxns"]
hmdb = H.get("metHMDBID")
rxset = set(rxns)
# mode: ids <base ids...> | names <regex...>
mode = sys.argv[1]
args = sys.argv[2:]
def show_base(b):
    comps = defaultdict(list)
    nm = set()
    hm = set()
    for i, m in enumerate(mets):
        mm = re.match(r"^(?:([A-Za-z]+)_)?" + re.escape(b) + r"\[(\w+)\]$", m)
        if mm:
            comps[mm.group(2)].append(mm.group(1) or "")
            nm.add(names[i])
            if isinstance(hmdb, list):
                hm.add(hmdb[i])
    ex = [r for r in (f"EX_{b}[u]", f"EX_{b}[fe]", f"EX_{b}[d]") if r in rxset]
    blood = f"{b}[bc]" in set(mets)
    csf = f"{b}[csf]" in set(mets)
    print(f"{b}: names={sorted(nm)} hmdb={sorted(hm)} [bc]={blood} [csf]={csf} [u]-met={f'{b}[u]' in set(mets)} EX={ex}")
    print("     comps:", {k: (len(v), sorted(set(v))[:6]) for k, v in comps.items()})
if mode == "ids":
    for b in args:
        show_base(b)
else:
    for pat in args:
        rx = re.compile(pat, re.I)
        bases = set()
        for i, m in enumerate(mets):
            if rx.search(names[i]):
                mm = re.match(r"^(.*)\[(\w+)\]$", m)
                full = mm.group(1)
                parts = full.split("_", 1)
                bases.add((full, names[i]))
        # collapse organ prefixes
        simple = defaultdict(set)
        for full, n in bases:
            simple[n].add(full)
        print(f"== {pat}")
        for n, fulls in sorted(simple.items()):
            fl = sorted(fulls)
            short = sorted({f for f in fl if not re.match(r"^(Liver|Brain|Kidney|Heart|Muscle|Lung|Adipocytes|Spleen|Pancreas|Stomach|Colon|sIEC|Gall|RBC|Skin|Thyroidgland|Agland|Retina|Urinarybladder|Prostate|Testis|Platelet|Monocyte|Nkcells|CD4Tcells|Bcells|BBB|Breast|Cervix|Ovary|Uterus|Scord|Parathyroidglands|Adrenal|Micro)_", f)})
            print(f"   {n!r}: base ids={short[:8]} (n full={len(fl)})")
