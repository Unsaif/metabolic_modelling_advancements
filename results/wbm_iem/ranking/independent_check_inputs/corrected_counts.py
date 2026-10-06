import json, copy
P = json.load(open("/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/my_profiles.json"))
def add(prof, iem, tuples):
    s = {tuple(t) for t in prof.get(iem, [])}; s |= set(tuples); prof[iem] = sorted(s)
def count(p): return len([v for v in p.values() if v]), sum(len(v) for v in p.values())
hpo, hf = copy.deepcopy(P["hpo"]), copy.deepcopy(P["hpo_frequent"])
# (1) the four lost annotations
add(hpo, "DPYR", [("EX_56dura[u]", "Increased"), ("EX_thym[u]", "Increased")]); add(hf, "DPYR", [("EX_56dura[u]", "Increased"), ("EX_thym[u]", "Increased")])
add(hpo, "HLYS1", [("DM_citr_L[bc]", "Increased"), ("EX_citr_L[u]", "Increased")])
print("with lost annotations restored: hpo", count(hpo), "hpo_frequent", count(hf))
# (2) directed metabolite terms the v0.1 regex does not parse (children of mapped terms, or lactic acidosis)
h2, f2 = copy.deepcopy(hpo), copy.deepcopy(hf)
add(h2, "HMG", [("DM_glc_D[bc]", "Decreased")]); add(f2, "HMG", [("DM_glc_D[bc]", "Decreased")])   # Nonketotic (VF) / Recurrent (F) hypoglycemia
add(h2, "GA1", [("DM_glc_D[bc]", "Decreased")])                                                    # Fasting hypoglycemia (Occasional)
add(h2, "IVA", [("DM_lac_L[bc]", "Increased")]); add(f2, "IVA", [("DM_lac_L[bc]", "Increased")])   # Lactic acidosis (Frequent)
print("plus parser-missed directed terms: hpo", count(h2), "hpo_frequent", count(f2))
prot = json.load(open("/home/claude/mma/data/iem/iem_protocol_v0.2.json"))
panel = {r for p in prot for r, _ in p["biomarkers"]}
print("new readouts needed outside protocol panel:", sorted({r for v in list(h2.values()) + list(f2.values()) for r, _ in v} - panel - set(json.load(open('/home/claude/mma/data/iem/iem_ranking_profiles_v0.1.json'))['extra_readouts'])))
