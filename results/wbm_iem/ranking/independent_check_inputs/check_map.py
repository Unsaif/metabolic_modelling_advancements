"""Check each HPO map row against Harvey: metabolite identity (names/annotations), readout existence."""
import csv, json, re, collections
H = json.load(open("/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/harvey_ids.json"))
mets, names = H["mets"], H["metNames"]
rx = set(H["rxns"])
metset = set(mets)
info = {}
for i, m in enumerate(mets):
    mm = re.match(r"^(?:[A-Za-z]+_)?(.+)\[(\w+)\]$", m)
    info.setdefault(m, (names[i], H["metHMDBID"][i], H["metKEGGID"][i], H["metChEBIID"][i], H["metFormulas"][i]))
rows = list(csv.DictReader(open("/home/claude/mma/data/iem/hpo_term_readout_map_v0.1.tsv"), delimiter="\t"))
print(len(rows), "rows;", collections.Counter(r["decision"] for r in rows))
print("duplicate ids:", [k for k, v in collections.Counter(r["hpo_id"] for r in rows).items() if v > 1])
conv = {"blood": ("DM_", "[bc]"), "urine": ("EX_", "[u]"), "csf": ("DM_", "[csf]")}
for r in rows:
    if r["decision"] != "include":
        continue
    met, fl, ro = r["vmh_metabolite"], r["biofluid"], r["readout"]
    pre, comp = conv[fl]
    expect = f"{pre}{met}{comp}"
    ok_conv = expect == ro
    if ro.startswith("EX_"):
        exists = ro in rx
        # also check the [u] metabolite exists
        exists_met = f"{met}[u]" in metset
    else:
        exists = ro[3:] in metset
        exists_met = exists
    # name of the plain metabolite in any compartment, e.g. [bc], [u], [csf], [c]
    cands = [m for m in (f"{met}[bc]", f"{met}[u]", f"{met}[csf]", f"{met}[e]", f"{met}[c]") if m in info]
    nm = info[cands[0]] if cands else None
    print(f"{r['hpo_id']:11s} {r['hpo_term'][:55]:55s} -> {ro:20s} conv={ok_conv} exists={exists} met_in_comp={exists_met} DMrxn_already={'DM_'+met+comp in rx} | {nm}")
