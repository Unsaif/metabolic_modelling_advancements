import json, re, collections
H = json.load(open("/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/harvey_ids.json"))
mets, names, hmdb = H["mets"], H["metNames"], H["metHMDBID"]
rx = set(H["rxns"])
# base id -> compartments of body-fluid forms (no organ prefix)
fluid = collections.defaultdict(set); nameof = {}; hm = {}
for i, m in enumerate(mets):
    mm = re.match(r"^([^\[]+)\[(\w+)\]$", m)
    if not mm: continue
    base, comp = mm.groups()
    if comp in ("bc", "u", "csf", "e", "luSI", "fe", "luLI", "bp", "bpL"):
        fluid[base].add(comp)
    nameof.setdefault(base, names[i]); hm.setdefault(base, hmdb[i])
def search(pat):
    out = []
    for base, nm in nameof.items():
        if "_" in base and base.split("_")[0] in ORG: continue
        if re.search(pat, nm, re.I) or re.search(pat, base, re.I):
            out.append((base, nm, hm.get(base), sorted(fluid.get(base, [])), "EX_%s[u]" % base in rx))
    return out
ORG = set()
for m in mets:
    if "_" in m:
        ORG.add(m.split("_")[0])
# organ prefixes are like 'Liver_', 'sIEC_' ...; keep only prefix tokens that start uppercase or 'sIEC'
ORG = {o for o in ORG if o[:1].isupper() or o in ("sIEC",)}
for pat in [r"cystin", r"homocyst", r"hcys", r"methylmalon", r"^HC00900$", r"glutar(ate|ic)", r"hydroxyprolin|4hpro", r"triglycer|triacylglyc|^tag", r"folate|^fol$|5mthf|thf$", r"pyridox|pydx|pydam|pydxn", r"hydroxyisovaler|3hivac|CE2028", r"coproporph", r"aminolevul|5aop", r"lact(ate|ic)", r"^k$|^na1$|^ca2$|potassium|sodium|calcium", r"ammon|^nh4$", r"homovanil", r"hydroxyindol", r"dihydrothym|dihydroura|56dura|56dthm", r"N-acetyl.*tyros|actyr", r"hydroxyphenyl(lact|pyruv|acet)", r"succinylacet", r"isovaleryl", r"deoxycort|cortexolone", r"citrull", r"oxal", r"argininosucc"]:
    res = search(pat)
    print("==", pat)
    for r in res[:40]:
        print("   ", r)
