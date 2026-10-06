"""Independent loaders for the v0.2 check (written by the checking agent; not derived from the repo builder)."""
import json
import os
import pickle
import xml.etree.ElementTree as ET

SCR = "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad"
ORPHA = os.path.join(SCR, "dl_orphadata", "en_product4.xml")
OBO = os.path.join(SCR, "dl_hpo", "hp.obo")
HARVEY = "/home/claude/mma/external/COBRA.models/mat/Harvey_1_03d.mat"
REPO = "/home/claude/mma"
CACHE = os.path.join(SCR, "independent_check_v02", "cache.pkl")


def parse_orphanet(path=ORPHA):
    tree = ET.parse(path)
    root = tree.getroot()
    out = {}
    # Walk HPODisorderSetStatus -> Disorder explicitly
    n_disorder_elems = 0
    for st in root.iter("HPODisorderSetStatus"):
        for d in st.findall("Disorder"):
            n_disorder_elems += 1
            code = d.find("OrphaCode").text.strip()
            name = d.find("Name").text
            dtype = d.find("DisorderType/Name")
            dgroup = d.find("DisorderGroup/Name")
            anns = []
            lst = d.find("HPODisorderAssociationList")
            if lst is not None:
                for a in lst.findall("HPODisorderAssociation"):
                    hid = a.find("HPO/HPOId").text.strip()
                    term = a.find("HPO/HPOTerm").text
                    fr = a.find("HPOFrequency/Name")
                    anns.append((hid, term, fr.text if fr is not None else ""))
            if code in out:
                out[code]["dup"] = out[code].get("dup", 0) + 1
                out[code]["annotations"] += anns
            else:
                out[code] = {"name": name, "type": dtype.text if dtype is not None else None,
                             "group": dgroup.text if dgroup is not None else None, "annotations": anns}
    return out, n_disorder_elems, root.attrib


def parse_obo(path=OBO):
    terms = {}
    cur = None
    in_term = False
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("["):
                in_term = (line.strip() == "[Term]")
                cur = None
                continue
            if not in_term:
                continue
            if line.startswith("id: "):
                cur = line[4:].strip()
                terms[cur] = {"name": None, "def": None, "is_a": [], "alt_id": [], "obsolete": False,
                              "replaced_by": [], "synonyms": [], "consider": []}
            elif cur is None:
                continue
            elif line.startswith("name: "):
                terms[cur]["name"] = line[6:].strip()
            elif line.startswith("def: "):
                terms[cur]["def"] = line[5:].strip()
            elif line.startswith("is_a: "):
                terms[cur]["is_a"].append(line[6:].split("!")[0].strip())
            elif line.startswith("alt_id: "):
                terms[cur]["alt_id"].append(line[8:].strip())
            elif line.startswith("is_obsolete: true"):
                terms[cur]["obsolete"] = True
            elif line.startswith("replaced_by: "):
                terms[cur]["replaced_by"].append(line[13:].strip())
            elif line.startswith("consider: "):
                terms[cur]["consider"].append(line[10:].strip())
            elif line.startswith("synonym: "):
                terms[cur]["synonyms"].append(line[9:].strip())
    return terms


def ancestors(tid, terms):
    seen = set()
    stack = [tid]
    while stack:
        x = stack.pop()
        for p in terms.get(x, {}).get("is_a", []):
            if p not in seen:
                seen.add(p)
                stack.append(p)
    return seen


def load_harvey(path=HARVEY):
    import scipy.io as sio
    d = sio.loadmat(path, squeeze_me=False, struct_as_record=False)
    key = [k for k in d if not k.startswith("__")][0]
    m = d[key][0, 0]

    def cells(arr):
        out = []
        for x in arr.ravel():
            if hasattr(x, "size") and x.size:
                out.append(str(x.ravel()[0]) if x.dtype.kind in "OU" and x.ndim > 0 else str(x))
            else:
                out.append("")
        return out

    fields = m._fieldnames
    res = {"fields": fields}
    for f in ("mets", "metNames", "rxns", "rxnNames", "metHMDBID", "metKEGGID", "metFormulas", "metCharges",
              "metHMDB", "metKEGG", "metCHEBIID", "subSystems", "grRules"):
        if f in fields:
            try:
                res[f] = cells(getattr(m, f))
            except Exception as e:  # noqa
                res[f] = ("ERR", repr(e))
    return res


def get_all(force=False):
    if os.path.exists(CACHE) and not force:
        with open(CACHE, "rb") as fh:
            return pickle.load(fh)
    orpha, nd, attrib = parse_orphanet()
    terms = parse_obo()
    harvey = load_harvey()
    data = {"orpha": orpha, "orpha_n_disorder_elems": nd, "orpha_attrib": attrib, "terms": terms, "harvey": harvey}
    with open(CACHE, "wb") as fh:
        pickle.dump(data, fh)
    return data


if __name__ == "__main__":
    d = get_all(force=True)
    print("orphanet disorders", len(d["orpha"]), "elements", d["orpha_n_disorder_elems"], d["orpha_attrib"])
    print("dups", [c for c, v in d["orpha"].items() if v.get("dup")])
    print("hpo terms", len(d["terms"]))
    h = d["harvey"]
    print("harvey fields", h["fields"])
    for f in ("mets", "metNames", "rxns", "rxnNames"):
        v = h.get(f)
        print(f, len(v) if isinstance(v, list) else v, (v[:3] if isinstance(v, list) else ""))
