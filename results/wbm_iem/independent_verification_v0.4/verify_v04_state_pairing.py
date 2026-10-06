"""Check that each minima run shares its healthy/disease states with the maxima run it is paired with in A5:
same IEM reaction sets and the same IEM-flux maxima (vmax_healthy, vmax_disease). Writes state_pairing_v04.json."""
import sys
from vcommon import dump, index_results, load_json

PAIRS = {"Harvey (v0.3 max vs B min)": ("results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json", "results/wbm_iem/Harvey_1_03d_iem_results_v0.4_min.json"),
         "Harvetta (A max vs C min)": ("results/wbm_iem/Harvetta_1_03d_iem_results_v0.4.json", "results/wbm_iem/Harvetta_1_03d_iem_results_v0.4_min.json")}


def main():
    out = {}
    for name, (pmax, pmin) in PAIRS.items():
        a, b = index_results(load_json(pmax)), index_results(load_json(pmin))
        same_rxns = sum(sorted(a[k]["iem_reactions"]) == sorted(b[k]["iem_reactions"]) for k in a)
        rel = []
        for k in a:
            for f in ("vmax_healthy", "vmax_disease"):
                x, y = a[k][f], b[k][f]
                if x is None or y is None:
                    rel.append((k[0], f, x, y, None))
                else:
                    rel.append((k[0], f, x, y, abs(x - y) / max(abs(x), abs(y), 1e-12)))
        finite = [r for r in rel if r[4] is not None]
        out[name] = dict(n_iems=len(a), same_iem_reaction_sets=same_rxns,
                         max_rel_diff_vmax=max(r[4] for r in finite), worst=sorted(finite, key=lambda r: -r[4])[:3],
                         n_rel_over_1e6=sum(r[4] > 1e-6 for r in finite), nonfinite=[r for r in rel if r[4] is None])
        print(name, {k: v for k, v in out[name].items() if k != "worst"}, out[name]["worst"])
    dump("state_pairing_v04.json", out)


if __name__ == "__main__":
    sys.exit(main())
