"""Benchmark paper Figure 1: study design and time order of the transfer study (3 October 2026).

All times are read from the repository: freeze manifests (results/study_freezes/), outcome download records
(data/fitness_browser_panel/outcomes_download*.json, fetch times recorded in the browser), run cards
(results/transfer_v1/<phase>/<org>/<arm>/card.json, 'created') and the commit times of the panel declaration
and selection. External time records in the Claude project are listed in the results documents.
Writes docs/paper/fig1_timeline.png and results/transfer_v1/timeline.json.
"""
from __future__ import annotations

import glob
import json
import os
import subprocess
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def t(s: str) -> datetime:
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def commit_time(rev: str) -> str:
    return subprocess.run(["git", "-C", ROOT, "log", "-1", "--format=%cI", rev], capture_output=True, text=True).stdout.strip()


def freeze(name: str) -> str:
    return json.load(open(os.path.join(ROOT, "results", "study_freezes", f"{name}.json")))["created_at"]


def downloads(name: str):
    files = json.load(open(os.path.join(ROOT, "data", "fitness_browser_panel", name)))["files"].values()
    return min(f["fetch_started"] for f in files), max(f["fetch_finished"] for f in files)


def runs(role: str, exclude=("M_noIonR6",)):
    created = [json.load(open(f))["created"] for f in glob.glob(os.path.join(ROOT, "results", "transfer_v1", role, "*", "*", "card.json"))
               if f.split(os.sep)[-2] not in exclude]
    return min(created), max(created)


def main() -> None:
    ev = {
        "panel_declared": commit_time("aa25b6d"), "panel_selected": commit_time("8d0d52b"),
        "development_runs": runs("development"), "method_freeze": freeze("transfer_v1_method"),
        "panel_A_inputs_freeze": freeze("transfer_v1_inputs_panel_A"),
        "panel_A_download": downloads("outcomes_download_2026-10-03.json"), "panel_A_runs": runs("evaluation_panel_A"),
        "replication_plan_freeze": freeze("transfer_v1_replication_plan"),
        "panel_B_inputs_freeze": freeze("transfer_v1_inputs_panel_B"),
        "panel_B_download": downloads("outcomes_download_panel_B_2026-10-03.json"), "panel_B_runs": runs("evaluation_panel_B"),
    }
    with open(os.path.join(ROOT, "results", "transfer_v1", "timeline.json"), "w") as fh:
        json.dump(ev, fh, indent=1)
        fh.write("\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 3.6))
    lanes = {"Development\n(4 organisms)": 2, "Panel A\n(6 new organisms)": 1, "Panel B\n(6 new, held back)": 0}
    grey, blue, red = "#7f7f7f", "#1f77b4", "#d62728"

    def span(lane, a, b, colour, text, alpha=0.35):
        ax.barh(lane, mdates.date2num(t(b)) - mdates.date2num(t(a)), left=mdates.date2num(t(a)), height=0.36,
                color=colour, alpha=alpha, edgecolor="none")
        width = mdates.date2num(t(b)) - mdates.date2num(t(a))
        if width < 25 / (24 * 60):   # narrow bar: label to the right
            ax.text(mdates.date2num(t(b)) + 2 / (24 * 60), lane, text, ha="left", va="center", fontsize=7)
        else:
            ax.text(mdates.date2num(t(a)) + width / 2, lane, text, ha="center", va="center", fontsize=7)

    def mark(lane, when, colour, label, marker="v", dy=0.3):
        x = mdates.date2num(t(when))
        ax.plot(x, lane + dy, marker, color=colour, ms=7)
        ax.text(x, lane + dy + 0.12, label, ha="center", va="bottom", fontsize=7, color=colour)

    span(2, *ev["development_runs"], grey, "arms developed and scored")
    mark(2, ev["panel_selected"], "black", "panel drawn\n(metadata only)", marker="o", dy=0.3)
    mark(2, ev["method_freeze"], "black", "method locked", marker="v")
    span(1, ev["method_freeze"], ev["panel_A_inputs_freeze"], blue, "inputs, blind\ncuration", alpha=0.15)
    mark(1, ev["panel_A_inputs_freeze"], blue, "inputs locked")
    span(1, *ev["panel_A_runs"], blue, "scored once")
    ax.plot(mdates.date2num(t(ev["panel_A_download"][0])), 1 - 0.3, "^", color=blue, ms=7)
    ax.text(mdates.date2num(t(ev["panel_A_download"][0])), 1 - 0.48, "fitness first\ndownloaded", ha="center", va="top", fontsize=7, color=blue)
    mark(0, ev["replication_plan_freeze"], red, "second-round\nplan locked")
    span(0, ev["replication_plan_freeze"], ev["panel_B_inputs_freeze"], red, "inputs, blind\ncuration", alpha=0.15)
    mark(0, ev["panel_B_inputs_freeze"], red, "inputs locked", dy=0.3)
    span(0, *ev["panel_B_runs"], red, "scored once")
    ax.plot(mdates.date2num(t(ev["panel_B_download"][0])), 0 - 0.3, "^", color=red, ms=7)
    ax.text(mdates.date2num(t(ev["panel_B_download"][0])), 0 - 0.48, "fitness first\ndownloaded", ha="center", va="top", fontsize=7, color=red)
    ax.set_yticks(list(lanes.values()))
    ax.set_yticklabels(list(lanes.keys()), fontsize=8)
    ax.set_ylim(-0.9, 2.9)
    x0, x1 = ax.get_xlim()
    ax.set_xlim(x0 - 12 / (24 * 60), x1)   # room for the first label
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    ax.tick_params(axis="x", labelsize=8, bottom=False, labelbottom=False)   # order, not clock times, in the paper
    ax.set_xlabel("3 October 2026: events in time order (exact times in Supplementary Table S1)", fontsize=8)
    ax.set_title("Every outcome was first accessed after its inputs were locked", fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "docs", "paper", "fig1_timeline.png"), dpi=200)
    print(json.dumps(ev, indent=1))


if __name__ == "__main__":
    main()
