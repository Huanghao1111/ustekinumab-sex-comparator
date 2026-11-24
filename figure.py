from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[0]
df = pd.read_csv(ROOT / "results/core_estimates_bh.csv")
df = df[(df["scenario"] == "all_ages") & df["outcome"].isin(["depression_pt", "anxiety_pt"])]
comparators = [f"C{i}" for i in range(1, 7)]
styles = {
    "faers": ("FAERS", "#174A7E", "o", -0.12),
    "canada": ("Canada Vigilance", "#D97732", "D", 0.12),
}
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "pdf.fonttype": 42, "svg.fonttype": "none"})
fig, axes = plt.subplots(1, 2, figsize=(183 / 25.4, 120 / 25.4))
fig.subplots_adjust(left=0.105, right=0.985, bottom=0.15, top=0.82, wspace=0.32)
y_base = np.arange(5, -1, -1, dtype=float)

for ax, outcome, title, label in zip(axes, ["depression_pt", "anxiety_pt"], ["Depression", "Anxiety"], ["a", "b"]):
    ax.set(xscale="log", xlim=(0.22, 60), ylim=(-0.65, 5.65), yticks=y_base, yticklabels=comparators)
    ax.spines[["top", "right"]].set_visible(False)
    ax.axhspan(1.58, 2.42, color="#EAF1F7", zorder=0)
    ax.axvline(1, color="#777777", linewidth=0.8, linestyle=(0, (3, 3)), zorder=1)
    for database, (_, color, marker, offset) in styles.items():
        part = df[(df["database"] == database) & (df["outcome"] == outcome)].set_index("comparator").loc[comparators]
        x = part["RORR"].to_numpy()
        ax.errorbar(x, y_base + offset, xerr=np.vstack([x - part["RORR_lci"], part["RORR_uci"] - x]),
                    fmt=marker, color=color, ecolor=color, markersize=4.1, markeredgecolor="white",
                    markeredgewidth=0.45, elinewidth=1, capsize=2, linestyle="none", zorder=3)
    ax.set_xticks([0.25, 0.5, 1, 2, 5, 10, 20, 50], ["0.25", "0.5", "1", "2", "5", "10", "20", "50"])
    ax.set_xlabel("Ratio of reporting odds ratios (95% CI)")
    ax.set_title(title, fontweight="bold", pad=7)
    ax.text(-0.155, 1.045, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")
    ax.tick_params(direction="out", length=3, width=0.7, pad=3)
    ax.get_yticklabels()[3].set_fontweight("bold")
    ax.get_yticklabels()[3].set_color("#174A7E")

handles = [Line2D([0], [0], marker=marker, color=color, markeredgecolor="white", markersize=4.5,
                  linewidth=1, label=label) for label, color, marker, _ in styles.values()]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.55, 0.975), ncol=2, frameon=False)
fig.savefig(ROOT / "results/Figure_1.pdf")
fig.savefig(ROOT / "results/Figure_1.png", dpi=600)
print("Wrote Figure_1.pdf and Figure_1.png")
