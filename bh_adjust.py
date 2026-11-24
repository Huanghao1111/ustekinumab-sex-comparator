from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[0]
source = ROOT / "results/core_estimates.csv"
out = ROOT / "results/core_estimates_bh.csv"
df = pd.read_csv(source)
p = df["p_value"].astype(float)
n = len(p)
order = p.sort_values(kind="mergesort").index
ranked = pd.Series(p.loc[order].to_numpy() * n / pd.Series(range(1, n + 1)).to_numpy(), index=order)
adjusted = ranked.iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
df["p_BH_72"] = adjusted.loc[df.index].to_numpy()
df["analysis_family_n"] = n
df.to_csv(out, index=False)
print(f"Wrote {out}")
