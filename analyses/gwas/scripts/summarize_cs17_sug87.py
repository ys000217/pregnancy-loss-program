#!/usr/bin/env python3
import shutil
from pathlib import Path
import pandas as pd

cauchy_p = Path("/mnt/d/gsMap/CS17_HESTA/CS17_E1S1_HESTA/cauchy_combination/CS17_E1S1_HESTA_RPL_sug87.Cauchy.csv.gz")
spatial_p = Path("/mnt/d/gsMap/CS17_HESTA/CS17_E1S1_HESTA/spatial_ldsc/CS17_E1S1_HESTA_RPL_sug87.csv.gz")
out = Path("/mnt/d/ONT/figure3/gwas/gsmap_CS17_sug87_results")
out.mkdir(parents=True, exist_ok=True)
res = Path("/mnt/d/ONT/analyses/gwas/results")
res.mkdir(parents=True, exist_ok=True)

shutil.copy2(cauchy_p, out / "CS17_RPL_sug87_cauchy.csv.gz")
shutil.copy2(spatial_p, out / "CS17_RPL_sug87_spatial_ldsc.csv.gz")
shutil.copy2(cauchy_p, res / "CS17_RPL_sug87_cauchy.csv.gz")

c = pd.read_csv(cauchy_p)
df = pd.read_csv(spatial_p)
z2 = df["z"] ** 2
print("==== CAUCHY ====")
print(c.to_string(index=False))
print("==== SPOT ====")
print("n_spots", len(df))
print("mean_z2", float(z2.mean()))
print("median_z2", float(z2.median()))
print("max_abs_z", float(df["z"].abs().max()))
print("min_p", float(df["p"].min()))
bonf = 0.05 / len(df)
print("bonferroni", bonf)
print("n_pass_bonf", int((df["p"] < bonf).sum()))
print("best_p_cauchy", float(c["p_cauchy"].min()), "celltype", c.loc[c["p_cauchy"].idxmin(), "annotation"])

summary = out / "CS17_RPL_sug87_summary.txt"
summary.write_text(
    "\n".join([
        f"best_p_cauchy={c['p_cauchy'].min()}",
        f"best_celltype={c.loc[c['p_cauchy'].idxmin(), 'annotation']}",
        f"n_spots={len(df)}",
        f"mean_z2={z2.mean()}",
        f"median_z2={z2.median()}",
        f"min_spot_p={df['p'].min()}",
        f"n_pass_bonferroni={int((df['p'] < bonf).sum())}",
    ]) + "\n",
    encoding="utf-8",
)
print("wrote", summary)
