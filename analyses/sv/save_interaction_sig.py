import pandas as pd, numpy as np

d = pd.read_csv(r"D:\ONT\interaction_results.tsv", sep="\t")
est = d[d.int_p.notna()].copy()
p = est.int_p.values
m = len(p)
order = np.argsort(p); sp = p[order]
cummin = np.minimum.accumulate(sp[::-1] * m / np.arange(m, 0, -1))[::-1]
q = np.empty(m); q[order] = cummin
p_star = float(sp[q < 0.10].max())
sig = est[est.int_p <= p_star].sort_values("int_p")
sig.to_csv(r"D:\ONT\interaction_significant.tsv", sep="\t", index=False)
print("n_sig =", len(sig), " p_star =", p_star)
print(sig[["site","gene_name","window","int_effect","int_p","n_case_car","n_ctrl_car"]].to_string(index=False))
