#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Full-cohort narrative: add 30 clustering suspects to abnormal, re-run Fig.1–5.

Does not overwrite plots/narrative (8–10 week set) or tables/gw8_10.

  analyses/burden/.venv/Scripts/python.exe scripts/run_narrative_all_suspect.py
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

MODULE = Path(__file__).resolve().parents[1]
SCRIPTS = Path(__file__).resolve().parent
PY = sys.executable
PHENO_SRC = Path("D:/ONT/figure2/sample_phenotype_648.tsv")
SUSPECT_TSV = Path(
    "D:/ONT/筛选高变CpG/clustering_output/NC_outlier_gestational_week/outlier30_samples.tsv"
)
OUT = MODULE / "tables" / "suspect_abn"
PLOT = MODULE / "plots" / "narrative_all_suspect"
LONGS = MODULE / "tables" / "epifactors"


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def build_phenotype() -> tuple[Path, list[str]]:
    suspects = [r["sample"].strip() for r in read_tsv(SUSPECT_TSV)]
    pheno = read_tsv(PHENO_SRC)
    by_id = {r["VCF_Sample_ID"]: r for r in pheno}
    lower = {k.lower(): k for k in by_id}
    matched, missing, already_ab = [], [], []
    for sid in suspects:
        key = sid if sid in by_id else lower.get(sid.lower())
        if key is None:
            missing.append(sid)
            continue
        matched.append(key)
        if by_id[key]["Condition"] == "abnormal":
            already_ab.append(key)
        by_id[key]["Condition"] = "abnormal"
        by_id[key]["Suspect_Abnormal"] = "1"
    for r in pheno:
        r.setdefault("Suspect_Abnormal", "0")
        if r["VCF_Sample_ID"] in matched:
            r["Suspect_Abnormal"] = "1"
    fields = list(pheno[0].keys())
    if "Suspect_Abnormal" not in fields:
        fields.append("Suspect_Abnormal")
    out_p = OUT / "phenotype.tsv"
    write_tsv(out_p, pheno, fields)
    n = Counter(r["Condition"] for r in pheno)
    notes = [
        f"suspects_listed\t{len(suspects)}",
        f"suspects_matched\t{len(matched)}",
        f"suspects_missing\t{len(missing)}\t{','.join(missing)}",
        f"suspects_already_abnormal\t{len(already_ab)}",
        f"n_abnormal\t{n.get('abnormal', 0)}",
        f"n_normal\t{n.get('normal', 0)}",
        f"n_control\t{n.get('control', 0)}",
        f"n_total\t{len(pheno)}",
        "source\t" + str(SUSPECT_TSV),
        "note\t30 NC-cluster suspects (weeks 7/11/12) relabeled abnormal; full 648, no GW filter",
    ]
    (OUT / "README.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")
    write_tsv(
        OUT / "suspect_ids.tsv",
        [{"VCF_Sample_ID": s, "Condition": "abnormal"} for s in matched],
        ["VCF_Sample_ID", "Condition"],
    )
    print("\n".join(notes))
    if missing:
        raise SystemExit(f"unmatched suspect IDs: {missing}")
    return out_p, matched


def run(cmd: list[str], extra_env: dict[str, str] | None = None) -> None:
    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)
    print(">", " ".join(cmd))
    subprocess.check_call(cmd, cwd=str(MODULE), env=env)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PLOT.mkdir(parents=True, exist_ok=True)
    pheno, _ = build_phenotype()
    env = {
        "BURDEN_NARRATIVE": "all_suspect",
        "BURDEN_PLOT": str(PLOT),
    }
    run(
        [
            PY, str(SCRIPTS / "compute_sv_locus_enrichment.py"),
            "--pheno", str(pheno),
            "--out", str(OUT),
            "--skip-all-pass",
        ],
        env,
    )
    run(
        [
            PY, str(SCRIPTS / "epifactors_catalog.py"),
            "--plot-only",
            "--from-longs", str(LONGS),
            "--derived-out", str(OUT / "epifactors"),
            "--pheno", str(pheno),
        ],
        env,
    )
    run([PY, str(SCRIPTS / "plot_narrative.py")], env)
    print(f"done -> {PLOT}")


if __name__ == "__main__":
    main()
