#!/usr/bin/env python3
"""Apply loose dual-filter SNV rules to cohort ANNOVAR + VCF; write summary + candidates."""
from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from snv_loose_filter_rules import (  # noqa: E402
    _parse_float,
    clinvar_ok_loose,
    consequence_class,
    decide_loose,
    function_ok_loose,
    pop_af_ok_loose,
    stratum_label,
)

ANNO = Path(r"D:/ONT/figure2/S956.snp.annovar.hg38_multianno.txt.gz")
VCF = Path(r"D:/ONT/figure2/WGS_ONT_Intersection_648samples.vcf.gz")
PHENO = Path(r"D:/ONT/figure2/sample_phenotype_648.tsv")
OUT_DIR = ROOT / "tables" / "loose_dual_filter"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def norm_chr(c: str) -> str:
    c = str(c).strip()
    return c[3:] if c.lower().startswith("chr") else c


def parse_af(x) -> float | None:
    if x is None or str(x).strip() in {"", ".", "NA", "na", "nan", "-"}:
        return None
    try:
        return float(x)
    except ValueError:
        return _parse_float(x)


def load_pheno(path: Path) -> dict[str, str]:
    pheno: dict[str, str] = {}
    with open(path, encoding="utf-8") as f:
        h = f.readline().rstrip("\n").split("\t")
        si, ci = h.index("VCF_Sample_ID"), h.index("Condition")
        for line in f:
            p = line.rstrip("\n").split("\t")
            c = p[ci].strip().lower()
            if c in {"abnormal", "case"}:
                c = "abnormal"
            elif c in {"normal", "normal_case"}:
                c = "normal"
            pheno[p[si]] = c
    return pheno


def main() -> None:
    pheno = load_pheno(PHENO)
    controls = {s for s, g in pheno.items() if g == "control"}
    cases = {s for s, g in pheno.items() if g == "abnormal"}
    normals = {s for s, g in pheno.items() if g == "normal"}
    print(
        f"PHENOTYPE abnormal={len(cases)} normal={len(normals)} "
        f"control={len(controls)} total={len(pheno)}",
        flush=True,
    )

    # Pass 1: ANNOVAR prefilter (function + clinvar + pop AF)
    sites: dict[tuple, dict] = {}
    funnel = Counter()
    with gzip.open(ANNO, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        need = ["CHROM", "POS", "REF", "ALT", "Func", "ExonicFunc", "CLNSIG", "GnomAD_genome_AF_ALL"]
        for c in need:
            if c not in idx:
                raise SystemExit(f"missing column {c}")
        optional = {
            "Gene": idx.get("Gene"),
            "GeneName": idx.get("GeneName"),
            "SIFT": idx.get("SIFT"),
            "Polyphen2_HDIV": idx.get("Polyphen2_HDIV"),
            "Polyphen2_HVAR": idx.get("Polyphen2_HVAR"),
            "CADD": idx.get("CADD"),
            "REVEL": idx.get("REVEL"),
        }
        for line in f:
            funnel["annovar_rows"] += 1
            p = line.rstrip("\n").split("\t")
            func = p[idx["Func"]]
            exonic = p[idx["ExonicFunc"]]
            if not function_ok_loose(func, exonic):
                funnel["drop_function"] += 1
                continue
            funnel["pass_function"] += 1
            clnsig = p[idx["CLNSIG"]]
            if not clinvar_ok_loose(clnsig):
                funnel["drop_clinvar_benign"] += 1
                continue
            funnel["pass_clinvar"] += 1
            af = parse_af(p[idx["GnomAD_genome_AF_ALL"]])
            if not pop_af_ok_loose(af):
                funnel["drop_pop_af"] += 1
                continue
            funnel["pass_pop_af"] += 1
            key = (norm_chr(p[idx["CHROM"]]), p[idx["POS"]], p[idx["REF"]], p[idx["ALT"]])
            gene = ""
            if optional["GeneName"] is not None:
                gene = p[optional["GeneName"]]
            if (not gene or gene in {".", "NA"}) and optional["Gene"] is not None:
                gene = p[optional["Gene"]]
            sites[key] = {
                "chrom": p[idx["CHROM"]],
                "pos": p[idx["POS"]],
                "ref": p[idx["REF"]],
                "alt": p[idx["ALT"]],
                "func": func,
                "exonic_func": exonic,
                "consequence": consequence_class(func, exonic),
                "clnsig": clnsig,
                "gnomad_af": "" if af is None else af,
                "gene": gene,
                "sift": p[optional["SIFT"]] if optional["SIFT"] is not None else "",
                "polyphen_hdiv": p[optional["Polyphen2_HDIV"]] if optional["Polyphen2_HDIV"] is not None else "",
                "polyphen_hvar": p[optional["Polyphen2_HVAR"]] if optional["Polyphen2_HVAR"] is not None else "",
                "cadd": p[optional["CADD"]] if optional["CADD"] is not None else "",
                "revel": p[optional["REVEL"]] if optional["REVEL"] is not None else "",
            }
    print(f"ANNOVAR prefilter sites (func+clinvar+popAF): {len(sites)}", flush=True)
    for k in [
        "annovar_rows",
        "drop_function",
        "pass_function",
        "drop_clinvar_benign",
        "pass_clinvar",
        "drop_pop_af",
        "pass_pop_af",
    ]:
        print(f"  FUNNEL\t{k}\t{funnel[k]}", flush=True)

    # Pass 2: VCF carriers among phenotype samples
    carriers: dict[tuple, set[str]] = defaultdict(set)
    n_vcf_hit = 0
    with gzip.open(VCF, "rt", encoding="utf-8", errors="replace") as f:
        samp_idx = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                vs = line.strip().split("\t")[9:]
                samp_idx = [(j, s) for j, s in enumerate(vs) if s in pheno]
                print(f"VCF phenotype samples indexed: {len(samp_idx)}", flush=True)
                continue
            p = line.rstrip("\n").split("\t")
            if p[6] != "PASS":
                continue
            key = (norm_chr(p[0]), p[1], p[3], p[4])
            if key not in sites:
                continue
            n_vcf_hit += 1
            assert samp_idx is not None
            gts = p[9:]
            for j, s in samp_idx:
                gt = gts[j].split(":")[0]
                if gt in {".", "./.", ".|."}:
                    continue
                alleles = gt.replace("|", "/").split("/")
                if any(a not in {"0", "."} for a in alleles):
                    carriers[key].add(s)
            if n_vcf_hit % 50000 == 0:
                print(f"  ... matched PASS sites in VCF: {n_vcf_hit}", flush=True)

    print(f"VCF matched prefilter sites: {n_vcf_hit}", flush=True)

    # Pass 3: apply control-absent + decide / stratum
    strata = Counter()
    cons_ct = Counter()
    group_carrier_sites = Counter()  # sites with >=1 carrier in group
    kept_rows: list[dict] = []
    case_only_genes: Counter = Counter()
    n_keep = 0
    n_drop_ctrl = 0
    n_no_carrier = 0

    for key, meta in sites.items():
        carr = carriers.get(key, set())
        ctrl_n = sum(1 for s in carr if s in controls)
        ab_n = sum(1 for s in carr if s in cases)
        norm_n = sum(1 for s in carr if s in normals)
        af = meta["gnomad_af"]
        af_f = None if af == "" else float(af)
        d = decide_loose(
            func=meta["func"],
            exonic=meta["exonic_func"],
            clnsig=meta["clnsig"],
            gnomad_af=af_f,
            control_carrier_n=ctrl_n,
            sift=meta["sift"],
            polyphen_hdiv=meta["polyphen_hdiv"],
            polyphen_hvar=meta["polyphen_hvar"],
            cadd=meta["cadd"],
            revel=meta["revel"],
        )
        if ctrl_n > 0:
            n_drop_ctrl += 1
        if not carr:
            n_no_carrier += 1
        if not d.keep:
            continue
        # optional: keep sites with zero cohort carriers? still "candidates" but empty
        n_keep += 1
        st = stratum_label(d)
        strata[st] += 1
        cons_ct[meta["consequence"]] += 1
        if ab_n:
            group_carrier_sites["abnormal"] += 1
        if norm_n:
            group_carrier_sites["normal"] += 1
        if ctrl_n:
            group_carrier_sites["control"] += 1  # should be 0 for kept
        if ab_n and ctrl_n == 0:
            for g in str(meta["gene"]).replace(";", ",").split(","):
                g = g.strip()
                if g and g not in {".", "NA"}:
                    case_only_genes[g] += 1
        sf = d.score
        kept_rows.append(
            {
                "chrom": meta["chrom"],
                "pos": meta["pos"],
                "ref": meta["ref"],
                "alt": meta["alt"],
                "gene": meta["gene"],
                "consequence": meta["consequence"],
                "clnsig": meta["clnsig"],
                "clinvar_tiers": ",".join(str(t) for t in sorted(d.clinvar_tiers)) if d.clinvar_tiers else "",
                "gnomad_af": meta["gnomad_af"],
                "n_carriers": len(carr),
                "n_abnormal": ab_n,
                "n_normal": norm_n,
                "n_control": ctrl_n,
                "stratum": st,
                "score_votes": sf.n_damaging_votes,
                "score_labels": "|".join(sf.labels),
                "sift": meta["sift"],
                "polyphen_hdiv": meta["polyphen_hdiv"],
                "cadd": meta["cadd"],
                "revel": meta["revel"],
            }
        )

    # Prefer sites seen in abnormal for a compact export; full kept is large
    kept_rows.sort(key=lambda r: (-r["n_abnormal"], -r["score_votes"], r["gene"], r["chrom"], r["pos"]))
    out_all = OUT_DIR / "loose_candidates.tsv"
    with open(out_all, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(kept_rows[0].keys()) if kept_rows else ["chrom"], delimiter="\t")
        w.writeheader()
        w.writerows(kept_rows)

    ab_present = [r for r in kept_rows if r["n_abnormal"] > 0]
    out_ab = OUT_DIR / "loose_candidates_in_abnormal.tsv"
    with open(out_ab, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(ab_present[0].keys()) if ab_present else ["chrom"], delimiter="\t")
        w.writeheader()
        w.writerows(ab_present)

    summary_path = OUT_DIR / "loose_filter_summary.tsv"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"n_abnormal\t{len(cases)}\n")
        f.write(f"n_normal\t{len(normals)}\n")
        f.write(f"n_control\t{len(controls)}\n")
        f.write(f"annovar_rows\t{funnel['annovar_rows']}\n")
        f.write(f"pass_function\t{funnel['pass_function']}\n")
        f.write(f"pass_clinvar\t{funnel['pass_clinvar']}\n")
        f.write(f"pass_pop_af_prefilter_sites\t{len(sites)}\n")
        f.write(f"vcf_matched_sites\t{n_vcf_hit}\n")
        f.write(f"drop_present_in_control\t{n_drop_ctrl}\n")
        f.write(f"prefilter_no_cohort_carrier\t{n_no_carrier}\n")
        f.write(f"kept_loose_candidates\t{n_keep}\n")
        f.write(f"kept_with_abnormal_carrier\t{len(ab_present)}\n")
        for k, v in cons_ct.most_common():
            f.write(f"consequence_{k}\t{v}\n")
        for k, v in strata.most_common():
            f.write(f"stratum_{k}\t{v}\n")

    print("RESULT", flush=True)
    print(f"  kept_loose_candidates\t{n_keep}", flush=True)
    print(f"  kept_with_abnormal_carrier\t{len(ab_present)}", flush=True)
    print(f"  drop_present_in_control (of prefilter)\t{n_drop_ctrl}", flush=True)
    print("CONSEQUENCE", dict(cons_ct), flush=True)
    print("STRATA_TOP", strata.most_common(12), flush=True)
    print("TOP_GENES_with_abnormal_and_control_absent", case_only_genes.most_common(20), flush=True)
    print(f"WROTE {out_all}", flush=True)
    print(f"WROTE {out_ab}", flush=True)
    print(f"WROTE {summary_path}", flush=True)


if __name__ == "__main__":
    main()
