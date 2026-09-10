#!/usr/bin/env python3
"""Scan ANNOVAR for Arnadottir pathogenic genes and exact SSV coordinates."""
from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from snv_loose_filter_rules import (  # noqa: E402
    _parse_float,
    clinvar_ok_loose,
    consequence_class,
    function_ok_loose,
    pop_af_ok_loose,
)

ANNO = Path(r"D:/ONT/figure2/S956.snp.annovar.hg38_multianno.txt.gz")
EXACT = ROOT / "tables" / "external_arnadottir2025" / "causal_ssv_grch38.tsv"
OUT = ROOT / "tables" / "external_arnadottir2025"

T9 = {
    "COL1A1", "CPLANE1", "CSNK2B", "DHCR7", "DNMT1", "FGFR2",
    "FOXP1", "LAMA4", "RANBP2", "ZMIZ1",
}
# Table 4 genes + paper extra constrained pLoF DNMs
CANON = {
    "ANKRD26", "BCL11A", "CACNA1G", "CADM3", "COL1A1", "CPLANE1", "CSNK2B",
    "CYFIP2", "DHCR7", "DMXL2", "DNAH7", "DNMT1", "EEF2", "FBLN1", "FGFR2",
    "FOXP1", "FZD2", "GRIA1", "IMPDH1", "KCNB1", "LAMA4", "RANBP2", "SCN2A",
    "SCN8A", "TCF7L2", "ZMIZ1", "DDX5", "ZMYM4", "HECTD1",
}
ALIASES = {
    "CPLANE1": {"CPLANE1", "C5ORF42", "JBTS17"},
    "CADM3": {"CADM3", "IGSF4B", "NECL1", "SYNCAM3", "TSLL1"},
    "CSNK2B": {"CSNK2B", "CSK2B"},
}
TOKEN_TO_CANON: dict[str, str] = {}
for g in CANON:
    TOKEN_TO_CANON[g.upper()] = g
    for a in ALIASES.get(g, set()):
        TOKEN_TO_CANON[a.upper()] = g


def gene_tokens(raw: str) -> set[str]:
    s = (raw or "").strip()
    if s in {"", ".", "NA", "na"}:
        return set()
    out = set()
    for part in s.replace("|", ";").split(";"):
        part = part.strip()
        if not part:
            continue
        sym = part.split("(")[0].strip()
        if sym:
            out.add(sym.upper())
    return out


def match_canon(raw_name: str, raw_gene: str) -> str | None:
    for tok in gene_tokens(raw_name) | gene_tokens(raw_gene):
        if tok in TOKEN_TO_CANON:
            return TOKEN_TO_CANON[tok]
    return None


def norm_chr(c: str) -> str:
    c = str(c).strip()
    if c.lower().startswith("chr"):
        return "chr" + c[3:]
    return "chr" + c


def load_exact() -> dict[tuple[str, int], list[dict]]:
    by: dict[tuple[str, int], list[dict]] = {}
    with EXACT.open(encoding="utf-8") as f:
        for rec in csv.DictReader(f, delimiter="\t"):
            if not rec.get("chrom") or not rec.get("pos"):
                continue
            key = (norm_chr(rec["chrom"]), int(rec["pos"]))
            by.setdefault(key, []).append(rec)
    return by


def main() -> None:
    exact = load_exact()
    exact_pos = set(exact)
    near_rows = []
    gene_rows = []
    n = 0
    with gzip.open(ANNO, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            n += 1
            p = line.rstrip("\n").split("\t")
            chrom = norm_chr(p[idx["CHROM"]])
            try:
                pos = int(p[idx["POS"]])
            except ValueError:
                continue
            ref, alt = p[idx["REF"]], p[idx["ALT"]]
            # exact / nearby coordinate audit
            for dpos in (0, -1, 1, -2, 2):
                k = (chrom, pos + dpos)
                if k in exact_pos:
                    near_rows.append(
                        {
                            "query_gene": ";".join(x["gene"] for x in exact[k]),
                            "query_hgvs": ";".join(x["hgvs_c"] for x in exact[k]),
                            "query_chrom": k[0],
                            "query_pos": k[1],
                            "query_ref": ";".join(x["ref"] for x in exact[k]),
                            "query_alt": ";".join(x["alt"] for x in exact[k]),
                            "annovar_chrom": chrom,
                            "annovar_pos": pos,
                            "annovar_ref": ref,
                            "annovar_alt": alt,
                            "pos_delta": -dpos,
                            "GeneName": p[idx["GeneName"]],
                            "Gene": p[idx["Gene"]],
                            "Func": p[idx["Func"]],
                            "ExonicFunc": p[idx["ExonicFunc"]],
                            "CLNSIG": p[idx["CLNSIG"]],
                            "AAChange": p[idx.get("AAChange", idx["ExonicFunc"])],
                            "FILTER": p[idx["FILTER"]],
                        }
                    )
                    break
            canon = match_canon(p[idx["GeneName"]], p[idx["Gene"]])
            if not canon:
                continue
            func, exonic = p[idx["Func"]], p[idx["ExonicFunc"]]
            cons = consequence_class(func, exonic)
            if not function_ok_loose(func, exonic):
                continue
            af = _parse_float(p[idx["GnomAD_genome_AF_ALL"]])
            gene_rows.append(
                {
                    "canon_gene": canon,
                    "euploid_likely_causal": "yes" if canon in T9 else "no",
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "GeneName": p[idx["GeneName"]],
                    "Func": func,
                    "ExonicFunc": exonic,
                    "consequence": cons,
                    "clnsig": p[idx["CLNSIG"]],
                    "clinvar_ok": "yes" if clinvar_ok_loose(p[idx["CLNSIG"]]) else "no",
                    "gnomad_af": "" if af is None else af,
                    "pop_af_lt1pct": "yes" if pop_af_ok_loose(af) else "no",
                    "AAChange": p[idx["AAChange"]],
                    "FILTER": p[idx["FILTER"]],
                    "sift": p[idx["SIFT"]],
                    "cadd": p[idx["CADD"]],
                    "revel": p[idx["REVEL"]],
                }
            )
            if n % 5_000_000 == 0:
                print(f"... annovar {n} gene_rows {len(gene_rows)}", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    gpath = OUT / "arnadottir_gene_annovar_sites.tsv"
    npath = OUT / "arnadottir_exact_ssv_annovar_nearby.tsv"
    fields_g = list(gene_rows[0].keys()) if gene_rows else ["canon_gene"]
    with gpath.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields_g, delimiter="\t")
        w.writeheader()
        w.writerows(gene_rows)
    fields_n = list(near_rows[0].keys()) if near_rows else ["query_gene"]
    with npath.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields_n, delimiter="\t")
        w.writeheader()
        w.writerows(near_rows)
    print("ANNOVAR_ROWS", n)
    print("GENE_SITES", len(gene_rows), "WROTE", gpath)
    print("EXACT_NEARBY", len(near_rows), "WROTE", npath)


if __name__ == "__main__":
    main()
