#!/usr/bin/env python3
"""Map Arnadottir Table 4 / Table 9 MANE HGVS to GRCh38 VCF alleles via Ensembl."""
from __future__ import annotations

import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(r"D:/ONT/analyses/burden/tables/external_arnadottir2025")
OUT = ROOT / "causal_ssv_grch38.tsv"

T9_GENES = {
    "COL1A1",
    "CPLANE1",
    "CSNK2B",
    "DHCR7",
    "DNMT1",
    "FGFR2",
    "FOXP1",
    "LAMA4",
    "RANBP2",
    "ZMIZ1",
}

HGVS_FALLBACK = {
    "NM_001360.2:c.964-1G>C": ["NM_001360.3:c.964-1G>C"],
    "NM_001174116.2:c.5763dup": ["NM_001174116.3:c.5763dup", "NM_015263.4:c.5763dup"],
    "NM_018897.2:c.6435del": ["NM_018897.4:c.6435del", "NM_018897.3:c.6435del"],
    "NM_018897.2:c.5656C>T": ["NM_018897.4:c.5656C>T", "NM_018897.3:c.5656C>T"],
}


def nc_to_chr(acc: str) -> str | None:
    m = re.match(r"NC_0000(\d+)\.\d+$", acc)
    if not m:
        return None
    n = int(m.group(1))
    if n == 23:
        return "chrX"
    if n == 24:
        return "chrY"
    return f"chr{n}"


def parse_hgvs_field(text: str) -> list[str]:
    text = (text or "").strip()
    if not text or text == ".":
        return []
    parts: list[str] = []
    last_nm = None
    for chunk in re.split(r";\s*", text):
        compact = chunk.replace(" ", "")
        m = re.search(r"(NM_\d+\.\d+):c\.([^/\s]+)", compact)
        if m:
            last_nm = m.group(1)
            parts.append(f"{m.group(1)}:c.{m.group(2)}")
            continue
        m = re.search(r"^c\.([^/\s]+)", compact)
        if m and last_nm:
            parts.append(f"{last_nm}:c.{m.group(1)}")
    return parts


def recode(hgvs: str) -> dict | list:
    q = urllib.parse.quote(hgvs, safe=":")
    url = (
        "https://rest.ensembl.org/variant_recoder/human/"
        + q
        + "?content-type=application/json"
    )
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            body = e.read().decode(errors="replace")[:400]
            if e.code in (429, 500, 502, 503) and attempt < 3:
                time.sleep(2 ** attempt)
                continue
            return {"error": f"HTTP {e.code} {body}"}
        except Exception as e:  # noqa: BLE001
            if attempt < 3:
                time.sleep(2 ** attempt)
                continue
            return {"error": str(e)}
    return {"error": "failed"}


def recode_with_fallback(hgvs: str) -> dict | list:
    payload = recode(hgvs)
    if not (isinstance(payload, dict) and payload.get("error")):
        return payload
    for alt in HGVS_FALLBACK.get(hgvs, []):
        payload2 = recode(alt)
        if not (isinstance(payload2, dict) and payload2.get("error")):
            return payload2
    return payload


def spdi_to_vcf(spdi: str) -> tuple[str, int, str, str] | None:
    if not spdi.startswith("NC_"):
        return None
    chrom_acc, pos0, ref, alt = spdi.split(":")
    chrom = nc_to_chr(chrom_acc)
    if not chrom:
        return None
    return chrom, int(pos0) + 1, ref, alt


def pick_coords(payload) -> dict:
    out = {
        "spdi": "",
        "hgvsg": "",
        "chrom": "",
        "pos": "",
        "ref": "",
        "alt": "",
        "error": "",
    }
    if isinstance(payload, dict) and payload.get("error"):
        out["error"] = payload["error"]
        return out
    if not isinstance(payload, list) or not payload:
        out["error"] = "empty recoder"
        return out
    inner = next(iter(payload[0].values()))
    spdis = inner.get("spdi") or []
    hgvsg = inner.get("hgvsg") or []
    out["spdi"] = ";".join(spdis[:4])
    out["hgvsg"] = ";".join([x for x in hgvsg if x.startswith("NC_")][:3])
    for s in spdis:
        v = spdi_to_vcf(s)
        if v:
            out["chrom"], out["pos"], out["ref"], out["alt"] = v
            return out
    out["error"] = "no GRCh38 SPDI"
    return out


def main() -> None:
    rows_in = []
    with (ROOT / "pathogenic_ssv_unique.tsv").open(encoding="utf-8") as f:
        for rec in csv.DictReader(f, delimiter="\t"):
            gene = rec["Gene"].strip()
            for i, hgvs in enumerate(
                parse_hgvs_field(rec["Effect on protein (MANE transcript)"]), 1
            ):
                rows_in.append(
                    {
                        "gene": gene,
                        "allele_n": i,
                        "hgvs_c": hgvs,
                        "genotype_class": rec.get("Genotype", ""),
                        "karyotype": rec.get("Fetal karyotype", ""),
                        "euploid_likely_causal": "yes" if gene in T9_GENES else "no",
                        "classification": rec.get("Classification", ""),
                    }
                )

    mapped = []
    for rec in rows_in:
        print("MAP", rec["gene"], rec["hgvs_c"], flush=True)
        rec.update(pick_coords(recode_with_fallback(rec["hgvs_c"])))
        mapped.append(rec)
        time.sleep(0.15)

    fields = [
        "gene",
        "allele_n",
        "hgvs_c",
        "genotype_class",
        "karyotype",
        "euploid_likely_causal",
        "classification",
        "chrom",
        "pos",
        "ref",
        "alt",
        "spdi",
        "hgvsg",
        "error",
    ]
    with OUT.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(mapped)
    print("WROTE", OUT, "N", len(mapped), "mapped", sum(1 for x in mapped if x["chrom"]))


if __name__ == "__main__":
    main()
