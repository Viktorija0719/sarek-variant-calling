#!/usr/bin/env python3
"""
Merge per-sample Sarek small-variant VCFs (SNVs/indels) from multiple callers into a per-sample table.

What this script does
- Looks for annotated VCFs under: annotation/<caller>/<sample>/
- Reads each VCF, merges records by (CHROM, POS, REF, ALT)
- Emits one row per unique variant, adding:
    * which callers called it (Called_by)
    * per-caller GT/DP/AD/AF/GQ (if present)
    * VarLociraptor probabilities (INFO_PROB_* when present)
    * VEP “best transcript” summary (Gene/Consequence/IMPACT/Transcript/HGVSc/HGVSp)
    * Extra VEP CSQ fields (by default: CLIN_SIG + SpliceAI DS/DP fields), prefixed with VEP_

Outputs
- TSV per sample (in --out-tsv-dir)
- One Excel workbook with one sheet per sample (unless --no-xlsx)

Example
  python3 scripts/merge_snv_vcfs.py \
    --annotation-dir results/annotation \
    --callers deepvariant,freebayes,strelka \
    --out-xlsx merged_snv.xlsx \
    --out-tsv-dir merged_snv

Notes
- The “extra CSQ fields” are taken from VEP CSQ FORMAT header.
  If a field isn't present in a given VCF, the column is still created but left empty for that row.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from fnmatch import fnmatch
from typing import Dict, Iterable, List, Optional, Tuple

try:
    from openpyxl import Workbook
except Exception:
    Workbook = None


# ----------------------------
# Helpers
# ----------------------------

def read_gz_lines(path: str) -> Iterable[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            yield line.rstrip("\n")


def parse_info_field(info_str: str) -> Dict[str, str]:
    if not info_str or info_str == ".":
        return {}
    out: Dict[str, str] = {}
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = "true"
    return out


def parse_format_sample(fmt: str, sample: str) -> Dict[str, str]:
    if not fmt or fmt == ".":
        return {}
    keys = fmt.split(":")
    vals = sample.split(":") if sample and sample != "." else []
    out = {}
    for i, k in enumerate(keys):
        out[k] = vals[i] if i < len(vals) else "."
    return out


def safe_sheet_name(name: str) -> str:
    name = re.sub(r"[\[\]\:\*\?\/\\]", "_", name)
    return name[:31]


# ----------------------------
# VEP CSQ / snpEff ANN parsing
# ----------------------------

def parse_csq_header(header_line: str) -> List[str]:
    # ##INFO=<ID=CSQ,...,Description="... Format: Allele|Consequence|IMPACT|...">
    m = re.search(r"Format:\s*([^\">]+)", header_line)
    if not m:
        return []
    return m.group(1).strip().strip('"').split("|")


def parse_csq_entries(csq_value: str) -> List[List[str]]:
    if not csq_value or csq_value == ".":
        return []
    return [entry.split("|") for entry in csq_value.split(",")]


def parse_ann_header(header_line: str) -> List[str]:
    # ##INFO=<ID=ANN,...,Description="... 'Allele | Annotation | Annotation_Impact | Gene_Name | ...' ">
    m = re.search(r"Allele\s*\|\s*Annotation", header_line)
    if not m:
        return []
    # take everything between first quote after Description= and the last quote
    m2 = re.search(r"Description=\"([^\"]+)\"", header_line)
    if not m2:
        return []
    desc = m2.group(1)
    # heuristic: find the first "Allele |" and split on '|'
    if "Allele" not in desc or "|" not in desc:
        return []
    part = desc[desc.find("Allele"):]
    fields = [x.strip() for x in part.split("|")]
    return fields


def parse_ann_entries(ann_value: str) -> List[List[str]]:
    if not ann_value or ann_value == ".":
        return []
    return [entry.split("|") for entry in ann_value.split(",")]


def impact_rank(impact: str) -> int:
    # VEP IMPACT: HIGH > MODERATE > LOW > MODIFIER
    order = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
    return order.get((impact or "").upper(), 0)


def choose_best_csq(csq_entries: List[List[str]], csq_fields: List[str]) -> Optional[List[str]]:
    if not csq_entries or not csq_fields:
        return None
    idx = {name: i for i, name in enumerate(csq_fields)}
    best = None
    best_key = None
    for e in csq_entries:
        impact = e[idx["IMPACT"]] if "IMPACT" in idx and idx["IMPACT"] < len(e) else ""
        canonical = e[idx["CANONICAL"]] if "CANONICAL" in idx and idx["CANONICAL"] < len(e) else ""
        # Higher IMPACT first, prefer canonical YES
        key = (impact_rank(impact), 1 if canonical == "YES" else 0)
        if best is None or key > best_key:
            best = e
            best_key = key
    return best


def csq_get(entry: List[str], csq_fields: List[str], name: str) -> str:
    if not entry or not csq_fields:
        return "."
    try:
        i = csq_fields.index(name)
    except ValueError:
        return "."
    return entry[i] if i < len(entry) and entry[i] != "" else "."


def resolve_requested_csq_fields(csq_fields: List[str], requested: List[str]) -> Tuple[List[str], List[str]]:
    """
    Returns:
      (resolved_fields, missing_explicit_fields)

    - explicit names (no wildcard) are always kept (even if missing) so columns stay stable.
    - patterns with wildcards (*) are expanded if they match.
    """
    resolved: List[str] = []
    missing: List[str] = []

    for r in requested:
        r = r.strip()
        if not r:
            continue
        has_wild = "*" in r or "?" in r or "[" in r
        if has_wild:
            matches = [f for f in csq_fields if fnmatch(f, r)]
            for m in matches:
                if m not in resolved:
                    resolved.append(m)
        else:
            if r in csq_fields:
                if r not in resolved:
                    resolved.append(r)
            else:
                # keep explicit field as "missing" so the column exists
                if r not in resolved:
                    resolved.append(r)
                missing.append(r)
    return resolved, missing


def spliceai_max_from_entry(entry: List[str], csq_fields: List[str]) -> Tuple[str, str]:
    """
    Compute max DS and max DP across AG/AL/DG/DL if SpliceAI fields exist.
    Returns: (ds_max, dp_max) as strings ('.' if unavailable).
    """
    def fval(name: str) -> Optional[float]:
        v = csq_get(entry, csq_fields, name)
        if v in (".", "", None):
            return None
        try:
            return float(v)
        except Exception:
            return None

    ds_vals = [
        fval("SpliceAI_pred_DS_AG"),
        fval("SpliceAI_pred_DS_AL"),
        fval("SpliceAI_pred_DS_DG"),
        fval("SpliceAI_pred_DS_DL"),
    ]
    dp_vals = [
        fval("SpliceAI_pred_DP_AG"),
        fval("SpliceAI_pred_DP_AL"),
        fval("SpliceAI_pred_DP_DG"),
        fval("SpliceAI_pred_DP_DL"),
    ]
    ds = max([x for x in ds_vals if x is not None], default=None)
    dp = max([x for x in dp_vals if x is not None], default=None)
    return (f"{ds:.4g}" if ds is not None else ".", f"{dp:.4g}" if dp is not None else ".")


# ----------------------------
# VCF reader (minimal)
# ----------------------------

@dataclass
class VcfData:
    csq_fields: List[str]
    ann_fields: List[str]
    records: Dict[Tuple[str, int, str, str], dict]  # key -> record dict


def read_vcf(path: str) -> VcfData:
    csq_fields: List[str] = []
    ann_fields: List[str] = []
    records: Dict[Tuple[str, int, str, str], dict] = {}

    header_cols: List[str] = []
    for line in read_gz_lines(path):
        if line.startswith("##INFO=<ID=CSQ"):
            csq_fields = parse_csq_header(line)
            continue
        if line.startswith("##INFO=<ID=ANN"):
            ann_fields = parse_ann_header(line)
            continue
        if line.startswith("#CHROM"):
            header_cols = line.lstrip("#").split("\t")
            continue
        if line.startswith("#"):
            continue

        cols = line.split("\t")
        if len(cols) < 8:
            continue

        chrom = cols[0]
        pos = int(cols[1])
        _id = cols[2]
        ref = cols[3]
        alts = cols[4].split(",") if cols[4] and cols[4] != "." else ["."]
        qual = cols[5]
        flt = cols[6]
        info = parse_info_field(cols[7])

        fmt = cols[8] if len(cols) > 8 else ""
        sample = cols[9] if len(cols) > 9 else ""
        sample_dict = parse_format_sample(fmt, sample)

        for alt in alts:
            key = (chrom, pos, ref, alt)
            # store per-ALT record
            records[key] = {
                "ID": _id,
                "QUAL": qual,
                "FILTER": flt,
                "INFO": info,
                "FORMAT": fmt,
                "SAMPLE": sample_dict,
            }

    return VcfData(csq_fields=csq_fields, ann_fields=ann_fields, records=records)


# ----------------------------
# Input discovery
# ----------------------------

def discover_samples(annotation_dir: str) -> List[str]:
    samples = set()
    if not os.path.isdir(annotation_dir):
        return []
    for caller in os.listdir(annotation_dir):
        p = os.path.join(annotation_dir, caller)
        if not os.path.isdir(p):
            continue
        for sample in os.listdir(p):
            sp = os.path.join(p, sample)
            if os.path.isdir(sp):
                samples.add(sample)
    return sorted(samples)


def pick_small_variant_vcf(sample: str, caller: str, annotation_dir: str) -> Optional[str]:
    """
    Prefer the most annotated VCF in the sample/caller dir:
      * *varlociraptor_snpEff_VEP.ann.vcf.gz
      * *snpEff_VEP.ann.vcf.gz
      * *VEP.ann.vcf.gz
      * *snpEff.ann.vcf.gz
      * *.vcf.gz
    """
    base = os.path.join(annotation_dir, caller, sample)
    if not os.path.isdir(base):
        return None

    patterns = [
        f"{sample}.{caller}.*varlociraptor_snpEff_VEP.ann.vcf.gz",
        f"{sample}.{caller}*varlociraptor_snpEff_VEP.ann.vcf.gz",
        f"{sample}.{caller}*snpEff_VEP.ann.vcf.gz",
        f"{sample}.{caller}*VEP.ann.vcf.gz",
        f"{sample}.{caller}*snpEff.ann.vcf.gz",
        f"{sample}.{caller}*.vcf.gz",
    ]

    files = sorted(os.listdir(base))
    for pat in patterns:
        for fn in files:
            if fnmatch(fn, pat):
                return os.path.join(base, fn)
    # last resort: any vcf.gz
    for fn in files:
        if fn.endswith(".vcf.gz"):
            return os.path.join(base, fn)
    return None


# ----------------------------
# Merging + row construction
# ----------------------------

DEFAULT_EXTRA_CSQ = [
    "CLIN_SIG",
    "SpliceAI_pred_SYMBOL",
    "SpliceAI_pred_DS_*",
    "SpliceAI_pred_DP_*",
    # Examples for future plugins (only if those fields exist in CSQ):
    # "MM_*", "Mastermind_*", "dbNSFP_*"
]

# INFO keys to keep as columns (pattern matched; avoids dumping huge fields like CSQ/ANN)
DEFAULT_KEEP_INFO = [
    "PROB*",
    "HINTS",
    # Examples for optional custom tracks:
    "ClinVar*",
    "MASTER*",
    "Mastermind*",
    "SpliceAI*",
    "gnomAD*",
    "CADD*",
    "dbNSFP*",
]


def build_rows_for_sample(
    sample: str,
    caller_to_vcf_path: Dict[str, Optional[str]],
    extra_csq_requested: List[str],
    keep_info_patterns: List[str],
) -> Tuple[List[str], List[Dict[str, str]]]:

    # Read caller VCFs
    caller_vcfs: Dict[str, VcfData] = {}
    for caller, path in caller_to_vcf_path.items():
        if not path:
            continue
        caller_vcfs[caller] = read_vcf(path)

    if not caller_vcfs:
        return [], []

    # Choose CSQ schema from first available caller (Sarek should be consistent across callers)
    first_v = next(iter(caller_vcfs.values()))
    csq_fields = first_v.csq_fields or []
    # Resolve requested CSQ fields against schema
    resolved_extra_csq, _missing = resolve_requested_csq_fields(csq_fields, extra_csq_requested)

    # Build union of variant keys
    all_keys = set()
    for v in caller_vcfs.values():
        all_keys |= set(v.records.keys())
    all_keys = sorted(all_keys, key=lambda x: (x[0], x[1], x[2], x[3]))

    # INFO keys to keep (pattern matched across all caller VCFs)
    info_keys_keep = set()
    for _caller, _vcf in caller_vcfs.items():
        for _rec in _vcf.records.values():
            _info = _rec.get("INFO") or {}
            for _k in _info.keys():
                if any(fnmatch(_k, pat) for pat in keep_info_patterns):
                    info_keys_keep.add(_k)
    info_keys_keep = sorted(info_keys_keep)

    rows: List[Dict[str, str]] = []

    for (chrom, pos, ref, alt) in all_keys:
        called_by = []
        filters = set()
        best_qual = None

        per_caller = {}  # caller -> record dict
        for caller, vcf in caller_vcfs.items():
            rec = vcf.records.get((chrom, pos, ref, alt))
            if not rec:
                continue
            per_caller[caller] = rec
            called_by.append(caller)
            filters.add(rec.get("FILTER") or ".")
            try:
                q = float(rec.get("QUAL")) if rec.get("QUAL") not in (None, ".", "") else None
            except Exception:
                q = None
            if q is not None:
                best_qual = q if best_qual is None else max(best_qual, q)

        # Pick one record as "representative" for INFO/CSQ
        rep = next(iter(per_caller.values()))
        info = rep.get("INFO") or {}

        # VEP CSQ summary
        csq_entries = parse_csq_entries(info.get("CSQ", "."))
        best_csq = choose_best_csq(csq_entries, csq_fields) if csq_entries and csq_fields else None

        gene = "."
        consequence = "."
        impact = "."
        transcript = "."
        hgvsc = "."
        hgvsp = "."

        if best_csq:
            gene = csq_get(best_csq, csq_fields, "SYMBOL")
            gene = gene if gene != "." else (csq_get(best_csq, csq_fields, "Gene") if "Gene" in csq_fields else ".")
            consequence = csq_get(best_csq, csq_fields, "Consequence")
            impact = csq_get(best_csq, csq_fields, "IMPACT")
            transcript = csq_get(best_csq, csq_fields, "Feature")
            hgvsc = csq_get(best_csq, csq_fields, "HGVSc")
            hgvsp = csq_get(best_csq, csq_fields, "HGVSp")

        # Extra CSQ fields (including SpliceAI / ClinVar if present)
        extra_vals = {}
        if best_csq:
            for f in resolved_extra_csq:
                extra_vals[f] = csq_get(best_csq, csq_fields, f)
            # derived SpliceAI max
            ds_max, dp_max = spliceai_max_from_entry(best_csq, csq_fields)
        else:
            for f in resolved_extra_csq:
                extra_vals[f] = "."
            ds_max, dp_max = ".", "."

        # Selected INFO fields (pattern matched; avoids dumping CSQ/ANN)
        kept_info = {k: info.get(k, ".") for k in info_keys_keep}

        row = OrderedDict()
        row["Sample"] = sample
        row["VariantKey"] = f"{chrom}:{pos}:{ref}:{alt}"
        row["CHROM"] = chrom
        row["POS"] = pos
        row["REF"] = ref
        row["ALT"] = alt
        row["QUAL_max"] = f"{best_qual:.4g}" if isinstance(best_qual, float) else (best_qual if best_qual is not None else ".")
        row["FILTER_union"] = ";".join(sorted(filters)) if filters else "PASS"
        row["Called_by"] = ",".join(sorted(set(called_by))) if called_by else "."

        row["Gene"] = gene
        row["Consequence"] = consequence
        row["IMPACT"] = impact
        row["Transcript"] = transcript
        row["HGVSc"] = hgvsc
        row["HGVSp"] = hgvsp

        # Extra CSQ fields as VEP_* columns
        for f in resolved_extra_csq:
            row[f"VEP_{f}"] = extra_vals.get(f, ".")
        # Derived SpliceAI summaries
        row["VEP_SpliceAI_DS_MAX"] = ds_max
        row["VEP_SpliceAI_DP_MAX"] = dp_max

        # Per-caller columns (fixed set)
        for caller in sorted(caller_to_vcf_path.keys()):
            rec = per_caller.get(caller)
            if not rec:
                row[f"{caller}_GT"] = "."
                row[f"{caller}_DP"] = "."
                row[f"{caller}_AD"] = "."
                row[f"{caller}_AF"] = "."
                row[f"{caller}_GQ"] = "."
                continue
            sd = rec.get("SAMPLE") or {}
            row[f"{caller}_GT"] = sd.get("GT", ".")
            row[f"{caller}_DP"] = sd.get("DP", ".")
            row[f"{caller}_AD"] = sd.get("AD", ".")
            row[f"{caller}_AF"] = sd.get("AF", ".")
            # some callers use GQ, others use RGQ, or no GQ
            row[f"{caller}_GQ"] = sd.get("GQ", sd.get("RGQ", "."))

        # INFO columns (stable)
        for k in info_keys_keep:
            row[f"INFO_{k}"] = kept_info.get(k, ".")

        rows.append(row)

    # Column order (stable)
    base_cols = [
        "Sample", "VariantKey", "CHROM", "POS", "REF", "ALT",
        "QUAL_max", "FILTER_union", "Called_by",
        "Gene", "Consequence", "IMPACT", "Transcript", "HGVSc", "HGVSp",
    ]
    extra_cols = [f"VEP_{f}" for f in resolved_extra_csq] + ["VEP_SpliceAI_DS_MAX", "VEP_SpliceAI_DP_MAX"]
    per_caller_cols = []
    for caller in sorted(caller_to_vcf_path.keys()):
        per_caller_cols += [f"{caller}_GT", f"{caller}_DP", f"{caller}_AD", f"{caller}_AF", f"{caller}_GQ"]
    info_cols = sorted({k for r in rows for k in r.keys() if k.startswith("INFO_")})

    cols = base_cols + extra_cols + per_caller_cols + info_cols
    # add any leftovers (shouldn't happen, but keep safe)
    leftovers = []
    for r in rows:
        for k in r.keys():
            if k not in cols and k not in leftovers:
                leftovers.append(k)
    cols += leftovers

    return cols, rows


# ----------------------------
# Writers
# ----------------------------

def write_tsv(path: str, cols: List[str], rows: List[Dict[str, str]]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, ".")) for c in cols) + "\n")


def write_excel(path: str, per_sample: Dict[str, Tuple[List[str], List[Dict[str, str]]]]) -> None:
    if Workbook is None:
        raise RuntimeError("openpyxl is not installed. Install it with: pip install openpyxl")

    wb = Workbook(write_only=True)

    for sample, (cols, rows) in per_sample.items():
        ws = wb.create_sheet(title=safe_sheet_name(sample))
        if not rows:
            ws.append(["(no variants found)"])
            continue
        ws.append(cols)
        for r in rows:
            ws.append([r.get(c, ".") for c in cols])

    wb.save(path)


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotation-dir", default="annotation", help="Path to annotation/ (default: annotation)")
    ap.add_argument("--callers", default="deepvariant,freebayes,strelka", help="Comma-separated callers to merge")
    ap.add_argument("--samples", default="", help="Optional comma-separated sample list (default: auto-discover)")
    ap.add_argument("--out-xlsx", default="merged_small_variant_tables.xlsx", help="Output Excel workbook")
    ap.add_argument("--out-tsv-dir", default="merged_small_variant_tsv", help="Output TSV directory (per sample)")
    ap.add_argument("--no-xlsx", action="store_true", help="Skip Excel, only TSVs")
    ap.add_argument(
        "--extra-csq",
        default=",".join(DEFAULT_EXTRA_CSQ),
        help="Comma-separated CSQ field names/patterns to include (default includes CLIN_SIG + SpliceAI)",
    )
    ap.add_argument(
        "--keep-info",
        default=",".join(DEFAULT_KEEP_INFO),
        help="Comma-separated INFO key patterns to include as columns (default keeps PROB*/HINTS plus common optional tracks)",
    )
    args = ap.parse_args()

    annotation_dir = args.annotation_dir
    callers = [c.strip() for c in args.callers.split(",") if c.strip()]
    requested_extra_csq = [x.strip() for x in args.extra_csq.split(",") if x.strip()]
    keep_info_patterns = [x.strip() for x in args.keep_info.split(",") if x.strip()]

    if args.samples.strip():
        samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    else:
        samples = discover_samples(annotation_dir)

    if not samples:
        raise SystemExit("No samples found under annotation/<caller>/<sample>/")

    per_sample_out: Dict[str, Tuple[List[str], List[Dict[str, str]]]] = {}

    for sample in samples:
        caller_to_vcf = {c: pick_small_variant_vcf(sample, c, annotation_dir) for c in callers}
        if all(v is None for v in caller_to_vcf.values()):
            continue

        print(f"[{sample}]")
        for c, v in caller_to_vcf.items():
            print(f"  {c}: {v if v else '(missing)'}")

        cols, rows = build_rows_for_sample(sample, caller_to_vcf, requested_extra_csq, keep_info_patterns)
        per_sample_out[sample] = (cols, rows)

        tsv_path = os.path.join(args.out_tsv_dir, f"{sample}.small_variants_merged.tsv")
        write_tsv(tsv_path, cols, rows)
        print(f"  wrote TSV: {tsv_path} ({len(rows)} rows)")

    if not args.no_xlsx:
        write_excel(args.out_xlsx, per_sample_out)
        print(f"Wrote Excel: {args.out_xlsx}")

    print("Done.")


if __name__ == "__main__":
    main()
