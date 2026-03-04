#!/usr/bin/env python3
"""
Merge per-sample Sarek SV VCFs (e.g. manta, tiddit) into AnnotSV-ish wide tables.

Updates in this rewrite:
- Includes extra VEP CSQ fields (default: CLIN_SIG + SpliceAI DS/DP fields, when present)
- Adds derived SpliceAI max columns (VEP_SpliceAI_DS_MAX / VEP_SpliceAI_DP_MAX)
- Keeps old behavior: one "full" row per SV + optional transcript "split" rows

Example
  python3 scripts/merge_sv_vcfs.py \
    --annotation-dir results/annotation \
    --callers manta,tiddit \
    --emit-split \
    --out-xlsx merged_sv.xlsx \
    --out-tsv-dir merged_sv
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from fnmatch import fnmatch
from typing import Dict, Iterable, List, Optional, Tuple

try:
    from openpyxl import Workbook
except Exception:
    Workbook = None


# ----------------------------
# Utilities
# ----------------------------

def read_gz_lines(path: str) -> Iterable[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            yield line.rstrip("\n")


def parse_info(info_str: str) -> Dict[str, str]:
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


def safe_sheet_name(name: str) -> str:
    name = re.sub(r"[\[\]\:\*\?\/\\]", "_", name)
    return name[:31]


# ----------------------------
# VEP CSQ parsing
# ----------------------------

def parse_csq_header(header_line: str) -> List[str]:
    m = re.search(r"Format:\s*([^\">]+)", header_line)
    if not m:
        return []
    return m.group(1).strip().strip('"').split("|")


def parse_csq_entries(csq_value: str) -> List[List[str]]:
    if not csq_value or csq_value == ".":
        return []
    return [entry.split("|") for entry in csq_value.split(",")]


def impact_rank(impact: str) -> int:
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


def resolve_requested_csq_fields(csq_fields: List[str], requested: List[str]) -> List[str]:
    resolved: List[str] = []
    for r in requested:
        r = r.strip()
        if not r:
            continue
        has_wild = "*" in r or "?" in r or "[" in r
        if has_wild:
            for f in csq_fields:
                if fnmatch(f, r) and f not in resolved:
                    resolved.append(f)
        else:
            if r not in resolved:
                resolved.append(r)
    return resolved


def spliceai_max_from_entry(entry: List[str], csq_fields: List[str]) -> Tuple[str, str]:
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
class SvRec:
    chrom: str
    pos: int
    _id: str
    ref: str
    alt: str
    qual: str
    flt: str
    info: Dict[str, str]


@dataclass
class VcfSvData:
    csq_fields: List[str]
    records: List[SvRec]


def read_sv_vcf(path: str) -> VcfSvData:
    csq_fields: List[str] = []
    records: List[SvRec] = []

    for line in read_gz_lines(path):
        if line.startswith("##INFO=<ID=CSQ"):
            csq_fields = parse_csq_header(line)
            continue
        if line.startswith("#CHROM"):
            continue
        if line.startswith("#"):
            continue

        cols = line.split("\t")
        if len(cols) < 8:
            continue

        chrom, pos, _id, ref, alt, qual, flt, info_str = cols[:8]
        info = parse_info(info_str)

        records.append(SvRec(
            chrom=chrom,
            pos=int(pos),
            _id=_id,
            ref=ref,
            alt=alt,
            qual=qual,
            flt=flt,
            info=info
        ))

    return VcfSvData(csq_fields=csq_fields, records=records)


# ----------------------------
# Discovery + file picking
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


def pick_sv_vcf(sample: str, caller: str, annotation_dir: str) -> Optional[str]:
    base = os.path.join(annotation_dir, caller, sample)
    if not os.path.isdir(base):
        return None

    patterns = [
        f"{sample}.{caller}.*germline.varlociraptor_snpEff_VEP.ann.vcf.gz",
        f"{sample}.{caller}.*varlociraptor_snpEff_VEP.ann.vcf.gz",
        f"{sample}.{caller}.*diploid_sv_snpEff_VEP.ann.vcf.gz",
        f"{sample}.{caller}.*_snpEff_VEP.ann.vcf.gz",
        f"{sample}.{caller}.*VEP.ann.vcf.gz",
        f"{sample}.{caller}*.vcf.gz",
    ]

    files = sorted(os.listdir(base))
    for pat in patterns:
        for fn in files:
            if fnmatch(fn, pat):
                return os.path.join(base, fn)

    for fn in files:
        if fn.endswith(".vcf.gz"):
            return os.path.join(base, fn)
    return None


# ----------------------------
# Row builders
# ----------------------------

DEFAULT_EXTRA_CSQ = [
    "CLIN_SIG",
    "SpliceAI_pred_SYMBOL",
    "SpliceAI_pred_DS_*",
    "SpliceAI_pred_DP_*",
]

DEFAULT_KEEP_INFO = [
    "SVTYPE",
    "SVLEN",
    "END",
    "CIPOS",
    "CIEND",
    "MATEID",
    "EVENT",
    "HINTS",
    "PROB*",
    "ClinVar*",
    "Mastermind*",
    "MASTER*",
]

BASE_COLS = [
    "Sample",
    "AnnotSV_ID",
    "SV_chrom",
    "SV_start",
    "SV_end",
    "SV_length",
    "SV_type",
    "SVLEN",
    "QUAL",
    "FILTER",
    "called_by",
    "Annotation_mode",   # full / split
    "Gene_name",
    "Transcript",
    "Consequence",
    "IMPACT",
]


def sv_key(rec: SvRec) -> Tuple[str, int, str, str]:
    return (rec.chrom, rec.pos, rec.ref, rec.alt)


def derive_sv_end(info: Dict[str, str], pos: int) -> int:
    for k in ("END", "SVEND", "END2"):
        if k in info:
            try:
                return int(float(info[k]))
            except Exception:
                pass
    # if no END, use POS as fallback
    return pos


def derive_svtype(info: Dict[str, str], alt: str) -> str:
    if "SVTYPE" in info and info["SVTYPE"]:
        return info["SVTYPE"]
    # fallback: symbolic ALT like <DEL>
    m = re.match(r"<([^>]+)>", alt or "")
    return m.group(1) if m else "."


def build_rows_for_sample(
    sample: str,
    caller_to_vcf: Dict[str, Optional[str]],
    extra_csq_requested: List[str],
    keep_info_patterns: List[str],
    emit_split: bool,
) -> Tuple[List[str], List[Dict[str, str]]]:

    # load SV VCFs
    caller_data: Dict[str, VcfSvData] = {}
    for caller, path in caller_to_vcf.items():
        if path:
            caller_data[caller] = read_sv_vcf(path)

    if not caller_data:
        return [], []

    # use CSQ schema from first vcf
    first = next(iter(caller_data.values()))
    csq_fields = first.csq_fields or []
    resolved_extra = resolve_requested_csq_fields(csq_fields, extra_csq_requested)

    # INFO keys to keep (pattern matched across all SV VCFs)
    info_keys_keep = set()
    for _caller, _data in caller_data.items():
        for _rec in _data.records:
            for _k in (_rec.info or {}).keys():
                if any(fnmatch(_k, pat) for pat in keep_info_patterns):
                    info_keys_keep.add(_k)
    info_keys_keep = sorted(info_keys_keep)

    # merge by SV key
    by_key: Dict[Tuple[str, int, str, str], Dict[str, SvRec]] = defaultdict(dict)
    for caller, data in caller_data.items():
        for rec in data.records:
            by_key[sv_key(rec)][caller] = rec

    rows: List[Dict[str, str]] = []

    for key in sorted(by_key.keys(), key=lambda x: (x[0], x[1], x[2], x[3])):
        per = by_key[key]
        callers = sorted(per.keys())
        rep = per[callers[0]]

        chrom, pos, ref, alt = rep.chrom, rep.pos, rep.ref, rep.alt
        info = rep.info
        sv_end = derive_sv_end(info, pos)
        svtype = derive_svtype(info, alt)
        svlen = info.get("SVLEN", ".")
        try:
            sv_length = abs(int(float(svlen))) if svlen not in (".", "", None) else abs(sv_end - pos)
        except Exception:
            sv_length = abs(sv_end - pos)

        ann_id = f"{sample}_{chrom}_{pos}_{sv_end}_{svtype}"

        # CSQ
        csq_entries = parse_csq_entries(info.get("CSQ", "."))
        best = choose_best_csq(csq_entries, csq_fields) if csq_entries and csq_fields else None

        # best summary fields
        gene = csq_get(best, csq_fields, "SYMBOL") if best else "."
        transcript = csq_get(best, csq_fields, "Feature") if best else "."
        consequence = csq_get(best, csq_fields, "Consequence") if best else "."
        impact = csq_get(best, csq_fields, "IMPACT") if best else "."

        # extra csq
        extra_vals = {}
        if best:
            for f in resolved_extra:
                extra_vals[f] = csq_get(best, csq_fields, f)
            ds_max, dp_max = spliceai_max_from_entry(best, csq_fields)
        else:
            for f in resolved_extra:
                extra_vals[f] = "."
            ds_max, dp_max = ".", "."

        # FULL row
        full = {
            "Sample": sample,
            "AnnotSV_ID": ann_id,
            "SV_chrom": chrom,
            "SV_start": pos,
            "SV_end": sv_end,
            "SV_length": sv_length,
            "SV_type": svtype,
            "SVLEN": svlen,
            "QUAL": rep.qual,
            "FILTER": rep.flt,
            "called_by": ",".join(callers),
            "Annotation_mode": "full",
            "Gene_name": gene,
            "Transcript": transcript,
            "Consequence": consequence,
            "IMPACT": impact,
        }
        for f in resolved_extra:
            full[f"VEP_{f}"] = extra_vals.get(f, ".")
        full["VEP_SpliceAI_DS_MAX"] = ds_max
        full["VEP_SpliceAI_DP_MAX"] = dp_max

        # Selected INFO fields
        for _k in info_keys_keep:
            full[f"INFO_{_k}"] = info.get(_k, ".")

        rows.append(full)

        # SPLIT rows: per transcript
        if emit_split and csq_entries and csq_fields:
            idx = {n: i for i, n in enumerate(csq_fields)}
            for e in csq_entries:
                gene2 = csq_get(e, csq_fields, "SYMBOL")
                transcript2 = csq_get(e, csq_fields, "Feature")
                consequence2 = csq_get(e, csq_fields, "Consequence")
                impact2 = csq_get(e, csq_fields, "IMPACT")

                split = dict(full)
                split["Annotation_mode"] = "split"
                split["Gene_name"] = gene2
                split["Transcript"] = transcript2
                split["Consequence"] = consequence2
                split["IMPACT"] = impact2

                # extra values per transcript row
                for f in resolved_extra:
                    split[f"VEP_{f}"] = csq_get(e, csq_fields, f)
                ds_max2, dp_max2 = spliceai_max_from_entry(e, csq_fields)
                split["VEP_SpliceAI_DS_MAX"] = ds_max2
                split["VEP_SpliceAI_DP_MAX"] = dp_max2

                rows.append(split)

    # Column order
    extra_cols = [f"VEP_{f}" for f in resolved_extra] + ["VEP_SpliceAI_DS_MAX", "VEP_SpliceAI_DP_MAX"]
    info_cols = [f"INFO_{k}" for k in info_keys_keep]
    cols = BASE_COLS + extra_cols + info_cols
    # any leftover keys
    leftovers = sorted({k for r in rows for k in r.keys()} - set(cols))
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
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def write_xlsx(path: str, per_sample: Dict[str, Tuple[List[str], List[Dict[str, str]]]]) -> None:
    if Workbook is None:
        raise RuntimeError("openpyxl is not installed. Install it with: pip install openpyxl")

    wb = Workbook(write_only=True)
    for sample, (cols, rows) in per_sample.items():
        ws = wb.create_sheet(title=safe_sheet_name(sample))
        if not rows:
            ws.append(["(no SVs found)"])
            continue
        ws.append(cols)
        for r in rows:
            ws.append([r.get(c, "") for c in cols])
    wb.save(path)


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotation-dir", default="annotation", help="Path to annotation/ (default: annotation)")
    ap.add_argument("--callers", default="manta,tiddit", help="Comma-separated SV callers (default: manta,tiddit)")
    ap.add_argument("--samples", default="", help="Optional comma-separated sample list (default: auto-discover)")
    ap.add_argument("--emit-split", action="store_true", help="Emit transcript split rows (recommended)")
    ap.add_argument("--out-xlsx", default="merged_sv_tables.xlsx", help="Output Excel workbook")
    ap.add_argument("--out-tsv-dir", default="merged_sv_tsv", help="Output TSV directory (per sample)")
    ap.add_argument("--no-xlsx", action="store_true", help="Skip Excel, only TSVs")
    ap.add_argument(
        "--extra-csq",
        default=",".join(DEFAULT_EXTRA_CSQ),
        help="Comma-separated CSQ field names/patterns to include (default includes CLIN_SIG + SpliceAI)",
    )
    ap.add_argument(
        "--keep-info",
        default=",".join(DEFAULT_KEEP_INFO),
        help="Comma-separated INFO key patterns to include as columns (default keeps common SV + optional tracks)",
    )
    args = ap.parse_args()

    annotation_dir = args.annotation_dir
    callers = [c.strip() for c in args.callers.split(",") if c.strip()]
    requested_extra = [x.strip() for x in args.extra_csq.split(",") if x.strip()]
    keep_info_patterns = [x.strip() for x in args.keep_info.split(",") if x.strip()]

    if args.samples.strip():
        samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    else:
        samples = discover_samples(annotation_dir)

    if not samples:
        raise SystemExit("No samples found under annotation/<caller>/<sample>/")

    per_sample_out: Dict[str, Tuple[List[str], List[Dict[str, str]]]] = {}

    for sample in samples:
        caller_to_vcf = {c: pick_sv_vcf(sample, c, annotation_dir) for c in callers}
        if all(v is None for v in caller_to_vcf.values()):
            continue

        print(f"[{sample}]")
        for c, v in caller_to_vcf.items():
            print(f"  {c}: {v if v else '(missing)'}")

        cols, rows = build_rows_for_sample(sample, caller_to_vcf, requested_extra, keep_info_patterns, emit_split=args.emit_split)
        per_sample_out[sample] = (cols, rows)

        tsv_path = os.path.join(args.out_tsv_dir, f"{sample}.sv_merged.tsv")
        write_tsv(tsv_path, cols, rows)
        print(f"  wrote TSV: {tsv_path} ({len(rows)} rows)")

    if not args.no_xlsx:
        write_xlsx(args.out_xlsx, per_sample_out)
        print(f"Wrote Excel: {args.out_xlsx}")

    print("Done.")


if __name__ == "__main__":
    main()