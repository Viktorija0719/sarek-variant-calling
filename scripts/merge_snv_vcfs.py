#!/usr/bin/env python3
"""
Merge per-sample annotated small-variant VCFs (SNVs/indels) from multiple callers into
wide per-sample tables with richer interpretation-oriented columns.

What this version adds
- ClinVar fields from INFO (ClinVar_CLNSIG / CLNDN / CLNREVSTAT / ALLELEID / CLNVC)
- Richer VEP fields: gene symbol, gene ID, transcript, MANE, canonical, HGVS, exon/intron,
  variant class, dbSNP/COSMIC IDs, gnomAD max AF/pop, SIFT, PolyPhen, splice scores
- Summary genotype / zygosity / depth / counts / allele balance columns
- Blank interpretation-ready columns (inheritance notes, phenotype match, final classification,
  interpretation comment)
- Keeps per-caller GT/DP/AD/AF/GQ columns and merges identical variants into a single row with
  Called_by listing all callers that support the variant.

Expected layout
- annotation/<caller>/<sample>/...annotated.vcf.gz

Outputs
- TSV per sample
- Excel workbook with one sheet per sample (unless --no-xlsx)
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import re
from collections import OrderedDict
from dataclasses import dataclass
from fnmatch import fnmatch
from typing import Dict, Iterable, List, Optional, Tuple

try:
    from openpyxl import Workbook
except Exception:
    Workbook = None


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
    return {k: vals[i] if i < len(vals) else "." for i, k in enumerate(keys)}


def safe_sheet_name(name: str) -> str:
    name = re.sub(r"[\[\]\:\*\?\/\\]", "_", name)
    return name[:31]




def infer_biological_sample(sample_name: str, caller: str) -> str:
    suffix = f"_{caller}"
    if sample_name.endswith(suffix):
        return sample_name[: -len(suffix)]
    return sample_name


def candidate_sample_dir_names(sample: str, caller: str) -> List[str]:
    names = [sample]
    suffixed = f"{sample}_{caller}"
    if suffixed not in names:
        names.append(suffixed)
    return names


def parse_csq_header(header_line: str) -> List[str]:
    m = re.search(r"Format:\s*([^\">]+)", header_line)
    if not m:
        return []
    return m.group(1).strip().strip('"').split("|")


def parse_csq_entries(csq_value: str) -> List[List[str]]:
    if not csq_value or csq_value == ".":
        return []
    return [entry.split("|") for entry in csq_value.split(",")]


@dataclass
class VcfData:
    csq_fields: List[str]
    records: Dict[Tuple[str, int, str, str], dict]


def read_vcf(path: str) -> VcfData:
    csq_fields: List[str] = []
    records: Dict[Tuple[str, int, str, str], dict] = {}

    for line in read_gz_lines(path):
        if line.startswith("##INFO=<ID=CSQ"):
            csq_fields = parse_csq_header(line)
            continue
        if line.startswith("#"):
            continue

        cols = line.split("\t")
        if len(cols) < 8:
            continue

        chrom = cols[0]
        pos = int(cols[1])
        vid = cols[2]
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
            records[key] = {
                "ID": vid,
                "QUAL": qual,
                "FILTER": flt,
                "INFO": info,
                "FORMAT": fmt,
                "SAMPLE": sample_dict,
            }

    return VcfData(csq_fields=csq_fields, records=records)


def discover_samples(annotation_dir: str, callers: Optional[List[str]] = None) -> List[str]:
    samples = set()
    if not os.path.isdir(annotation_dir):
        return []
    caller_dirs = callers or sorted(os.listdir(annotation_dir))
    for caller in caller_dirs:
        p = os.path.join(annotation_dir, caller)
        if not os.path.isdir(p):
            continue
        for sample_name in os.listdir(p):
            sp = os.path.join(p, sample_name)
            if os.path.isdir(sp):
                samples.add(infer_biological_sample(sample_name, caller))
    return sorted(samples)


def pick_small_variant_vcf(sample: str, caller: str, annotation_dir: str) -> Optional[str]:
    for sample_dir_name in candidate_sample_dir_names(sample, caller):
        base = os.path.join(annotation_dir, caller, sample_dir_name)
        if not os.path.isdir(base):
            continue

        patterns = [
            f"{sample_dir_name}.HG002.{caller}.tp_snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.HG002.{caller}.tp_snpEff_VEP.ann.vcf.gz",
            f"{sample_dir_name}.HG002.{caller}.tp_VEP.ann.vcf.gz",
            f"{sample_dir_name}.HG002.{caller}.tp_snpEff.ann.vcf.gz",

            f"{sample_dir_name}.{caller}.*tp_snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*tp_snpEff_VEP.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*tp_VEP.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*tp_snpEff.ann.vcf.gz",

            f"{sample_dir_name}.{caller}.*varlociraptor_snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*varlociraptor_snpEff_VEP.ann.vcf.gz",
            f"{sample_dir_name}.{caller}*snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.{caller}*snpEff_VEP.ann.vcf.gz",
            f"{sample_dir_name}.{caller}*VEP.ann.vcf.gz",
            f"{sample_dir_name}.{caller}*snpEff.ann.vcf.gz",
            f"{sample_dir_name}.{caller}*.vcf.gz",
        ]

        files = sorted(os.listdir(base))
        for pat in patterns:
            for fn in files:
                if fnmatch(fn, pat):
                    return os.path.join(base, fn)
        for fn in files:
            if fn.endswith('.vcf.gz'):
                return os.path.join(base, fn)
    return None

    patterns = [
        f"{sample}.{caller}.*varlociraptor_snpEff_VEP_snpSift.ann.vcf.gz",
        f"{sample}.{caller}.*varlociraptor_snpEff_VEP.ann.vcf.gz",
        f"{sample}.{caller}*snpEff_VEP_snpSift.ann.vcf.gz",
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
    for fn in files:
        if fn.endswith(".vcf.gz"):
            return os.path.join(base, fn)
    return None


def impact_rank(impact: str) -> int:
    order = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
    return order.get((impact or "").upper(), 0)


def csq_get(entry: Optional[List[str]], csq_fields: List[str], name: str) -> str:
    if not entry or not csq_fields:
        return "."
    try:
        i = csq_fields.index(name)
    except ValueError:
        return "."
    return entry[i] if i < len(entry) and entry[i] != "" else "."


def choose_best_csq(csq_entries: List[List[str]], csq_fields: List[str]) -> Optional[List[str]]:
    if not csq_entries or not csq_fields:
        return None
    best = None
    best_key = None
    for e in csq_entries:
        key = (
            1 if csq_get(e, csq_fields, "MANE_PLUS_CLINICAL") not in ("", ".") else 0,
            1 if csq_get(e, csq_fields, "MANE_SELECT") not in ("", ".") else 0,
            1 if csq_get(e, csq_fields, "MANE") not in ("", ".") else 0,
            1 if csq_get(e, csq_fields, "CANONICAL") == "YES" else 0,
            impact_rank(csq_get(e, csq_fields, "IMPACT")),
        )
        if best is None or key > best_key:
            best = e
            best_key = key
    return best


def resolve_requested_csq_fields(csq_fields: List[str], requested: List[str]) -> List[str]:
    resolved: List[str] = []
    for r in requested:
        r = r.strip()
        if not r:
            continue
        has_wild = any(ch in r for ch in "*?[")
        if has_wild:
            for f in csq_fields:
                if fnmatch(f, r) and f not in resolved:
                    resolved.append(f)
        else:
            if r not in resolved:
                resolved.append(r)
    return resolved


def _to_float(s: str) -> Optional[float]:
    if s in (None, "", "."):
        return None
    try:
        return float(s)
    except Exception:
        return None


def parse_counts_from_sample(sd: Dict[str, str]) -> Tuple[str, str, str, str]:
    dp = sd.get("DP", ".")
    ad = sd.get("AD", ".")
    af = sd.get("AF", ".")

    ref_count = "."
    alt_count = "."
    if ad not in (".", ""):
        parts = ad.split(",")
        if len(parts) >= 2:
            ref_count = parts[0]
            alt_count = parts[1]

    ab = "."
    alt_n = _to_float(alt_count)
    ref_n = _to_float(ref_count)
    af_n = _to_float(af.split(",")[0] if af not in (".", "") else ".")
    if af_n is not None:
        ab = f"{af_n:.4g}"
    elif alt_n is not None and ref_n is not None and (alt_n + ref_n) > 0:
        ab = f"{alt_n/(alt_n+ref_n):.4g}"

    return dp, ref_count, alt_count, ab


def normalize_zygosity(gt: str) -> str:
    if gt in (".", "", None):
        return "."
    g = gt.replace("|", "/")
    if g in ("0/1", "1/0"):
        return "het"
    if g == "1/1":
        return "hom_alt"
    if g == "0/0":
        return "hom_ref"
    if "." in g:
        return "partial_missing"
    alleles = [x for x in g.split("/") if x != "."]
    if not alleles:
        return "."
    uniq = set(alleles)
    if uniq == {"0"}:
        return "hom_ref"
    if len(uniq) == 1 and "0" not in uniq:
        return "hom_alt"
    if "0" in uniq and len(uniq) > 1:
        return "het"
    return "other"


def derive_hgvsg(chrom: str, pos: int, ref: str, alt: str) -> str:
    if ref == alt:
        return f"{chrom}:g.{pos}{ref}="
    if len(ref) == 1 and len(alt) == 1:
        return f"{chrom}:g.{pos}{ref}>{alt}"
    if len(ref) > len(alt) and alt == ref[0]:
        deleted = ref[1:]
        start = pos + 1
        end = pos + len(deleted)
        return f"{chrom}:g.{start}_{end}del{deleted}" if len(deleted) > 1 else f"{chrom}:g.{start}del{deleted}"
    if len(alt) > len(ref) and ref == alt[0]:
        inserted = alt[1:]
        return f"{chrom}:g.{pos}_{pos+1}ins{inserted}"
    return f"{chrom}:g.{pos}_{pos+len(ref)-1}delins{alt}"


def split_existing_variation(existing: str) -> Tuple[str, str, str]:
    if existing in (".", "", None):
        return ".", ".", "."
    items = [x for x in re.split(r"[&,]", existing) if x]
    rs = sorted({x for x in items if x.startswith("rs")})
    cosmic = sorted({x for x in items if x.startswith("COS")})
    other = sorted({x for x in items if x not in rs and x not in cosmic})
    return (";".join(rs) if rs else ".", ";".join(cosmic) if cosmic else ".", ";".join(other) if other else ".")


def spliceai_max_from_entry(entry: Optional[List[str]], csq_fields: List[str]) -> Tuple[str, str, str]:
    if not entry:
        return ".", ".", "."
    ds_names = [
        "SpliceAI_pred_DS_AG",
        "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_DG",
        "SpliceAI_pred_DS_DL",
    ]
    dp_names = [
        "SpliceAI_pred_DP_AG",
        "SpliceAI_pred_DP_AL",
        "SpliceAI_pred_DP_DG",
        "SpliceAI_pred_DP_DL",
    ]
    ds_vals = [(_to_float(csq_get(entry, csq_fields, n)), n) for n in ds_names]
    dp_vals = [(_to_float(csq_get(entry, csq_fields, n)), n) for n in dp_names]
    ds_present = [(v, n) for v, n in ds_vals if v is not None]
    dp_present = [(v, n) for v, n in dp_vals if v is not None]
    ds = max(ds_present, default=(None, "."), key=lambda x: (-1 if x[0] is None else x[0]))
    dp = max(dp_present, default=(None, "."), key=lambda x: (-1 if x[0] is None else x[0]))
    symbol = csq_get(entry, csq_fields, "SpliceAI_pred_SYMBOL")
    return (f"{ds[0]:.4g}" if ds[0] is not None else ".", f"{dp[0]:.4g}" if dp[0] is not None else ".", symbol)


DEFAULT_EXTRA_CSQ = [
    "SYMBOL", "Gene", "Feature", "HGNC_ID",
    "MANE", "MANE_SELECT", "MANE_PLUS_CLINICAL", "CANONICAL",
    "Consequence", "IMPACT", "HGVSc", "HGVSp", "EXON", "INTRON",
    "VARIANT_CLASS", "Existing_variation", "MAX_AF", "MAX_AF_POPS",
    "SIFT", "PolyPhen", "CLIN_SIG", "GENE_PHENO",
    "SpliceAI_pred_SYMBOL", "SpliceAI_pred_DS_*", "SpliceAI_pred_DP_*",
]

DEFAULT_KEEP_INFO = [
    "ClinVar_CLNSIG", "ClinVar_CLNDN", "ClinVar_CLNREVSTAT", "ClinVar_ALLELEID", "ClinVar_CLNVC",
    "PROB*", "HINTS",
]


def build_rows_for_sample(
    sample: str,
    caller_to_vcf_path: Dict[str, Optional[str]],
    extra_csq_requested: List[str],
    keep_info_patterns: List[str],
) -> Tuple[List[str], List[Dict[str, str]]]:
    caller_vcfs: Dict[str, VcfData] = {}
    for caller, path in caller_to_vcf_path.items():
        if path:
            caller_vcfs[caller] = read_vcf(path)

    if not caller_vcfs:
        return [], []

    first = next(iter(caller_vcfs.values()))
    csq_fields = first.csq_fields or []
    resolved_extra_csq = resolve_requested_csq_fields(csq_fields, extra_csq_requested)

    all_keys = set()
    for v in caller_vcfs.values():
        all_keys |= set(v.records.keys())
    all_keys = sorted(all_keys, key=lambda x: (x[0], x[1], x[2], x[3]))

    info_keys_keep = set()
    for vcf in caller_vcfs.values():
        for rec in vcf.records.values():
            info = rec.get("INFO") or {}
            for key in info.keys():
                if any(fnmatch(key, pat) for pat in keep_info_patterns):
                    info_keys_keep.add(key)
    info_keys_keep = sorted(info_keys_keep)

    caller_order = sorted(caller_to_vcf_path.keys())
    rows: List[Dict[str, str]] = []

    for chrom, pos, ref, alt in all_keys:
        called_by: List[str] = []
        filters = set()
        best_qual = None
        per_caller = {}

        for caller, vcf in caller_vcfs.items():
            rec = vcf.records.get((chrom, pos, ref, alt))
            if not rec:
                continue
            per_caller[caller] = rec
            called_by.append(caller)
            filters.add(rec.get("FILTER") or ".")
            q = _to_float(rec.get("QUAL"))
            if q is not None:
                best_qual = q if best_qual is None else max(best_qual, q)

        rep = next(iter(per_caller.values()))
        info = rep.get("INFO") or {}
        csq_entries = parse_csq_entries(info.get("CSQ", "."))
        best_csq = choose_best_csq(csq_entries, csq_fields) if csq_entries and csq_fields else None

        overall_gts = []
        overall_dps = []
        overall_ref_counts = []
        overall_alt_counts = []
        overall_abs = []

        row = OrderedDict()
        row["Sample"] = sample
        row["CHROM"] = chrom
        row["POS"] = pos
        row["REF"] = ref
        row["ALT"] = alt
        row["VariantKey"] = f"{chrom}:{pos}:{ref}:{alt}"
        row["HGVSg"] = derive_hgvsg(chrom, pos, ref, alt)
        row["QUAL_max"] = f"{best_qual:.4g}" if isinstance(best_qual, float) else "."
        row["FILTER_union"] = ";".join(sorted(filters)) if filters else "."
        row["Called_by"] = ",".join(sorted(set(called_by), key=lambda x: caller_order.index(x) if x in caller_order else 999))

        gt_summary_parts = []
        for caller in caller_order:
            rec = per_caller.get(caller)
            row[f"{caller}_SampleID"] = f"{sample}_{caller}"
            if not rec:
                row[f"{caller}_GT"] = "."
                row[f"{caller}_Zygosity"] = "."
                row[f"{caller}_DP"] = "."
                row[f"{caller}_RefCount"] = "."
                row[f"{caller}_AltCount"] = "."
                row[f"{caller}_AB_or_VAF"] = "."
                row[f"{caller}_AF"] = "."
                row[f"{caller}_GQ"] = "."
                continue
            sd = rec.get("SAMPLE") or {}
            gt = sd.get("GT", ".")
            zyg = normalize_zygosity(gt)
            dp, ref_count, alt_count, ab = parse_counts_from_sample(sd)
            af = sd.get("AF", ".")
            gq = sd.get("GQ", sd.get("RGQ", "."))

            row[f"{caller}_GT"] = gt
            row[f"{caller}_Zygosity"] = zyg
            row[f"{caller}_DP"] = dp
            row[f"{caller}_RefCount"] = ref_count
            row[f"{caller}_AltCount"] = alt_count
            row[f"{caller}_AB_or_VAF"] = ab
            row[f"{caller}_AF"] = af
            row[f"{caller}_GQ"] = gq

            if gt not in (".", ""):
                overall_gts.append(gt)
                gt_summary_parts.append(f"{caller}:{gt}")
            if dp not in (".", ""):
                overall_dps.append(_to_float(dp))
            if ref_count not in (".", ""):
                overall_ref_counts.append(_to_float(ref_count))
            if alt_count not in (".", ""):
                overall_alt_counts.append(_to_float(alt_count))
            if ab not in (".", ""):
                overall_abs.append(_to_float(ab))

        uniq_gt = sorted(set(overall_gts))
        row["Genotype"] = uniq_gt[0] if len(uniq_gt) == 1 else (";".join(gt_summary_parts) if gt_summary_parts else ".")
        uniq_zyg = sorted(set(normalize_zygosity(gt) for gt in overall_gts if gt not in (".", "")))
        row["Zygosity"] = uniq_zyg[0] if len(uniq_zyg) == 1 else (";".join(uniq_zyg) if uniq_zyg else ".")
        row["Depth_max"] = f"{max([x for x in overall_dps if x is not None]):.0f}" if any(x is not None for x in overall_dps) else "."
        row["RefCount_max"] = f"{max([x for x in overall_ref_counts if x is not None]):.0f}" if any(x is not None for x in overall_ref_counts) else "."
        row["AltCount_max"] = f"{max([x for x in overall_alt_counts if x is not None]):.0f}" if any(x is not None for x in overall_alt_counts) else "."
        row["AB_or_VAF_max"] = f"{max([x for x in overall_abs if x is not None]):.4g}" if any(x is not None for x in overall_abs) else "."

        # Rich annotation fields
        row["Gene_symbol"] = csq_get(best_csq, csq_fields, "SYMBOL")
        row["Gene_ID"] = csq_get(best_csq, csq_fields, "Gene")
        row["HGNC_ID"] = csq_get(best_csq, csq_fields, "HGNC_ID")
        row["Transcript_ID"] = csq_get(best_csq, csq_fields, "Feature")
        row["MANE_Select"] = csq_get(best_csq, csq_fields, "MANE_SELECT")
        row["MANE_Plus_Clinical"] = csq_get(best_csq, csq_fields, "MANE_PLUS_CLINICAL")
        row["MANE"] = csq_get(best_csq, csq_fields, "MANE")
        row["Canonical"] = csq_get(best_csq, csq_fields, "CANONICAL")
        row["Consequence"] = csq_get(best_csq, csq_fields, "Consequence")
        row["Impact"] = csq_get(best_csq, csq_fields, "IMPACT")
        row["HGVSc"] = csq_get(best_csq, csq_fields, "HGVSc")
        row["HGVSp"] = csq_get(best_csq, csq_fields, "HGVSp")
        row["Exon"] = csq_get(best_csq, csq_fields, "EXON")
        row["Intron"] = csq_get(best_csq, csq_fields, "INTRON")
        row["Variant_class"] = csq_get(best_csq, csq_fields, "VARIANT_CLASS")

        existing = csq_get(best_csq, csq_fields, "Existing_variation")
        dbsnp_ids, cosmic_ids, other_ids = split_existing_variation(existing)
        row["dbSNP_IDs"] = dbsnp_ids
        row["COSMIC_IDs"] = cosmic_ids
        row["Other_existing_IDs"] = other_ids
        row["Existing_variation_raw"] = existing

        row["gnomAD_max_AF"] = csq_get(best_csq, csq_fields, "MAX_AF")
        row["gnomAD_max_pop"] = csq_get(best_csq, csq_fields, "MAX_AF_POPS")
        row["SIFT"] = csq_get(best_csq, csq_fields, "SIFT")
        row["PolyPhen"] = csq_get(best_csq, csq_fields, "PolyPhen")
        row["Gene_phenotype_flag"] = csq_get(best_csq, csq_fields, "GENE_PHENO")

        ds_max, dp_max, splice_symbol = spliceai_max_from_entry(best_csq, csq_fields)
        row["SpliceAI_symbol"] = splice_symbol
        row["SpliceAI_DS_MAX"] = ds_max
        row["SpliceAI_DP_MAX"] = dp_max

        # ClinVar preferring explicit INFO annotations if available
        row["ClinVar_significance"] = info.get("ClinVar_CLNSIG", csq_get(best_csq, csq_fields, "CLIN_SIG"))
        row["ClinVar_condition"] = info.get("ClinVar_CLNDN", ".")
        row["ClinVar_review_strength"] = info.get("ClinVar_CLNREVSTAT", ".")
        allele_id = info.get("ClinVar_ALLELEID", ".")
        row["ClinVar_ID"] = f"ALLELEID:{allele_id}" if allele_id not in (".", "") else "."
        row["ClinVar_variant_type"] = info.get("ClinVar_CLNVC", ".")

        # Blank/placeholder interpretation columns
        row["Inheritance_relevant_notes"] = "."
        row["Phenotype_match"] = "."
        row["Gene_disease_validity"] = "."
        row["Final_classification"] = "."
        row["Interpretation_comment"] = "."

        for k in info_keys_keep:
            row[f"INFO_{k}"] = info.get(k, ".")

        for f in resolved_extra_csq:
            row[f"VEP_{f}"] = csq_get(best_csq, csq_fields, f)

        rows.append(row)

    base_cols = [
        "Sample", "VariantKey", "CHROM", "POS", "REF", "ALT", "HGVSg",
        "Genotype", "Zygosity", "FILTER_union", "QUAL_max",
        "Depth_max", "RefCount_max", "AltCount_max", "AB_or_VAF_max",
        "Called_by",
        "Gene_symbol", "Gene_ID", "HGNC_ID", "Transcript_ID",
        "MANE_Select", "MANE_Plus_Clinical", "MANE", "Canonical",
        "Consequence", "Impact", "HGVSc", "HGVSp", "Exon", "Intron",
        "Variant_class",
        "ClinVar_significance", "ClinVar_condition", "ClinVar_review_strength",
        "ClinVar_ID", "ClinVar_variant_type",
        "dbSNP_IDs", "COSMIC_IDs", "Other_existing_IDs", "Existing_variation_raw",
        "gnomAD_max_AF", "gnomAD_max_pop", "SIFT", "PolyPhen",
        "SpliceAI_symbol", "SpliceAI_DS_MAX", "SpliceAI_DP_MAX",
        "Gene_phenotype_flag",
        "Inheritance_relevant_notes", "Phenotype_match", "Gene_disease_validity",
        "Final_classification", "Interpretation_comment",
    ]

    per_caller_cols = []
    for caller in caller_order:
        per_caller_cols += [
            f"{caller}_GT", f"{caller}_Zygosity", f"{caller}_DP", f"{caller}_RefCount",
            f"{caller}_AltCount", f"{caller}_AB_or_VAF", f"{caller}_AF", f"{caller}_GQ",
        ]

    extra_cols = [f"VEP_{f}" for f in resolved_extra_csq]
    info_cols = [f"INFO_{k}" for k in info_keys_keep]

    cols = base_cols + per_caller_cols + extra_cols + info_cols
    leftovers = sorted({k for r in rows for k in r.keys()} - set(cols))
    cols += leftovers
    return cols, rows


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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotation-dir", default="annotation", help="Path to annotation/ (default: annotation)")
    ap.add_argument("--callers", default="deepvariant,freebayes,strelka,haplotypecaller", help="Comma-separated callers to merge")
    ap.add_argument("--samples", default="", help="Optional comma-separated sample list (default: auto-discover)")
    ap.add_argument("--out-xlsx", default="merged_snv.xlsx", help="Output Excel workbook")
    ap.add_argument("--out-tsv-dir", default="merged_snv", help="Output TSV directory (per sample)")
    ap.add_argument("--no-xlsx", action="store_true", help="Skip Excel, only TSVs")
    ap.add_argument("--extra-csq", default=",".join(DEFAULT_EXTRA_CSQ), help="Comma-separated CSQ field names/patterns to include")
    ap.add_argument("--keep-info", default=",".join(DEFAULT_KEEP_INFO), help="Comma-separated INFO key patterns to include as columns")
    args = ap.parse_args()

    callers = [c.strip() for c in args.callers.split(",") if c.strip()]
    requested_extra = [x.strip() for x in args.extra_csq.split(",") if x.strip()]
    keep_info_patterns = [x.strip() for x in args.keep_info.split(",") if x.strip()]

    if args.samples.strip():
        samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    else:
        samples = discover_samples(args.annotation_dir, callers)

    if not samples:
        raise SystemExit("No samples found under annotation/<caller>/<sample>/")

    per_sample_out: Dict[str, Tuple[List[str], List[Dict[str, str]]]] = {}

    for sample in samples:
        caller_to_vcf = {c: pick_small_variant_vcf(sample, c, args.annotation_dir) for c in callers}
        if all(v is None for v in caller_to_vcf.values()):
            continue

        print(f"[{sample}]")
        for c, v in caller_to_vcf.items():
            print(f"  {c}: {v if v else '(missing)'}")

        cols, rows = build_rows_for_sample(sample, caller_to_vcf, requested_extra, keep_info_patterns)
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
