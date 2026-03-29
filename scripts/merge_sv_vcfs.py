#!/usr/bin/env python3
"""
Merge per-sample annotated structural-variant VCFs into wide per-sample tables with richer
interpretation-oriented columns.

What this version adds
- ClinVar-related INFO fields when present
- Richer VEP fields: gene symbol, gene ID, transcript, MANE, canonical, consequence, impact,
  exon/intron, HGVS, variant class, IDs, gnomAD AF/pop, SIFT/PolyPhen, splice scores
- FORMAT-derived genotype / zygosity / depth summaries when available
- Blank interpretation-ready columns
- Keeps old behavior: one full row per SV + optional transcript split rows
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
    fmt: str
    sample_dict: Dict[str, str]


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
        if line.startswith("#"):
            continue

        cols = line.split("\t")
        if len(cols) < 8:
            continue

        chrom, pos, _id, ref, alt, qual, flt, info_str = cols[:8]
        info = parse_info(info_str)
        fmt = cols[8] if len(cols) > 8 else ""
        sample = cols[9] if len(cols) > 9 else ""
        sample_dict = parse_format_sample(fmt, sample)

        records.append(SvRec(
            chrom=chrom,
            pos=int(pos),
            _id=_id,
            ref=ref,
            alt=alt,
            qual=qual,
            flt=flt,
            info=info,
            fmt=fmt,
            sample_dict=sample_dict,
        ))

    return VcfSvData(csq_fields=csq_fields, records=records)


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


def pick_sv_vcf(sample: str, caller: str, annotation_dir: str) -> Optional[str]:
    for sample_dir_name in candidate_sample_dir_names(sample, caller):
        base = os.path.join(annotation_dir, caller, sample_dir_name)
        if not os.path.isdir(base):
            continue
        patterns = [
            f"{sample_dir_name}.HG002.{caller}.tp-comp_snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.HG002.{caller}.tp-comp_snpEff_VEP.ann.vcf.gz",
            f"{sample_dir_name}.HG002.{caller}.tp_snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.HG002.{caller}.tp_snpEff_VEP.ann.vcf.gz",

            f"{sample_dir_name}.{caller}.*tp-comp_snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*tp-comp_snpEff_VEP.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*tp_snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*tp_snpEff_VEP.ann.vcf.gz",

            f"{sample_dir_name}.{caller}.*germline.varlociraptor_snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*germline.varlociraptor_snpEff_VEP.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*varlociraptor_snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*varlociraptor_snpEff_VEP.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*diploid_sv_snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*diploid_sv_snpEff_VEP.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*_snpEff_VEP_snpSift.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*_snpEff_VEP.ann.vcf.gz",
            f"{sample_dir_name}.{caller}.*VEP.ann.vcf.gz",
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


def sv_key(rec: SvRec) -> Tuple[str, int, str, str]:
    return (rec.chrom, rec.pos, rec.ref, rec.alt)


def derive_sv_end(info: Dict[str, str], pos: int) -> int:
    for k in ("END", "SVEND", "END2"):
        if k in info:
            try:
                return int(float(info[k]))
            except Exception:
                pass
    return pos


def derive_svtype(info: Dict[str, str], alt: str) -> str:
    if "SVTYPE" in info and info["SVTYPE"]:
        return info["SVTYPE"]
    m = re.match(r"<([^>]+)>", alt or "")
    return m.group(1) if m else "."


def derive_hgvsg_sv(chrom: str, start: int, end: int, svtype: str) -> str:
    if start == end:
        return f"{chrom}:g.{start}{svtype.lower()}"
    return f"{chrom}:g.{start}_{end}{svtype.lower()}"


DEFAULT_EXTRA_CSQ = [
    "SYMBOL", "Gene", "Feature", "HGNC_ID",
    "MANE", "MANE_SELECT", "MANE_PLUS_CLINICAL", "CANONICAL",
    "Consequence", "IMPACT", "HGVSc", "HGVSp", "EXON", "INTRON",
    "VARIANT_CLASS", "Existing_variation", "MAX_AF", "MAX_AF_POPS",
    "SIFT", "PolyPhen", "CLIN_SIG", "GENE_PHENO",
    "SpliceAI_pred_SYMBOL", "SpliceAI_pred_DS_*", "SpliceAI_pred_DP_*",
]

DEFAULT_KEEP_INFO = [
    "SVTYPE", "SVLEN", "END", "CIPOS", "CIEND", "MATEID", "EVENT",
    "HINTS", "PROB*",
    "ClinVar_CLNSIG", "ClinVar_CLNDN", "ClinVar_CLNREVSTAT", "ClinVar_ALLELEID", "ClinVar_CLNVC",
]

CALLER_ORDER = {"manta": 0, "tiddit": 1}


def build_rows_for_sample(
    sample: str,
    caller_to_vcf: Dict[str, Optional[str]],
    extra_csq_requested: List[str],
    keep_info_patterns: List[str],
    emit_split: bool,
) -> Tuple[List[str], List[Dict[str, str]]]:
    caller_data: Dict[str, VcfSvData] = {}
    for caller, path in caller_to_vcf.items():
        if path:
            caller_data[caller] = read_sv_vcf(path)

    if not caller_data:
        return [], []

    first = next(iter(caller_data.values()))
    csq_fields = first.csq_fields or []
    resolved_extra = resolve_requested_csq_fields(csq_fields, extra_csq_requested)

    info_keys_keep = set()
    for data in caller_data.values():
        for rec in data.records:
            for key in (rec.info or {}).keys():
                if any(fnmatch(key, pat) for pat in keep_info_patterns):
                    info_keys_keep.add(key)
    info_keys_keep = sorted(info_keys_keep)

    by_key: Dict[Tuple[str, int, str, str], Dict[str, SvRec]] = defaultdict(dict)
    for caller, data in caller_data.items():
        for rec in data.records:
            by_key[sv_key(rec)][caller] = rec

    caller_order = sorted(caller_to_vcf.keys(), key=lambda x: CALLER_ORDER.get(x, 999))
    rows: List[Dict[str, str]] = []

    for key in sorted(by_key.keys(), key=lambda x: (x[0], x[1], x[2], x[3])):
        per = by_key[key]
        callers = sorted(per.keys(), key=lambda x: caller_order.index(x) if x in caller_order else 999)
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
        csq_entries = parse_csq_entries(info.get("CSQ", "."))
        best = choose_best_csq(csq_entries, csq_fields) if csq_entries and csq_fields else None

        # summary genotype/depth across callers
        overall_gts = []
        overall_dps = []
        for caller in caller_order:
            rec = per.get(caller)
            if rec and rec.sample_dict.get("GT", ".") not in (".", ""):
                overall_gts.append(rec.sample_dict.get("GT", "."))
            if rec and rec.sample_dict.get("DP", ".") not in (".", ""):
                overall_dps.append(_to_float(rec.sample_dict.get("DP", ".")))

        uniq_gt = sorted(set(overall_gts))
        gt_summary = uniq_gt[0] if len(uniq_gt) == 1 else (";".join(f"{c}:{per[c].sample_dict.get('GT','.') }" for c in callers) if callers else ".")
        uniq_zyg = sorted(set(normalize_zygosity(gt) for gt in overall_gts if gt not in (".", "")))
        zyg_summary = uniq_zyg[0] if len(uniq_zyg) == 1 else (";".join(uniq_zyg) if uniq_zyg else ".")

        full = {
            "Sample": sample,
            "AnnotSV_ID": ann_id,
            "SV_chrom": chrom,
            "SV_start": pos,
            "SV_end": sv_end,
            "SV_length": sv_length,
            "SV_type": svtype,
            "SVLEN": svlen,
            "REF": ref,
            "ALT": alt,
            "HGVSg": derive_hgvsg_sv(chrom, pos, sv_end, svtype),
            "Genotype": gt_summary,
            "Zygosity": zyg_summary,
            "Depth_max": f"{max([x for x in overall_dps if x is not None]):.0f}" if any(x is not None for x in overall_dps) else ".",
            "QUAL": rep.qual,
            "FILTER": rep.flt,
            "called_by": ",".join(callers),
            "Annotation_mode": "full",
            "Gene_symbol": csq_get(best, csq_fields, "SYMBOL"),
            "Gene_ID": csq_get(best, csq_fields, "Gene"),
            "HGNC_ID": csq_get(best, csq_fields, "HGNC_ID"),
            "Transcript_ID": csq_get(best, csq_fields, "Feature"),
            "MANE_Select": csq_get(best, csq_fields, "MANE_SELECT"),
            "MANE_Plus_Clinical": csq_get(best, csq_fields, "MANE_PLUS_CLINICAL"),
            "MANE": csq_get(best, csq_fields, "MANE"),
            "Canonical": csq_get(best, csq_fields, "CANONICAL"),
            "Consequence": csq_get(best, csq_fields, "Consequence"),
            "Impact": csq_get(best, csq_fields, "IMPACT"),
            "HGVSc": csq_get(best, csq_fields, "HGVSc"),
            "HGVSp": csq_get(best, csq_fields, "HGVSp"),
            "Exon": csq_get(best, csq_fields, "EXON"),
            "Intron": csq_get(best, csq_fields, "INTRON"),
            "Variant_class": csq_get(best, csq_fields, "VARIANT_CLASS"),
        }

        existing = csq_get(best, csq_fields, "Existing_variation")
        dbsnp_ids, cosmic_ids, other_ids = split_existing_variation(existing)
        full["dbSNP_IDs"] = dbsnp_ids
        full["COSMIC_IDs"] = cosmic_ids
        full["Other_existing_IDs"] = other_ids
        full["Existing_variation_raw"] = existing
        full["gnomAD_max_AF"] = csq_get(best, csq_fields, "MAX_AF")
        full["gnomAD_max_pop"] = csq_get(best, csq_fields, "MAX_AF_POPS")
        full["SIFT"] = csq_get(best, csq_fields, "SIFT")
        full["PolyPhen"] = csq_get(best, csq_fields, "PolyPhen")
        full["Gene_phenotype_flag"] = csq_get(best, csq_fields, "GENE_PHENO")

        ds_max, dp_max, splice_symbol = spliceai_max_from_entry(best, csq_fields)
        full["SpliceAI_symbol"] = splice_symbol
        full["SpliceAI_DS_MAX"] = ds_max
        full["SpliceAI_DP_MAX"] = dp_max

        full["ClinVar_significance"] = info.get("ClinVar_CLNSIG", csq_get(best, csq_fields, "CLIN_SIG"))
        full["ClinVar_condition"] = info.get("ClinVar_CLNDN", ".")
        full["ClinVar_review_strength"] = info.get("ClinVar_CLNREVSTAT", ".")
        allele_id = info.get("ClinVar_ALLELEID", ".")
        full["ClinVar_ID"] = f"ALLELEID:{allele_id}" if allele_id not in (".", "") else "."
        full["ClinVar_variant_type"] = info.get("ClinVar_CLNVC", ".")

        full["Inheritance_relevant_notes"] = "."
        full["Phenotype_match"] = "."
        full["Gene_disease_validity"] = "."
        full["Final_classification"] = "."
        full["Interpretation_comment"] = "."

        for caller in caller_order:
            rec = per.get(caller)
            if not rec:
                full[f"{caller}_GT"] = "."
                full[f"{caller}_Zygosity"] = "."
                full[f"{caller}_DP"] = "."
                full[f"{caller}_AD"] = "."
                full[f"{caller}_AF"] = "."
                full[f"{caller}_GQ"] = "."
            else:
                sd = rec.sample_dict or {}
                full[f"{caller}_GT"] = sd.get("GT", ".")
                full[f"{caller}_Zygosity"] = normalize_zygosity(sd.get("GT", "."))
                full[f"{caller}_DP"] = sd.get("DP", ".")
                full[f"{caller}_AD"] = sd.get("AD", ".")
                full[f"{caller}_AF"] = sd.get("AF", ".")
                full[f"{caller}_GQ"] = sd.get("GQ", sd.get("RGQ", "."))

        for _k in info_keys_keep:
            full[f"INFO_{_k}"] = info.get(_k, ".")

        for f in resolved_extra:
            full[f"VEP_{f}"] = csq_get(best, csq_fields, f)

        rows.append(full)

        if emit_split and csq_entries and csq_fields:
            for e in csq_entries:
                split = dict(full)
                split["Annotation_mode"] = "split"
                split["Gene_symbol"] = csq_get(e, csq_fields, "SYMBOL")
                split["Gene_ID"] = csq_get(e, csq_fields, "Gene")
                split["HGNC_ID"] = csq_get(e, csq_fields, "HGNC_ID")
                split["Transcript_ID"] = csq_get(e, csq_fields, "Feature")
                split["MANE_Select"] = csq_get(e, csq_fields, "MANE_SELECT")
                split["MANE_Plus_Clinical"] = csq_get(e, csq_fields, "MANE_PLUS_CLINICAL")
                split["MANE"] = csq_get(e, csq_fields, "MANE")
                split["Canonical"] = csq_get(e, csq_fields, "CANONICAL")
                split["Consequence"] = csq_get(e, csq_fields, "Consequence")
                split["Impact"] = csq_get(e, csq_fields, "IMPACT")
                split["HGVSc"] = csq_get(e, csq_fields, "HGVSc")
                split["HGVSp"] = csq_get(e, csq_fields, "HGVSp")
                split["Exon"] = csq_get(e, csq_fields, "EXON")
                split["Intron"] = csq_get(e, csq_fields, "INTRON")
                split["Variant_class"] = csq_get(e, csq_fields, "VARIANT_CLASS")
                ds2, dp2, sym2 = spliceai_max_from_entry(e, csq_fields)
                split["SpliceAI_symbol"] = sym2
                split["SpliceAI_DS_MAX"] = ds2
                split["SpliceAI_DP_MAX"] = dp2
                for f in resolved_extra:
                    split[f"VEP_{f}"] = csq_get(e, csq_fields, f)
                rows.append(split)

    base_cols = [
        "Sample", "AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type", "SVLEN",
        "REF", "ALT", "HGVSg", "Genotype", "Zygosity", "Depth_max",
        "QUAL", "FILTER", "called_by", "Annotation_mode",
        "Gene_symbol", "Gene_ID", "HGNC_ID", "Transcript_ID",
        "MANE_Select", "MANE_Plus_Clinical", "MANE", "Canonical",
        "Consequence", "Impact", "HGVSc", "HGVSp", "Exon", "Intron", "Variant_class",
        "ClinVar_significance", "ClinVar_condition", "ClinVar_review_strength",
        "ClinVar_ID", "ClinVar_variant_type",
        "dbSNP_IDs", "COSMIC_IDs", "Other_existing_IDs", "Existing_variation_raw",
        "gnomAD_max_AF", "gnomAD_max_pop", "SIFT", "PolyPhen",
        "SpliceAI_symbol", "SpliceAI_DS_MAX", "SpliceAI_DP_MAX", "Gene_phenotype_flag",
        "Inheritance_relevant_notes", "Phenotype_match", "Gene_disease_validity",
        "Final_classification", "Interpretation_comment",
    ]

    per_caller_cols = []
    for caller in caller_order:
        per_caller_cols += [f"{caller}_GT", f"{caller}_Zygosity", f"{caller}_DP", f"{caller}_AD", f"{caller}_AF", f"{caller}_GQ"]

    extra_cols = [f"VEP_{f}" for f in resolved_extra]
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotation-dir", default="annotation", help="Path to annotation/ (default: annotation)")
    ap.add_argument("--callers", default="manta,tiddit", help="Comma-separated SV callers")
    ap.add_argument("--samples", default="", help="Optional comma-separated sample list (default: auto-discover)")
    ap.add_argument("--emit-split", action="store_true", help="Emit transcript split rows")
    ap.add_argument("--out-xlsx", default="merged_sv.xlsx", help="Output Excel workbook")
    ap.add_argument("--out-tsv-dir", default="merged_sv", help="Output TSV directory (per sample)")
    ap.add_argument("--no-xlsx", action="store_true", help="Skip Excel, only TSVs")
    ap.add_argument("--extra-csq", default=",".join(DEFAULT_EXTRA_CSQ), help="Comma-separated CSQ field names/patterns to include")
    ap.add_argument("--keep-info", default=",".join(DEFAULT_KEEP_INFO), help="Comma-separated INFO key patterns to include")
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
        caller_to_vcf = {c: pick_sv_vcf(sample, c, args.annotation_dir) for c in callers}
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
