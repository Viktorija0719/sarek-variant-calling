#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

CALLER_ORDER_SMALL = {
    "deepvariant": 0,
    "freebayes": 1,
    "strelka": 2,
    "haplotypecaller": 3,
}

CALLER_ORDER_STRUCTURAL = {
    "manta": 0,
    "tiddit": 1,
}


def make_unique_id(sample: str, caller: str) -> str:
    return f"{sample}_{caller}"


def find_first_per_sample(pattern: str) -> Dict[str, Path]:
    matches = sorted(PROJECT_DIR.glob(pattern))
    result: Dict[str, Path] = {}

    for path in matches:
        if path.name.endswith(".tbi"):
            continue

        sample = path.parent.name
        if sample not in result:
            result[sample] = path
        else:
            print(
                f"WARN: multiple matches for sample={sample} pattern={pattern}. "
                f"Keeping {result[sample]} and ignoring {path}",
                file=sys.stderr,
            )
    return result


def find_joint_haplotypecaller_vcf() -> Optional[Path]:
    candidates = sorted(
        PROJECT_DIR.glob(
            "results/variant_calling/haplotypecaller/**/joint_germline_recalibrated.vcf.gz"
        )
    )
    if not candidates:
        return None
    if len(candidates) > 1:
        print(
            f"WARN: multiple joint HaplotypeCaller VCFs found. Using {candidates[0]}",
            file=sys.stderr,
        )
    return candidates[0]


def get_vcf_samples_with_bcftools(vcf_path: Path) -> Optional[List[str]]:
    if shutil.which("bcftools") is None:
        return None

    try:
        res = subprocess.run(
            ["bcftools", "query", "-l", str(vcf_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        samples = [line.strip() for line in res.stdout.splitlines() if line.strip()]
        return samples or None
    except subprocess.CalledProcessError as e:
        print(f"WARN: bcftools query failed for {vcf_path}: {e}", file=sys.stderr)
        return None


def rel_or_abs(path: Path, absolute: bool) -> str:
    return str(path.resolve()) if absolute else str(path.relative_to(PROJECT_DIR))


def resolve_haplotypecaller_subsample(
    sample: str,
    joint_sample_set: Set[str],
) -> Optional[str]:
    """
    Resolve the subsample name for a multisample joint HaplotypeCaller VCF.

    Matching order:
    1. exact match: sample
    2. duplicated pattern: sample_sample
    3. otherwise None
    """
    if sample in joint_sample_set:
        return sample

    duplicated = f"{sample}_{sample}"
    if duplicated in joint_sample_set:
        return duplicated

    return None


def write_small_samplesheet(
    out_csv: Path,
    deepvariant: Dict[str, Path],
    freebayes: Dict[str, Path],
    strelka: Dict[str, Path],
    joint_hc_vcf: Optional[Path],
    absolute_paths: bool,
) -> None:
    rows: List[Tuple[str, str, str, str]] = []

    discovered_small_samples: Set[str] = set()
    discovered_small_samples.update(deepvariant.keys())
    discovered_small_samples.update(freebayes.keys())
    discovered_small_samples.update(strelka.keys())

    for sample in sorted(discovered_small_samples):
        if sample in deepvariant:
            rows.append(
                (
                    make_unique_id(sample, "deepvariant"),
                    rel_or_abs(deepvariant[sample], absolute_paths),
                    "deepvariant",
                    "",
                )
            )
        if sample in freebayes:
            rows.append(
                (
                    make_unique_id(sample, "freebayes"),
                    rel_or_abs(freebayes[sample], absolute_paths),
                    "freebayes",
                    "",
                )
            )
        if sample in strelka:
            rows.append(
                (
                    make_unique_id(sample, "strelka"),
                    rel_or_abs(strelka[sample], absolute_paths),
                    "strelka",
                    "",
                )
            )

    if joint_hc_vcf is not None:
        joint_samples = get_vcf_samples_with_bcftools(joint_hc_vcf)

        if joint_samples is None:
            print(
                "WARN: bcftools unavailable or joint VCF sample query failed; "
                "falling back to sample_sample for HaplotypeCaller subsample values.",
                file=sys.stderr,
            )
            for sample in sorted(discovered_small_samples):
                rows.append(
                    (
                        make_unique_id(sample, "haplotypecaller"),
                        rel_or_abs(joint_hc_vcf, absolute_paths),
                        "haplotypecaller",
                        f"{sample}_{sample}",
                    )
                )
        else:
            joint_sample_set = set(joint_samples)

            for sample in sorted(discovered_small_samples):
                subsample = resolve_haplotypecaller_subsample(sample, joint_sample_set)
                if subsample is None:
                    print(
                        f"WARN: could not resolve HaplotypeCaller subsample for {sample}. "
                        f"Tried '{sample}' and '{sample}_{sample}'. Skipping row.",
                        file=sys.stderr,
                    )
                    continue

                rows.append(
                    (
                        make_unique_id(sample, "haplotypecaller"),
                        rel_or_abs(joint_hc_vcf, absolute_paths),
                        "haplotypecaller",
                        subsample,
                    )
                )

    rows.sort(
        key=lambda x: (
            x[0].rsplit("_", 1)[0] if "_" in x[0] else x[0],
            CALLER_ORDER_SMALL.get(x[2], 999),
            x[1],
        )
    )

    with out_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", "test_vcf", "caller", "subsample"])
        writer.writerows(rows)


def write_structural_samplesheet(
    out_csv: Path,
    manta: Dict[str, Path],
    tiddit: Dict[str, Path],
    absolute_paths: bool,
) -> None:
    rows: List[Tuple[str, str, str]] = []

    discovered_struct_samples: Set[str] = set()
    discovered_struct_samples.update(manta.keys())
    discovered_struct_samples.update(tiddit.keys())

    for sample in sorted(discovered_struct_samples):
        if sample in manta:
            rows.append(
                (
                    make_unique_id(sample, "manta"),
                    rel_or_abs(manta[sample], absolute_paths),
                    "manta",
                )
            )
        if sample in tiddit:
            rows.append(
                (
                    make_unique_id(sample, "tiddit"),
                    rel_or_abs(tiddit[sample], absolute_paths),
                    "tiddit",
                )
            )

    rows.sort(
        key=lambda x: (
            x[0].rsplit("_", 1)[0] if "_" in x[0] else x[0],
            CALLER_ORDER_STRUCTURAL.get(x[2], 999),
            x[1],
        )
    )

    with out_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", "test_vcf", "caller"])
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate nf-core/variantbenchmarking samplesheets from nf-core/sarek results."
    )
    parser.add_argument(
        "--project-dir",
        default=".",
        help="Path to sarek-variant-calling project directory (default: current directory)",
    )
    parser.add_argument(
        "--small-out",
        default="benchmark_small_samplesheet.csv",
        help="Output CSV for small variants",
    )
    parser.add_argument(
        "--structural-out",
        default="benchmark_structural_samplesheet.csv",
        help="Output CSV for structural variants",
    )
    parser.add_argument(
        "--absolute-paths",
        action="store_true",
        help="Write absolute paths instead of project-relative paths",
    )
    args = parser.parse_args()

    global PROJECT_DIR
    PROJECT_DIR = Path(args.project_dir).resolve()

    results_dir = PROJECT_DIR / "results" / "variant_calling"
    if not results_dir.exists():
        raise SystemExit(f"ERROR: expected directory not found: {results_dir}")

    deepvariant = find_first_per_sample(
        "results/variant_calling/deepvariant/*/*.deepvariant.vcf.gz"
    )
    freebayes = find_first_per_sample(
        "results/variant_calling/freebayes/*/*.freebayes.filtered.vcf.gz"
    )
    strelka = find_first_per_sample(
        "results/variant_calling/strelka/*/*.strelka.variants.vcf.gz"
    )
    manta = find_first_per_sample(
        "results/variant_calling/manta/*/*.manta.diploid_sv.vcf.gz"
    )
    tiddit = find_first_per_sample(
        "results/variant_calling/tiddit/*/*.tiddit.vcf.gz"
    )
    joint_hc_vcf = find_joint_haplotypecaller_vcf()

    small_out = PROJECT_DIR / args.small_out
    structural_out = PROJECT_DIR / args.structural_out

    write_small_samplesheet(
        out_csv=small_out,
        deepvariant=deepvariant,
        freebayes=freebayes,
        strelka=strelka,
        joint_hc_vcf=joint_hc_vcf,
        absolute_paths=args.absolute_paths,
    )

    write_structural_samplesheet(
        out_csv=structural_out,
        manta=manta,
        tiddit=tiddit,
        absolute_paths=args.absolute_paths,
    )

    print(f"Wrote: {small_out}")
    print(f"Wrote: {structural_out}")
    print()
    print("Discovered small-variant inputs:")
    print(f"  deepvariant      : {len(deepvariant)}")
    print(f"  freebayes        : {len(freebayes)}")
    print(f"  strelka          : {len(strelka)}")
    print(f"  haplotypecaller  : {'1 joint VCF' if joint_hc_vcf else 'not found'}")
    print()
    print("Discovered structural-variant inputs:")
    print(f"  manta            : {len(manta)}")
    print(f"  tiddit           : {len(tiddit)}")


if __name__ == "__main__":
    PROJECT_DIR = Path(".").resolve()
    main()
