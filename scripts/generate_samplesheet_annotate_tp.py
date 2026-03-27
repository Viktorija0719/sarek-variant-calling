#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

CALLER_ORDER = {
    "deepvariant": 0,
    "freebayes": 1,
    "haplotypecaller": 2,
    "strelka": 3,
    "manta": 4,
    "tiddit": 5,
}


def make_sample_name(patient: str, caller: str) -> str:
    return f"{patient}_{caller}"


def rel_or_abs(path: Path, absolute: bool) -> str:
    return str(path.resolve()) if absolute else str(path.relative_to(PROJECT_DIR))


def find_one_per_patient(pattern: str) -> Dict[str, Path]:
    """
    Find files matching a pattern and map them by patient/sample directory name.
    For paths like:
      results_variantbenchmarking/small/small/ERR031932_deepvariant/benchmarks/rtgtools/...tp.vcf.gz
    patient becomes ERR031932 (taken from directory stem before final _caller suffix).
    """
    matches = sorted(PROJECT_DIR.glob(pattern))
    result: Dict[str, Path] = {}

    for path in matches:
        if path.name.endswith(".tbi"):
            continue

        parent_dir = path.parents[2].name  # e.g. ERR031932_deepvariant
        if "_" not in parent_dir:
            print(
                f"WARN: could not parse patient from directory name: {parent_dir}",
                file=sys.stderr,
            )
            continue

        patient = parent_dir.rsplit("_", 1)[0]

        if patient not in result:
            result[patient] = path
        else:
            print(
                f"WARN: multiple matches for patient={patient} pattern={pattern}. "
                f"Keeping {result[patient]} and ignoring {path}",
                file=sys.stderr,
            )
    return result


def write_samplesheet(
    out_csv: Path,
    deepvariant: Dict[str, Path],
    freebayes: Dict[str, Path],
    haplotypecaller: Dict[str, Path],
    strelka: Dict[str, Path],
    manta: Dict[str, Path],
    tiddit: Dict[str, Path],
    absolute_paths: bool,
) -> None:
    rows: List[Tuple[str, str, str, str]] = []

    patients: Set[str] = set()
    patients.update(deepvariant.keys())
    patients.update(freebayes.keys())
    patients.update(haplotypecaller.keys())
    patients.update(strelka.keys())
    patients.update(manta.keys())
    patients.update(tiddit.keys())

    for patient in sorted(patients):
        if patient in deepvariant:
            rows.append(
                (
                    patient,
                    make_sample_name(patient, "deepvariant"),
                    "deepvariant",
                    rel_or_abs(deepvariant[patient], absolute_paths),
                )
            )
        if patient in freebayes:
            rows.append(
                (
                    patient,
                    make_sample_name(patient, "freebayes"),
                    "freebayes",
                    rel_or_abs(freebayes[patient], absolute_paths),
                )
            )
        if patient in haplotypecaller:
            rows.append(
                (
                    patient,
                    make_sample_name(patient, "haplotypecaller"),
                    "haplotypecaller",
                    rel_or_abs(haplotypecaller[patient], absolute_paths),
                )
            )
        if patient in strelka:
            rows.append(
                (
                    patient,
                    make_sample_name(patient, "strelka"),
                    "strelka",
                    rel_or_abs(strelka[patient], absolute_paths),
                )
            )
        if patient in manta:
            rows.append(
                (
                    patient,
                    make_sample_name(patient, "manta"),
                    "manta",
                    rel_or_abs(manta[patient], absolute_paths),
                )
            )
        if patient in tiddit:
            rows.append(
                (
                    patient,
                    make_sample_name(patient, "tiddit"),
                    "tiddit",
                    rel_or_abs(tiddit[patient], absolute_paths),
                )
            )

    rows.sort(key=lambda x: (x[0], CALLER_ORDER.get(x[2], 999), x[3]))

    with out_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["patient", "sample", "variantcaller", "vcf"])
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate samplesheet_annotate_tp.csv from results_variantbenchmarking outputs."
    )
    parser.add_argument(
        "--project-dir",
        default=".",
        help="Path to project directory (default: current directory)",
    )
    parser.add_argument(
        "--out",
        default="samplesheet_annotate_tp.csv",
        help="Output CSV filename",
    )
    parser.add_argument(
        "--absolute-paths",
        action="store_true",
        help="Write absolute paths instead of project-relative paths",
    )
    args = parser.parse_args()

    global PROJECT_DIR
    PROJECT_DIR = Path(args.project_dir).resolve()

    results_dir = PROJECT_DIR / "results_variantbenchmarking"
    if not results_dir.exists():
        raise SystemExit(f"ERROR: expected directory not found: {results_dir}")

    # small: tp.vcf.gz
    deepvariant = find_one_per_patient(
        "results_variantbenchmarking/small/small/*_deepvariant/benchmarks/rtgtools/*.deepvariant.tp.vcf.gz"
    )
    freebayes = find_one_per_patient(
        "results_variantbenchmarking/small/small/*_freebayes/benchmarks/rtgtools/*.freebayes.tp.vcf.gz"
    )
    haplotypecaller = find_one_per_patient(
        "results_variantbenchmarking/small/small/*_haplotypecaller/benchmarks/rtgtools/*.haplotypecaller.tp.vcf.gz"
    )
    strelka = find_one_per_patient(
        "results_variantbenchmarking/small/small/*_strelka/benchmarks/rtgtools/*.strelka.tp.vcf.gz"
    )

    # structural: tp-comp.vcf.gz
    manta = find_one_per_patient(
        "results_variantbenchmarking/structural/structural/*_manta/benchmarks/truvari/*.manta.tp-comp.vcf.gz"
    )
    tiddit = find_one_per_patient(
        "results_variantbenchmarking/structural/structural/*_tiddit/benchmarks/truvari/*.tiddit.tp-comp.vcf.gz"
    )

    out_csv = PROJECT_DIR / args.out

    write_samplesheet(
        out_csv=out_csv,
        deepvariant=deepvariant,
        freebayes=freebayes,
        haplotypecaller=haplotypecaller,
        strelka=strelka,
        manta=manta,
        tiddit=tiddit,
        absolute_paths=args.absolute_paths,
    )

    print(f"Wrote: {out_csv}")
    print()
    print("Discovered small TP VCFs:")
    print(f"  deepvariant      : {len(deepvariant)}")
    print(f"  freebayes        : {len(freebayes)}")
    print(f"  haplotypecaller  : {len(haplotypecaller)}")
    print(f"  strelka          : {len(strelka)}")
    print()
    print("Discovered structural TP-comp VCFs:")
    print(f"  manta            : {len(manta)}")
    print(f"  tiddit           : {len(tiddit)}")


if __name__ == "__main__":
    PROJECT_DIR = Path(".").resolve()
    main()
