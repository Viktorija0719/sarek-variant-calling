#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def rel_from_project(path: Path) -> str:
    """
    Always write project-relative paths, e.g.
    results/preprocessing/converted/cram_to_bam/S5139Nr29/S5139Nr29.bam
    """
    return str(path.relative_to(PROJECT_DIR))


def find_bam_pairs(pattern: str) -> Dict[str, Tuple[Path, Path]]:
    """
    Find BAM/BAI pairs and map them by sample directory name.

    Expected layout:
      results/preprocessing/converted/cram_to_bam/S5139Nr29/S5139Nr29.bam
      results/preprocessing/converted/cram_to_bam/S5139Nr29/S5139Nr29.bam.bai
    """
    matches = sorted(PROJECT_DIR.glob(pattern))
    result: Dict[str, Tuple[Path, Path]] = {}

    for bam_path in matches:
        if bam_path.name.endswith(".bam.bai"):
            continue

        sample = bam_path.parent.name
        bai_path = bam_path.with_name(bam_path.name + ".bai")

        if not bai_path.exists():
            print(
                f"WARN: missing BAI for sample={sample}, BAM={bam_path}. Skipping.",
                file=sys.stderr,
            )
            continue

        if sample not in result:
            result[sample] = (bam_path, bai_path)
        else:
            print(
                f"WARN: multiple BAMs found for sample={sample}. "
                f"Keeping {result[sample][0]} and ignoring {bam_path}",
                file=sys.stderr,
            )

    return result


def write_samplesheet(
    out_csv: Path,
    bam_pairs: Dict[str, Tuple[Path, Path]],
) -> None:
    rows: List[Tuple[str, str, str, str]] = []

    for sample in sorted(bam_pairs):
        bam_path, bai_path = bam_pairs[sample]
        rows.append(
            (
                sample,                  # patient
                sample,                  # sample
                rel_from_project(bam_path),
                rel_from_project(bai_path),
            )
        )

    with out_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["patient", "sample", "bam", "bai"])
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate BAM samplesheet from Sarek results/preprocessing/converted/cram_to_bam."
    )
    parser.add_argument(
        "--project-dir",
        default=".",
        help="Path to sarek project directory (default: current directory)",
    )
    parser.add_argument(
        "--out",
        default="samplesheet_2.csv",
        help="Output CSV filename",
    )
    args = parser.parse_args()

    global PROJECT_DIR
    PROJECT_DIR = Path(args.project_dir).resolve()

    converted_dir = PROJECT_DIR / "results" / "preprocessing" / "converted" / "cram_to_bam"
    if not converted_dir.exists():
        raise SystemExit(f"ERROR: expected directory not found: {converted_dir}")

    bam_pairs = find_bam_pairs("results/preprocessing/converted/cram_to_bam/*/*.bam")

    if not bam_pairs:
        raise SystemExit(
            "ERROR: no BAM/BAI pairs found under results/preprocessing/converted/cram_to_bam"
        )

    out_csv = PROJECT_DIR / args.out
    write_samplesheet(out_csv=out_csv, bam_pairs=bam_pairs)

    print(f"Wrote: {out_csv}")
    print(f"Discovered BAM samples: {len(bam_pairs)}")
    print("Paths written as project-relative paths starting with results/")


if __name__ == "__main__":
    PROJECT_DIR = Path(".").resolve()
    main()
