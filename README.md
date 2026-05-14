# sarek-variant-calling on RTU HPC (PBS)

Pipeline for germline variant calling, benchmarking, annotation, and post-processing on the RTU HPC cluster using [nf-core/sarek](https://nf-co.re/sarek/3.8.1/), [Viktorija0719/germline-nextflow](https://github.com/Viktorija0719/germline-nextflow), and [nf-core/variantbenchmarking](https://nf-co.re/variantbenchmarking/1.4.0/).

Designed for **RTU HPC with PBS**, uses **shared RSU reference resources** under `/home/groups/rsu/dauksaite_v/`, and avoids hardcoding personal paths where possible.

---

## Overview

The full pipeline (`run_all.pbs`) runs 8 sequential steps:

| Step | Description |
|------|-------------|
| 1 | `nf-core/sarek` — FASTQ → BAM, CRAM, initial VCF calling |
| 2 | `generate_bam_samplesheet.py` — build `samplesheet_bam.csv` from Sarek BAM outputs |
| 3 | `Viktorija0719/germline-nextflow` — DeepVariant + Strelka2 small variant calling matching `pipeline_nextflow.txt` |
| 4 | `generate_variantbenchmarking_samplesheets.py` — prepare benchmark inputs |
| 5A | `nf-core/variantbenchmarking` — small variant benchmarking |
| 5B | `nf-core/variantbenchmarking` — structural variant benchmarking |
| 6 | `generate_samplesheet_annotate_tp.py` — prepare TP annotation inputs |
| 7 | `nf-core/sarek` — annotate TP VCFs (VEP / snpEff) |
| 8 | `merge_snv_vcfs.py` + `merge_sv_vcfs.py` — merge annotated callsets to Excel/TSV |

The workflow is intended for **whole-exome sequencing (WES)**.

---

## Requirements

### RTU HPC modules

- `java/jdk-21.0.2`
- `python/3.9.19`
- `singularity`

### Shared reference resources

Expected under:

```
/home/groups/rsu/dauksaite_v/
```

Includes iGenomes reference files, truth VCFs, benchmark BED files, annotation resources, and VEP/snpEff caches.

---

## Installation

### 1. Install Nextflow on the login node

```bash
module load java/jdk-21.0.2

wget -qO- https://get.nextflow.io | bash
mkdir -p ~/bin
mv nextflow ~/bin
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

nextflow info
```

### 2. Pull the workflows

```bash
nextflow pull nf-core/sarek -r 3.8.1
nextflow pull nf-core/variantbenchmarking -r 1.4.0
nextflow pull Viktorija0719/germline-nextflow
```

### 3. Create the Python environment

```bash
module purge
module load python/3.9.19

python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
pip install -r requirements.txt
```

---

## Input files

### `samplesheet.csv` (Step 1 — Sarek FASTQ input)

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
ERR031932,XY,0,ERR031932,L001,/path/to/ERR031932_1.fastq.gz,/path/to/ERR031932_2.fastq.gz
ERR031933,XY,0,ERR031933,L001,/path/to/ERR031933_1.fastq.gz,/path/to/ERR031933_2.fastq.gz
```

See the [nf-core/sarek samplesheet documentation](https://nf-co.re/sarek/3.8.1/docs/usage/#input-sample-sheet-configurations).

### `params/params_bam.yaml` (Step 3 — germline-nextflow)

Must point `input` to the BAM samplesheet generated in Step 2:

```yaml
input: samplesheet_bam.csv
variant_target_bed: "/path/to/your/target_regions.bed"
genome: GATK.GRCh38
outdir: results_germline
# ... other germline-nextflow params
```

### `conf/match_pipeline_nextflow.config` (Step 3)

Applies DeepVariant + Strelka2 settings matching `pipeline_nextflow.txt`:

- `no_intervals = true` — one DeepVariant job per sample on full BED
- `wes = true` — WES model
- `deepvariant_args = '--postprocess_variants_extra_args=only_keep_pass=true'`
- `combine_dv_strelka_enable = true` — merges DV + Strelka2 VCFs
- `cpus = 64` for DeepVariant (`--num_shards=64`)
- `ext.args = '--exome --callContinuousVf chrM'` for Strelka2

Reduce `cpus` to `24` or `16` if 64-core nodes are unavailable on the cluster.

---

## Running

Submit the full pipeline:

```bash
qsub run_all.pbs
```

The test samples used in this project come from:

- **Project:** `PRJNA60113`
- **Description:** Exome sequencing of GBR (British from England and Scotland) HapMap population

---

## Outputs

### Intermediate

| File | Produced by |
|------|-------------|
| `samplesheet_bam.csv` | Step 2 — BAM paths from Sarek cram_to_bam output |
| `benchmark_small_samplesheet.csv` | Step 4 |
| `benchmark_structural_samplesheet.csv` | Step 4 |
| `samplesheet_annotate_tp.csv` | Step 6 |

### Final merged callsets (Step 8)

- `merged_snv.xlsx` + `merged_snv/*.tsv` — one sheet/TSV per sample
- `merged_sv.xlsx` + `merged_sv/*.tsv` — one sheet/TSV per sample

---

## Logs

All logs land in `logs/`:

| Log file | Step |
|----------|------|
| `report.html`, `timeline.html`, `trace.txt` | Step 1 (Sarek) |
| `germline_report.html`, `germline_timeline.html`, `germline_trace.txt` | Step 3 (germline-nextflow) |
| `variantbenchmarking_small_*.html/txt` | Step 5A |
| `variantbenchmarking_structural_*.html/txt` | Step 5B |
| `annotate_tp_*.html/txt` | Step 7 |

PBS job output: `*.o<jobid>` in the working directory.

---

## Troubleshooting

### Check job status

```bash
qstat
```

### Resume failed runs

All steps use `-resume`, so reruns reuse completed Nextflow work.

### Paths

- User caches: `/home_beegfs/$USER/.nextflow`, `/home_beegfs/$USER/.singularity`
- Work directories: `/home_beegfs/$USER/nxf_work/<project>*`
- Shared RSU resources: `/home/groups/rsu/dauksaite_v/`

---

## Acknowledgements

- [nf-core/sarek](https://nf-co.re/sarek/3.8.1/)
- [Viktorija0719/germline-nextflow](https://github.com/Viktorija0719/germline-nextflow)
- [nf-core/variantbenchmarking](https://nf-co.re/variantbenchmarking/1.4.0/)
