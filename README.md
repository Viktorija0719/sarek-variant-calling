# nf-core/sarek on RTU HPC (PBS)

[nf-core/sarek](https://nf-co.re/sarek/3.8.1/) is a workflow designed to detect variants on whole genome or targeted sequencing data. Sarek can also handle tumour / normal pairs and could include additional relapses.

This setup keeps shared RSU resources in `/home/groups/rsu/dauksaite_v/` and avoids hardcoding personal `/home_beegfs/<user>/…` paths (uses project-relative paths + `$USER` in the PBS script).

## 1) Install Nextflow (login node)

Requires: `java/jdk-21.0.2`.

```bash
module load java/jdk-21.0.2

wget -qO- https://get.nextflow.io | bash
mkdir -p ~/bin && mv nextflow ~/bin
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

nextflow info
````

## 2) Pull the pipeline

```bash
nextflow pull nf-core/sarek
nextflow pull nf-core/variantbenchmarking

module purge
module load python/3.9.19

python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
pip install -r requirements.txt

```

## Quick test run (recommended)

To test that everything works on RTU HPC, submit the PBS script:

```bash
qsub /path/to/run_all.pbs
```

The test samples used in this project come from:

**Project:** PRJNA60113

Exome sequencing of (GBR) British from England and Scotland HapMap population



## Run germline variant calling

1. Create `samplesheet.csv` in your working directory following:
   [https://nf-co.re/sarek/3.8.1/docs/usage/#input-sample-sheet-configurations](https://nf-co.re/sarek/3.8.1/docs/usage/#input-sample-sheet-configurations)

2. Submit the run:

```bash
qsub run_all.pbs
```

## Expected outputs

`results/` will contain:

* `preprocessing/` (alignment, markduplicates, BQSR, CRAM/BAM)
* `variant_calling/` (per-tool callsets + `haplotypecaller/joint_variant_calling/` if enabled)
* `annotation/` (VEP / snpEff annotated results per caller)
* `multiqc/multiqc_report.html`
* `pipeline_info/` (params JSON, DAGs, software versions)

Tool folders appear only for enabled tools (e.g., deepvariant/freebayes/haplotypecaller/manta/tiddit/cnvkit/indexcov/vep/snpeff).




**Outputs**

* `merged_sv.xlsx` + `merged_sv/*.tsv` (one sheet / TSV per sample)
* `merged_snv.xlsx` + `merged_snv/*.tsv` (one sheet / TSV per sample)






















# sarek-variant-calling on RTU HPC (PBS)

Pipeline for germline variant calling, benchmarking, annotation, and post-processing on the RTU HPC cluster using [nf-core/sarek](https://nf-co.re/sarek/3.8.1/) and [nf-core/variantbenchmarking](https://nf-co.re/variantbenchmarking/1.4.0/).

This setup is designed for **RTU HPC with PBS**, uses **shared RSU reference resources** under `/home/groups/rsu/dauksaite_v/`, and avoids hardcoding personal paths where possible.

## Overview

This project supports:

- germline variant calling from FASTQ input with `nf-core/sarek`
- benchmarking of small and structural variants with `nf-core/variantbenchmarking`
- annotation of selected VCFs
- merging annotated SNV/indel and SV callsets into per-sample tables

The workflow is intended for **whole-exome sequencing (WES)** but can be adapted for other targeted data.



## Requirements

### RTU HPC environment

* Java module:

  * `java/jdk-21.0.2`
* Python module:

  * `python/3.9.19`

### Shared reference resources

This setup expects shared RSU resources under:

```bash
/home/groups/rsu/dauksaite_v/
```

These include:

* iGenomes reference files
* truth VCFs and benchmark BED files
* annotation resources
* cache files for tools such as VEP and snpEff

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
```

### 3. Create the Python environment for downstream merging scripts

```bash
module purge
module load python/3.9.19

python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
pip install -r requirements.txt
```

## Quick test run

To test that the setup works on RTU HPC, submit:

```bash
qsub run_all.pbs
```

The test samples used in this project come from:

* **Project:** `PRJNA60113`
* **Description:** Exome sequencing of GBR (British from England and Scotland) HapMap population

## Input

### Main Sarek samplesheet

Create `samplesheet.csv` in the working directory.

Example:

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
ERR031932,XY,0,ERR031932,L001,/path/to/ERR031932_1.fastq.gz,/path/to/ERR031932_2.fastq.gz
ERR031933,XY,0,ERR031933,L001,/path/to/ERR031933_1.fastq.gz,/path/to/ERR031933_2.fastq.gz
```

See the official Sarek documentation for supported sample sheet formats:

* [nf-core/sarek input samplesheet documentation](https://nf-co.re/sarek/3.8.1/docs/usage/#input-sample-sheet-configurations)





## Final merged outputs

Post-processing produces:

* `merged_sv.xlsx`
* `merged_sv/*.tsv`
* `merged_snv.xlsx`
* `merged_snv/*.tsv`

These are organized as:

* **one Excel workbook per variant type**
* **one sheet per sample**
* **one TSV per sample**

## Notes

* Shared RSU resources are kept under `/home/groups/rsu/dauksaite_v/`.
* User-specific caches and work directories are created under `/home_beegfs/$USER/`.
* Project-level scripts use project-relative paths where possible.
* Large temporary work directories must be visible from compute nodes.

## Example user-specific cache and work paths

Used in PBS scripts:

```bash
NXF_HOME=/home_beegfs/$USER/.nextflow
NXF_SINGULARITY_CACHEDIR=/home_beegfs/$USER/.singularity
/home_beegfs/$USER/nxf_work/
```

## Troubleshooting

### Check job status

```bash
qstat
```

### Inspect logs

PBS output and workflow logs are useful for debugging:

* `*.o<jobid>`
* `logs/report.html`
* `logs/timeline.html`
* `logs/trace.txt`
* `.nextflow.log`

### Resume failed runs

All provided scripts use:

```bash
-resume
```

so reruns will reuse completed steps when possible.

## Acknowledgements

This project uses:

* [nf-core/sarek](https://nf-co.re/sarek/3.8.1/)
* [nf-core/variantbenchmarking](https://nf-co.re/variantbenchmarking/1.4.0/)







