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
nextflow pull Viktorija0719/sarek -r varlofix-3.8.1
```

## Quick test run (recommended)

To test that everything works on RTU HPC, submit the PBS script:

```bash
qsub /path/to/run_sarek.pbs
```

The test samples used in this project come from:

**Project:** PRJNA60113



## Run Sarek

1. Create `samplesheet.csv` in your working directory following:
   [https://nf-co.re/sarek/3.8.1/docs/usage/#input-sample-sheet-configurations](https://nf-co.re/sarek/3.8.1/docs/usage/#input-sample-sheet-configurations)

2. Submit the run:

```bash
qsub run_sarek.pbs
```

## Expected outputs

`results/` will contain:

* `preprocessing/` (alignment, markduplicates, BQSR, CRAM/BAM)
* `variant_calling/` (per-tool callsets + `haplotypecaller/joint_variant_calling/` if enabled)
* `annotation/` (VEP / snpEff annotated results per caller)
* `multiqc/multiqc_report.html`
* `pipeline_info/` (params JSON, DAGs, software versions)

Tool folders appear only for enabled tools (e.g., deepvariant/freebayes/haplotypecaller/manta/tiddit/cnvkit/indexcov/vep/snpeff).


## Combine annotated VCFs into convenient per-sample tables

After the pipeline finishes, the annotated VCFs are in `results/annotation/` (grouped by caller and sample). The scripts below merge them into **one table per sample** (Excel with one sheet per sample + per-sample TSVs):

* **`merged_sv`**: structural variants from **Manta + TIDDIT**
* **`merged_snv`**: small variants from **DeepVariant + FreeBayes + Strelka**

```bash
# From the project root 
module purge
module load python/3.9.19
source .venv/bin/activate

# 1) Merge SV callsets (Manta + TIDDIT)
python3 scripts/merge_sv_vcfs.py \
  --annotation-dir results/annotation \
  --callers manta,tiddit \
  --emit-split \
  --out-xlsx merged_sv.xlsx \
  --out-tsv-dir merged_sv

# 2) Merge SNV/indel callsets (DeepVariant + FreeBayes + Strelka)
python3 scripts/merge_snv_vcfs.py \
  --annotation-dir results/annotation \
  --callers deepvariant,freebayes,strelka \
  --out-xlsx merged_snv.xlsx \
  --out-tsv-dir merged_snv
```

**Outputs**

* `merged_sv.xlsx` + `merged_sv/*.tsv` (one sheet / TSV per sample)
* `merged_snv.xlsx` + `merged_snv/*.tsv` (one sheet / TSV per sample)
