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
```

## Run Sarek

1. Create `samplesheet.csv` in your working directory following:
   [https://nf-co.re/sarek/3.8.1/docs/usage/#input-sample-sheet-configurations](https://nf-co.re/sarek/3.8.1/docs/usage/#input-sample-sheet-configurations)

2. Submit the run:

```bash
qsub run_sarek.pbs
```

## Expected outputs (high level)

`results/` will contain:

* `preprocessing/` (alignment, markduplicates, BQSR, CRAM/BAM)
* `variant_calling/` (per-tool callsets + `haplotypecaller/joint_variant_calling/` if enabled)
* `annotation/` (VEP / snpEff annotated results per caller)
* `multiqc/multiqc_report.html`
* `pipeline_info/` (params JSON, DAGs, software versions)

Tool folders appear only for enabled tools (e.g., deepvariant/freebayes/haplotypecaller/manta/tiddit/cnvkit/indexcov/vep/snpeff).


