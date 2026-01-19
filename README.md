# JuxtaBAM

JuxtaBAM is a species and variant set-independent Python tool for cross-assay genotype concordance and sample identity validation using joint SNV genotyping. It supports various epigenomics and transcriptomic assays (i.e., ATAC-seq, DNase-seq, RNA-seq, etc.) and can work with different species and reference genomes (hg38, mm10, etc.).

## Requirements
- Python >= 3.8
- GATK >= v4.1.9.0
- bcftools
- samtools
- Slurm (optional)

## Running genotype comparison
python scripts/juxtabam_run.py [arguments]

## Generating plots
python scripts/juxtabam_plot.py [arguments]

## Notes
- The main pipeline script manages Slurm submission internally when enabled.
- Plotting is a separate step by design.
