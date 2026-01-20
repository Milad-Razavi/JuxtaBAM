# JuxtaBAM

JuxtaBAM (juxtaposing BAM files) is a Python-based tool for cross-assay genotype concordance and sample identity validation using joint SNV genotyping. It is designed to work across different sequencing assays such as RNA-seq, ATAC-seq, and DNase-seq, and is agnostic to species, reference genome, and variant sets.

The primary goal of JuxtaBAM is to verify sample identity across assays, detect swaps or mislabeling, and quantify genotype mismatch ratesk.

---

## Key features

* Joint SNV genotyping using GATK
* Cross-assay and within-assay genotype concordance
* Species and reference independent (hg38, mm10, others)
* Optional Slurm support for parallel execution
* Downstream visualization tools

---

## Requirements

* Python >= 3.9
* pandas
* GATK >= 4.1.9.0
* samtools
* bcftools
* Slurm (optional, for distributed execution)

All external tools must be available in your PATH, unless explicit paths are provided via command-line arguments.

---

## Input specification

JuxtaBAM requires a tab-separated (or comma-separated) input manifest describing the samples to be compared.

Required columns:

* sample_id: unique identifier for each sample
* biosample: biological source identifier (shared across assays for the same specimen)
* assay: assay name (for example RNA-seq, ATAC-seq)
* bam: path to a coordinate-sorted BAM file

Example input_info.tsv:

sample_id&emsp;biosample&emsp;assay&emsp;bam  
DNase_CH12LX_ENCFF668NSZ&emsp;CH12LX&emsp;DNase&emsp;/path/to/DNase_CH12LX_ENCFF668NSZ_chr1.bam
DNase_A20_ENCFF129FUZ&emsp;A20&emsp;DNase&emsp;/path/to/DNase_A20_ENCFF129FUZ_chr1.bam
RNA_CH12LX_ENCFF203XTH&emsp;CH12LX&emsp;RNA&emsp;/path/to/RNA_CH12LX_ENCFF203XTH_chr1.bam
RNA_A20_SRR19241698&emsp;A20&emsp;RNA&emsp;/path/to/RNA_A20_SRR19241698_chr1.bam


Each BAM file must:

* Be sorted
* Have a corresponding .bai index
* Pass samtools quickcheck
* Contain valid or fixable read group information

The reference genome file (.fa):
* Must have a .dict file in the same directory
* Must have a .fai index file in the same directory

Read more about the reference genome dict and index file: https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
---

## Running the pipeline

The main pipeline entry point is:

python scripts/juxtabam_run.py [arguments]

### Minimal example

python scripts/juxtabam_run.py 
--input_info input_info.tsv 
--reference mm10.fa 
--output_dir juxtabam_output

### Common options

* --input_info: input manifest TSV file (required)
* --reference: reference FASTA file (required). Must have .fai and .dict
* --output_dir: base output directory (required)
* --regions: genomic regions to analyze. Can be chr1 or chr15,chr16,chr17 . Or comma-separated list, BED file, or VCF file
* --threads: number of threads for GATK operations
* --java_mem: Java memory setting passed to GATK (default="-Xmx32G")


### Other options
#### Read depth and variant filters

* --min-dp: minimum read depth per SNP per sample (default: 20)
* --min-qd: minimum QD value for variant filtering (default: 2.0)
* --min-qual: minimum QUAL value for variant filtering (default: 30)

#### BAM preprocessing controls

* --skip_AddOrReplaceRG: skip AddOrReplaceReadGroups
* --skip_SplitNCigar: skip SplitNCigarReads
* --call_SplitNCigar_for_assay_types: comma-separated assay names for which SplitNCigarReads is applied

#### Tool paths

* --gatk_path: path to GATK executable (only if not in your system path)
* --samtools_path: path to samtools executable
* --bcftools_path: path to bcftools executable

---

## Slurm support

JuxtaBAM can optionally offload expensive steps to Slurm.

### Slurm-enabled options

* --slurm_haplotypecaller: run per-sample HaplotypeCaller jobs via sbatch
* --slurm_threads: CPUs per Slurm HaplotypeCaller job
* --slurm_ram: memory per Slurm HaplotypeCaller job
* --slurm_pairwise_bcftools: run pairwise bcftools queries on Slurm

When Slurm options are enabled, JuxtaBAM automatically generates sbatch scripts, submits jobs, and waits for completion. No manual Slurm scripting is required.

### Example Slurm run

python scripts/juxtabam_run.py 
--input_info data/example_mouse_ENCODE/input_info.tsv 
--reference data/mm10_data/mm10.fa 
--output_dir data/example_mouse_ENCODE/outputs 
--regions chr1 
--threads 6 
--min-dp 10 
--slurm_haplotypecaller 
--slurm_threads 1 
--slurm_ram 40G 
--slurm_pairwise_bcftools

The same command runs locally if Slurm flags are omitted.

---

## Output files

The main pipeline produces:

* gvcf/: per-sample GVCF files
* gendb/: GenomicsDB workspace
* joint/: joint VCFs and filtered SNP sets
* genotype_mismatch_rates.tsv: pairwise mismatch summary
* genotype_mismatch_rates.csv: CSV version of results
* nearest_neighbor_matches.tsv: closest genotype matches per sample
* run_info.txt: run provenance and parameters

Mismatch metrics include:

* total_sites: number of usable SNPs
* mismatches: number of discordant genotypes
* percent_mismatch: mismatch rate

Warnings are emitted if comparisons are based on few SNPs.

---

## Plotting and visualization

Plot generation is handled by a separate script:

python scripts/juxtabam_plot.py [arguments]

Required inputs:

* --genotype_mismatch_rates: output TSV from juxtabam_run.py
* --input_info: same input manifest used for the run

Optional parameters:

* --outdir: output directory for plots
* --min_annot_percent: annotation threshold for heatmaps
* --no_line_annotations: disable line plot annotations

Generated figures include:

* Cross-assay heatmaps
* Bubble heatmaps
* All-vs-all mismatch matrices
* SNV site count distributions
* Nearest-neighbor line plots

---

## Notes and recommendations

* Strongly recommended: Limit regions during testing to reduce runtime (e.g., --regions chr19 or --regions chr1,chr2,chr3) and excluding the sex chromosomes
* RNA-seq and other spliced assays should use SplitNCigarReads
* Interpret mismatch rates in the context of coverage and other assays
* Very low SNP counts may lead to unstable estimates
* If the pipeline fails halfway, make sure you remove the partially-generated corrupted files (or the entire gendb/ directory for a failure in the GenomicImportDB step), before running again

---

## License

MIT License. See the LICENSE file for details.

---

## Authors

Milad Razavi-Mohseni, Michael A Beer
