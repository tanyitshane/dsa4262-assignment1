# DSA4262 Assignment 1 – Long-read RNA-Seq pipeline

## Pipeline (Task 5)
- `nextflow/main.nf` – ALIGN (minimap2) → SAM2BAM (samtools) → QC → BAMBU
- `nextflow/bin/run_bambu.R` – Bambu step
- `environment.yml` – conda environment with tool versions

Run from a directory containing `fastq/` and `reference/`:

    # Scenario 1: with annotations
    nextflow run nextflow/main.nf -with-report report_scenario1.html --outdir results_scenario1 --useAnnotations true
    # Scenario 2: without annotations, reusing cached steps
    nextflow run nextflow/main.nf -resume -with-report report_scenario2.html --outdir results_scenario2 --useAnnotations false

QC thresholds (defaults): `--minReads 10000 --minMapRate 50`

## Other scripts (`scripts/`)
- `align.sh` – Task 3 alignment and BAM conversion
- `run_bambu.R` – Task 4 Bambu run
- `analysis_task4.R` – Task 4.2 answers and 4.4 candidates
- `plot_task4.R` – plotBambu figures (4.4)

## Results (`results/`)
Nextflow reports, QC summary, Task 4 results and plots.

Data (fastq, reference, BAM) are not included; they are available from the SG-NEx project (https://github.com/GoekeLab/sg-nex-data).
