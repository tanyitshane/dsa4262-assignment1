#!/bin/bash
set -euo pipefail
set -x
cd ~/assignment
REF=reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
T=$(nproc)
for fq in fastq/*.fastq.gz; do
  s=$(basename "$fq" .fastq.gz)
  if [[ $s == *directRNA* ]]; then FLAGS="-ax splice -uf -k14"; else FLAGS="-ax splice"; fi
  minimap2 $FLAGS -t $T $REF "$fq" > bam/$s.sam
  samtools view -b -@ $T bam/$s.sam | samtools sort -@ $T -o bam/$s.bam
  samtools index bam/$s.bam
  rm bam/$s.sam
done
