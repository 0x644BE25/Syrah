#!/bin/bash

######################################################
# QUANTIFY COUNTS
#
# GOAL: Use subread to quantify counts, then sort and
# generate gene (row) x bead (col) UMI count matrices.
######################################################

set -e

manifestFile=$1
source "$manifestFile"

# this is the step that requires subread
featureCounts -a "${gtf}" \
              -o "${batchName}_gene_assigned" \
              -R BAM "${batchName}_Aligned.sortedByCoord.out.bam" \
              -T "${nCores}";
              
samtools sort "${batchName}_Aligned.sortedByCoord.out.bam.featureCounts.bam" -o "${batchName}_assigned_sorted.bam"
samtools index "${batchName}_assigned_sorted.bam"
umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I "${batchName}_assigned_sorted.bam" -S "${writeDir}${batchName}_counts.tsv.gz"

              
if [ $doNonSyrah = true ]; then
  featureCounts -a "${gtf}" \
                -o "${batchName}_nonSyrah_gene_assigned" \
                -R BAM "${batchName}_nonSyrah_Aligned.sortedByCoord.out.bam" \
                -T "${nCores}";
                
  samtools sort "${batchName}_nonSyrah_Aligned.sortedByCoord.out.bam.featureCounts.bam" -o "${batchName}_nonSyrah_assigned_sorted.bam"
  samtools index "${batchName}_nonSyrah_assigned_sorted.bam"
  umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I "${batchName}_nonSyrah_assigned_sorted.bam" -S "${writeDir}${batchName}_nonSyrah_counts.tsv.gz"
fi
