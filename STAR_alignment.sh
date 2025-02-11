#!/bin/bash

######################################################
# STAR ALIGNMENT
#
# GOAL: Use STAR to align barcode and UMI tagged 
# read 2 data to reference genome. 
######################################################

set -e

manifestFile=$1
source "${manifestFile}"

echo ""
echo "Aligning ${writeDir}${batchName}_r2_barcode_tagged.fastq to STAR index at ${STARindex}"
STAR --runThreadN "${nCores}" \
     --genomeDir "${STARindex}" \
     --readFilesIn "${writeDir}${batchName}_r2_barcode_tagged.fastq" \
     --outFileNamePrefix="${batchName}_" \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate;
      
if [ $doNonSyrah = true ]; then
  echo ""
  echo "Aligning ${writeDir}${batchName}_r2_barcode_tagged_nonSyrah.fastq.gz to STAR index at ${STARindex}"
  STAR --runThreadN "${nCores}" \
       --readFilesCommand zcat \
       --genomeDir "${STARindex}" \
       --readFilesIn "${writeDir}${batchName}_r2_barcode_tagged_nonSyrah.fastq.gz" \
       --outFileNamePrefix="${batchName}_nonSyrah_" \
       --outFilterMultimapNmax 1 \
       --outSAMtype BAM SortedByCoordinate;
fi