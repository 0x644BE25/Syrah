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
if ! [[ $writeDir == */ ]] ; then
  writeDir="${writeDir}/"
fi
if ! [[ $syrahDir == */ ]] ; then
  syrahDir="${syrahDir}/"
fi


echo ""
echo "Aligning ${writeDir}intermediate_files/${batchName}_r2_barcode_tagged.fastq to STAR index at ${STARindex}"
STAR --runThreadN "${nCores}" \
     --genomeDir "${STARindex}" \
     --readFilesIn "${writeDir}intermediate_files/${batchName}_r2_barcode_tagged.fastq" \
     --outFileNamePrefix="${writeDir}intermediate_files/${batchName}_" \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate;
      
if [ $doNonSyrah = true ]; then
  echo ""
  echo "Aligning ${writeDir}intermediate_files/${batchName}_r2_barcode_tagged_nonSyrah.fastq to STAR index at ${STARindex}"
  STAR --runThreadN "${nCores}" \
       --genomeDir "${STARindex}" \
       --readFilesIn "${writeDir}intermediate_files/${batchName}_r2_barcode_tagged_nonSyrah.fastq" \
       --outFileNamePrefix="${writeDir}intermediate_files/${batchName}_nonSyrah_" \
       --outFilterMultimapNmax 1 \
       --outSAMtype BAM SortedByCoordinate;
fi
