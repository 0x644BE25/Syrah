#!/bin/bash

######################################################
# EXTRACT BEAD BARCODES NONSYRAH
#
# GOAL: Use UMI-tools and read 1 pattern determined
# by determine_version.R to extract barcode and UMI
# and append to read 2 sequence ID for the non-Syrah
# version of the data.
######################################################

set -e

manifestFile=$1
source "$manifestFile"

r1pattern=$(cat "${writeDir}${batchName}_r1_pattern.txt")
echo "Using read1 static pattern ${r1pattern} for non-Syrah barcode extraction"
                  
  umi_tools extract --extract-method=string \
                    --bc-pattern="${r1pattern}" \
                    --error-correct-cell \
                    --stdin="${read1fastq}" \
                    --stdout="${writeDir}${batchName}_r1_barcode_tagged_nonSyrah.fastq.gz" \
                    --read2-in="${read2fastq}" \
                    --read2-out="${writeDir}${batchName}_r2_barcode_tagged_nonSyrah.fastq.gz" \
                    --whitelist="${writeDir}${batchName}_barcode_whitelist_nonSyrah.txt";