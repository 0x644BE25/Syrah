#!/bin/bash

# ================= SETUP =======================

source MANIFEST.txt
echo "SYRAH PROCESSING LOG"
echo ""
echo "BAM directory: $BAMdir"
echo "Puck file: $puckFile"
echo "Write directory: $writeDir"
echo "Batch name: $batchName"
echo "Read 1 oligo version: $vs"
echo "Max CPU to use: $nCores"
echo "Minimum nUMI threshold: $minUMI"
echo "Maximum linker alignment distance: $maxLinkDist"
echo ""

mkdir -p "$writeDir""intermediate_files"


# ================= CREATE BEAD BARCODE MAP ==========

echo "Building barcode matching/de-forking map"
Rscript ./create_bead_barcode_mapping.R $puckFile $writeDir $batchName $nCores
echo ""

# ================= PREP READ 1 ======================

r1BAMs=$(ls $BAMdir*unmapped.bam)
echo "Merging read 1 BAM files:"
echo "$r1BAMs"
echo ""

n1Reads <- sum(sapply(r1Bams,function(x){ as.integer(system(paste0('samtools view -c ',bamDir,x),intern=TRUE)) }))
# R1 = UNMAPPED, NEED R1 42 NT SEQUENCE FOR BC/UMI RE-MAPPING
system(paste0('samtools merge ',newDir,batch,'_r1_merged.bam ',paste(paste0(bamDir,r1Bams),collapse=' ')))

read1count=344947432
read1count=0
for bam in $r1BAMs
do
  n=$(samtools view -c $bam)
  read1count=$((read1count+n))
done

samtools merge "${writeDir}intermediate_files/${batchName}_r1_merged.bam" $r1BAMs 
read1mergedCount=$(samtools view -c "${writeDir}intermediate_files/${batchName}_r1_merged.bam")
if [ $read1mergedCount != $read1count ]; then
  echo "ERROR: expected $read1count reads in merged read 1 file, got $read1mergedCount instead."
  echo "This happens sometimes, I'm not sure why. Try this step again."
else
  echo "Read 1 merged, $read1mergedCount total reads."
fi
samtools view -f 77 -o "${writeDir}intermediate_files/${batchName}_r1_merged_filtered.bam" "${writeDir}intermediate_files/${batchName}_r1_merged.bam"
samtools sort -n -o "${writeDir}intermediate_files/${batchName}_r1_merged_filtered_sorted.bam" "${writeDir}intermediate_files/${batchName}_r1_merged_filtered.bam"
echo "Read 1 filtered and sorted"

# ================ PREP READ 2 =======================

r2BAM=$(ls $BAMdir*bam | grep -v 'unmapped')
echo ""
echo "Prepping read 2 BAM file:"
echo "$r2BAM"
echo ""

samtools view -h -F 4 -O bam -o "${writeDir}intermediate_files/${batchName}_r2_filtered.bam" $r2BAM
samtools sort -n -O bam -o "${writeDir}intermediate_files/${batchName}_r2_filtered_sorted.bam" "${writeDir}intermediate_files/${batchName}_r2_filtered.bam"
samtools view -H "${writeDir}intermediate_files/${batchName}_r2_filtered_sorted.bam" | grep -v ^@PG > "${writeDir}intermediate_files/${batchName}_r2_header.bam"
samtools view "${writeDir}intermediate_files/${batchName}_r2_filtered_sorted.bam" > "${writeDir}intermediate_files/${batchName}_r2_reads_only.sam" 
cat "${writeDir}intermediate_files/${batchName}_r2_reads_only.sam" | cut -f1 > "${writeDir}intermediate_files/${batchName}_r2_qnames.txt" 

# ================= QNAME FILTER READ 1 ==============

echo "Filtering read 1 to by sequnce IDs present in read 2"
samtools view -N "${writeDir}intermediate_files/${batchName}_r2_qnames.txt" "${writeDir}intermediate_files/${batchName}_r1_merged_filtered_sorted.bam" | cut -f1,10 > "${writeDir}intermediate_files/${batchName}_r1_qname_filtered_qname_seq_only.txt"

# ================= MAKE CORRECTED BAM ===============

echo "correcting bead barcode (XC) and UMI (XM)"
Rscript ./correct_XC_XM.R $writeDir $batchName $vs $maxLinkDist $nCores

echo "geting read fates"
cat ${writeDir}intermediate_files/${batchName}_XC_XM_corrected.txt | cut -f 20 | sed 's/XX:Z://' > ${writeDir}${batchName}Syrah_test_fates.txt

echo "filter by linker alignment quality""
linkRange=$(seq -s '' 0 $maxLinkDist)
cat ${writeDir}intermediate_files/${batchName}_XC_XM_corrected.txt | grep "XX:Z:[$linkRange]d[89]s" | grep -v "_noBC" > ${writeDir}intermediate_files/${batchName}_XC_XM_corrected_filtered.txt

echo "add SAM headers"
cat ${writeDir}intermediate_files/${batchName}_r2_header.bam ${writeDir}intermediate_files/${batchName}_XC_XM_corrected_filtered.txt > ${writeDir}intermediate_files/${batchName}_XC_XM_corrected.sam
echo "converting to BAM"
samtools view -O bam -o ${writeDir}${batchName}_corrected.bam ${writeDir}intermediate_files/${batchName}_XC_XM_corrected.sam

echo "You're done!"
echo "Just use your preferred method (Picard's DigitalExpression works well) to convert"
echo "${writeDir}${batchName}_corrected.bam to a Digital Gene Expression Matrix"
echo "and get on to the fun part!"