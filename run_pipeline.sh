# 1. Create barcode map
# 2. Merge unaligned BAMs
# 3. Filter unaligned to R1 only
# 4. BC + UMI correction
# 5. Picard sort by qname
# 6. Picard sort aligned R2 by qname
# 7. Check all qnames match
# 8. Picard MergeBamAlignment

# ================= PARAMS ===========================

# COMPUTING PARAMS
nCores=64

# PATHS
bamDir='<path to BAM file directory>'
puckFile='<path to puck info txt>'
referenceFasta='<path to genome reference FASTA>'
writeDir='<path for writing data>'
tempDir='<path to TEMP directory>'
picardJar='<path to picard.jar>'

bamDir='/n/core/Bioinformatics/tools/slideseq-tools/runs/HNCLCBGXH/libraries/2021-3-27_L45629/'
puckFile='/n/analysis/slideseq/210113_17/puck.txt'
referenceFasta='/n/core/Bioinformatics/tools/slideseq-tools/ref/mm10.Ens_98/mm10.Ens_98.fasta'
writeDir='/home/cb2350/slide-seq/new_pipeline/TEST/'
tempDir='/home/cb2350/slide-seq/TEMP/'
picardJar='/n/apps/CentOS7/bin/picard.jar'

name=`echo "$bamDir" | tr -cs 'a-zA-Z0-9' '\n' | tail -n 1`
puck=`echo "$puckFile" | tr -cs 'a-zA-Z0-9_' '\n' | tail -n 3 | head -n 1`

# ================= PUCK BARCODES ====================

#Rscript generate_barcode_map.R "$nCores" "$puckFile" "$puck" "$writeDir"

# ================= UNALIGNED ========================

# merge flowcell lanes, filter to read 1 only
#samtools merge -n - "$bamDir"*"unmapped.bam"  | samtools view -f 77 -o "$writeDir""$name""_unaligned_r1.bam" -
#printf "\n## unaligned reads merged + filtered to r1 at\n  $writeDir""$name""_unaligned_r1.bam\n\n"

# correct bead barcode and UMI
#Rscript tag_r1_with_XC_XM.R "$nCores" "$writeDir" "$puck" "$name"

# convert SAM to BAM
#samtools view -O bam -o "$writeDir""$name""_unaligned_r1_tagged.bam" "$writeDir""$name""_unaligned_r1_tagged.sam" 
#printf "\n## corrected unaligned SAM converted to BAM\n\n"

# picard sort by queryname
# NOTE: Don't try doing this with samtools -- just use Picard for both BAMs
#java -jar "$picardJar" SortSam I="$writeDir""$name""_unaligned_r1_tagged.bam" O="$writeDir""$name""_unaligned_r1_tagged_sorted.bam" SORT_ORDER=queryname TMP_DIR="$tempDir"
#printf "\n## sorted tagged unaligned BAM at\n  $writeDir""$name""_unaligned_r1_tagged_sorted.bam\n\n"

# get queryname to filter r2
#samtools view "$writeDir""$name""_unaligned_r1_tagged_sorted.bam" | cut -f1 -d$'\t' > "$writeDir""$name""_r1_qnames.txt"

# ================= ALIGNED ==========================

samtools view -N "$writeDir""$name""_r1_qnames.txt" -o "$writeDir""$name""_r2_filtered.bam" "$bamDir""$name"".bam"
#printf "\n## aligned BAM to filtered to barcode-mapped reads at\n  $bamDir""$name"".bam\n\n"

java -jar "$picardJar" SortSam I="$writeDir""$name""_r2_filtered.bam" O="$writeDir""$name""_r2_filtered_sorted.bam" SORT_ORDER=queryname TMP_DIR="$tempDir"
#printf "\n## sorted aligned read 2 BAM at\n  $writeDir""$name""_r2_filtered_sorted.bam\n\n"

# ================= MERGE ============================

java -jar "$picardJar" MergeBamAlignment ALIGNED="$writeDir""$name""_r2_filtered_sorted.bam" UNMAPPED="$writeDir""$name""_unaligned_r1_tagged_sorted.bam" O="$writeDir""$name""_merged.bam" R="$referenceFasta" SORT_ORDER=queryname TMP_DIR="$tempDir"
printf "\n## merged aligned and unaligned $name reads. CONGRATS! Your finished BAM is at\n  $writeDir""$name""_merged.bam\n\n"

printf "\nRECOMMENDED NEXT STEPS: Use Picard's DigitalExpression to convert finished BAM to a DGE matrix and import into your favorite single-cell system.\nDon't forget to add the x,y coordinates from the puck file at $puckFile \n"