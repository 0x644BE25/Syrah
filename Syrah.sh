#!/bin/bash

######################################################
# SYRAH
#
# GOAL: String together the various bits, checking to
# see if we need to run all or can resume a previous 
# incomplete run.
######################################################

set -e

echo "SYRAH PROCESSING LOG"
i=1;
for manifestFile in "$@"
do
  resume=false
  source "$manifestFile"
  if ! [[ $writeDir == */ ]] ; then
    writeDir="${writeDir}/"
  fi
  if ! [[ $syrahDir == */ ]] ; then
    syrahDir="${syrahDir}/"
  fi

  
  
  echo ""
  echo $(date)
  echo "Manifest ${i}: $manifestFile"
  echo "Syrah code directory: $syrahDir"
  echo "Read 1 FASTQ: $read1fastq"
  echo "Read 2 FASTQ: $read2fastq"
  echo "Puck file: $puckFile"
  echo "Write directory: $writeDir"
  echo "Batch name: $batchName"
  echo "Read 1 oligo version: $read1format"
  echo "Max CPU to use: $nCores"
  echo "Minimum nUMI threshold: $minUMI"
  echo "Maximum linker alignment distance: $maxLinkerDistance"
  i=$((i + 1))
  
  if [ -f "${syrahDir}paths.txt" ]; then
    while IFS="" read -r p || [ -n "$p" ]
    do
      PATH=$p:$PATH
    done < "${syrahDir}"paths.txt
  fi
  
  if [ ! -d "${writeDir}intermediate_files" ]; then
    mkdir "${writeDir}intermediate_files"
  fi
  
  # DETERMINE READ 1 VERSION
  if [ ! -f "${writeDir}intermediate_files/${batchName}_r1_version.txt" ] || [ "$resume" = false ]; then
    resume=false
    Rscript "${syrahDir}determine_version.R" "${manifestFile}"
    echo "Nucleotide frequency plot at ${writeDir}intermediate_files/${batchName}_estimated_read_1_nucleotide_frequencies.png"
  else
    echo ""
    echo "Looks like read1 version already detected! Moving on..."
    echo "  (Delete files from previous runs or use 'resume=false' to change this behavior)"
  fi
  
  # DEFORK MAPPING
  if [ ! -f "${writeDir}intermediate_files/${batchName}_deduplication_map.txt" ] || [ "$resume" = false ]; then
    resume=false
    Rscript "${syrahDir}create_bead_deduplication_map.R" "${manifestFile}"
  else
    echo ""
    echo "Looks like bead deduplication maps were already computed! Moving on..."
    echo "  (Delete files from previous runs or use 'resume=false' to change this behavior)"
  fi

  # BARCODE MATCHING + CORRECTION WHITELIST
  if [ ! -f "${writeDir}intermediate_files/${batchName}_barcode_whitelist.txt"  ] || [ "$resume" = false ]; then
    resume=false
    Rscript "${syrahDir}generate_bead_barcode_whitelist.R" "${manifestFile}"
  else
    echo ""
    echo "Looks like bead barcode whitelists have already been generated! Moving on..."
    echo "  (Delete files from previous runs or use 'resume=false' to change this behavior)"
  fi
  

  # READ BARCODE EXTRACTION + CORRECTION
  if [ ! -f "${writeDir}intermediate_files/${batchName}_r2_barcode_tagged.fastq" ] || [ "$resume" = false ]; then
    resume=false
    Rscript "${syrahDir}extract_bead_barcodes.R" "${manifestFile}"
  else
    echo ""
    echo "Looks like we already added barcodes to read 2! Moving on..."
    echo "  (Delete files from previous runs or use 'resume=false' to change this behavior)"
  fi
  
  # STAR ALIGNMENT
  if [ ! -f "${writeDir}intermediate_files/${batchName}_Aligned.sortedByCoord.out.bam" ] || [ "$resume" = false ]; then
    resume=false
    bash "${syrahDir}STAR_alignment.sh" "${manifestFile}"
  else
    echo ""
    echo "Looks like you've already STAR aligned! Moving on..."
    echo "  (Delete files from previous runs or use 'resume=false' to change this behavior)"
  fi
  
  # READ QUANTIFICATION
  if [ ! -f "${writeDir}${batchName}_counts.tsv.gz" ] || [ "$resume" = false ]; then
    resume=false
    bash "${syrahDir}quantify_counts.sh" "${manifestFile}"
    Rscript "${syrahDir}graphical_outputs.R" "${manifestFile}"
    echo "${batchName} ALL DONE! Find your counts matrix at ${writeDir}${batchName}_counts.tsv.gz"
    echo "and a summary of results at ${writeDir}${batchName}_Syrah_results_summary.pdf"
    if [ "$doNonSyrah" = true ]; then
      echo "and for the non-Syrah version at ${writeDir}${batchName}_nonSyrah_results_summary.pdf"
    fi
  else
    echo ""
    echo "Looks like you've already made your count matrices, so this run has done nothing."
    echo ""
    echo "     (what are we even doing here)"
    echo ""
    echo "Unless this was your intended behavior, you should delete the files from old runs"
    echo "or make sure to set 'resume=false' in the manifest."
  fi

 i=$(($i+1))
 
done
