# Modified slide-seq pipeline v2
Modification for the slide-seq v2 analysis pipeline to correct for deletions in bead oligos. Now improved for user-friendliness, minimal computation, and maximized parallelism with the standard pipeline.

Correspondence to: cbrewster@stowers.org

## Required software

* [**Samtools**](https://www.htslib.org/)
* [**Picard**](https://broadinstitute.github.io/picard/)
* [**R**](https://www.r-project.org/) with packages
  + data.table
  + doParallel
  + dplyr
  + feather
  + foreach
  + parallel
  + parallelDist
  + readr
  + Rsamtools
  + stringdist
  + stringi
  
## Required data

### Puck info

This should be a `.txt` or `.csv` with bead sequences and spatial coordinates provided with the puck. 

### Unaligned BAM files

These should all be in the same directory and end in `unmapped.bam`. The last uninterrupted alphanumeric chunk of the directory path should be the puck name, i.e. `/n/core/Bioinformatics/tools/slideseq-tools/runs/HNCLCBGXH/libraries/2021-3-27_L45629/` will have the name **L45629**.

### Aligned BAM file

This should be read 2 only and in the same directory as the unaligned BAMs.

### Reference FASTA

This should be the same FASTA that was used for the aligned BAM file.

## Process

All steps are contained in `Syrah.sh`, but running individual commands within the script will allow you to perform different steps as the data becomes available. Here's how it works:

1. **SETUP** Open the `run_pipeline.sh` script and enter the necessary file and pathway info along with the number of cores to use. Setting `nCores` to 1 will simply run non-parallelized. If you want any finer control over distance metrics and whatnot, you will need to edit the PARAMETER section of the relevant R scripts.
2. **PUCK BARCODES** This step only requires the puck data, and the results can be reused for any sequencing done on the same puck.
  2.1 De-fork barcodes and create 1 Hamming or 1 deletion barcode matching lookup table using `generate_barcode_map.R`
3. **UNALIGNED READS** This can be done while waiting on sequence alignment. It is the most time-consuming portion.
  3.1 Merge unaligned BAMs and filter to read 1 only using Samtools.
  3.2 Tag read 1 with bead (XC) and molecular (XM) barcodes using `tag_r1_with_XC_XM.R` This is the time-consuming part.
  3.3 Convert resulting SAM to BAM with Samtools
  3.4 Sort by queryname using Picard (yes, it is important to use Picard and not Samtools for this).
  3.5 Extract a list of querynames to filter aligned reads.
4. **ALIGNED READS** You're almost there!
  4.1 Filter aligned reads using the queryname list from 3.5
  4.2 Sort by queryname using Picard (I haven't tested whether this is strictly necessary, if you do please let me know the result!)
5. **MERGE**
  5.1 Merge the BAMs from 3.4 and 4.2 using Picard.

You're done! Use the merged BAM to create a DGE matrix and get on to the downstream analysis.

## Thanks!

* For trying this out
* For finding the inevitable bugs
* And for letting me know how everything goes and what I can do to make this more useful :)
* Correspondence to: cbrewster@stowers.org
