# slide-seq_pipeline_modification
Modification for the slide-seq v2 analysis pipeline to correct for deletions in bead oligos.

# Required R packages

NOTE: Some of these may not actually be required and are relics from development stage. Let me know if you find a missing or unnecessary one!

* library(readr)
* library(reader)
* library(parallel)
* library(doParallel)
* library(data.table)
* library(foreach)
* library(feather)
* library(parallelDist)

# Instructions

This addition to the slide-seq v2 pipeline will allow you to create a dataset that has (attempted to) correct for synthesis errors in the bead capture oligos. NOTE: If you're still waiting on data, you can still get the first step of **Prep data** done while you wait.

# Gather data

* Combined (all lanes) unaligned read 1 SAM
* Aligned read 2 SAM
* Puck coordinates CSV
* Read 1 pattern

  + JJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJTCNNNNNNNNT for 44 nt read 1
  + JJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJNNNNNNNVVT for 42 nt read 1

# Prep data

* Run `generate_barcode_mapping.R` with PARAMS section customized to your situation.
  (I keep this separate from the main analysis, since the resulting feather file will work for other data from this puck, re-running the pipeline, etc.)
* Create headerless read 1 SAM 
   + sorted by query name
* Create headerless read 2 SAM 
   + Only mapped reads
   + sorted by query name
   + without ZP or ZS tags
* Create header-only read 2 SAM
  + no @PG lines

# Barcode + UMI correction

* Run `generate_corrected_aligned_SAM.R` with PARAMS section customized to your situation
* NOTE: Because read 2 SAM only has successfully aligned reads, there is a discrepancy between the reads in read 1 SAM and read 2 SAM. As chunks of reads are getting read in for processing, we need to read in a bigger chunk of read 1 to ensure we have all read 2 reads in there. The parameter **r1cushion** is a coefficient for how much bigger the read 1 chunk is, with half of the extra cushion effectively on each side of the read 2 chunk. You can either:
  + leave as is
  + decrease to improve runtime
  + increase if you're getting messages about missing reads
  + (You could also try using a read 2 SAM *not* filtered to mapped reads and setting this to 1. I haven't tried it, so no promises.)

# Re-assemble data

* `cat` together no-@PG read 2 header with corrected SAM
* Convert to BAM
* NOTE: Barcode mapping has already been performed, do not repeat.
* Use `DigitalExpression` to convert to DGE matrix, then proceed as usual
