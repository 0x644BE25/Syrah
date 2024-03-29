# Syrah: Slide-seqV2 pipeline augmentation
Modification for the slide-seqV2 analysis pipeline to correct for deletions in bead oligos and bead forking.

Correspondence to: cbrewster@stowers.org

## Required software

* [**Samtools**](https://www.htslib.org/)
* [**R**](https://www.r-project.org/) with the package [groundhog](https://cran.r-project.org/web/packages/groundhog/index.html) installed. Groundhog will then manage the rest of the required packages/versions for you, which include:
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
 
NOTE that groundhog will manage installation and loading of these packages to ensure version compatibility **without** affecting your currently installed versions of any packages. Learn more about Groundhog [here](https://groundhogr.com/).
  
## Required data

### Puck info

This should be a `.txt` or `.csv` with bead sequences and spatial coordinates provided with the puck. 

### Unaligned BAM files

These should all be in the same directory and end in `unmapped.bam`.

### Aligned BAM file

This should be read 2 only and in the same directory as the unaligned BAMs, ending in `.bam` and NOT `unmapped.bam`.

## Process
0. If you want to do several batches, prepare several manifold files (step 1) and pass them to `Syrah.sh` with something like `bash Syrah.sh manifold*txt` to loop over your batches. All files (intermediate or otherwise) are prefixed with the batch name, so MAKE SURE ALL BATCH NAMES ARE UNIQUE!
1. Open `MANIFOLD.txt` and fill out the following information:
    * BAMdir: path to the directory holing the read 1 and read 2 BAM files
    * puckFile: path to the puck info file
    * writeDir: path to the directory where `Syrah` should write files
    * batchName: a name for this batch of data
    * vs: version of bead oligo used in this data, either "vs1" or "vs2 (see Stickels et. al [Slide-seqV2 paper](https://pubmed.ncbi.nlm.nih.gov/33288904/)). You can check this using [queSyrah](https://simrcompbio.shinyapps.io/queSyrah/).
    * nCores: max number of CPU cores to use
    * maxLinkDist: maximum acceptable linker alignment distance, default=5
    * keepIntermed: whether to retain intermediate files (generally for troubleshooting), default=true
2. Run `bash Syrah.sh MANIFOLD.txt` (make sure `Syrah.sh` has [proper executable permissions](https://bash.cyberciti.biz/guide/Setting_up_permissions_on_a_script)). This is going to take **a while**, particularly when using few cores, so use of `nohup` is recommended.
3. You're done! Use the merged BAM to create a DGE matrix and get on to the downstream analysis.

NOTE: All analysis steps are contained in `Syrah.sh`, but running individual commands within the script will allow you to perform different steps as the data becomes available. You can generate a barcode matching/deforking map before even using the puck, and the percentage of barcodes in forks is a rough indicator of puck quality (low % = higher quality puck, high % = lower quality puck, usable range is ~35%, okay range is closer to 20%). You can also merge/filter/sort the read 1 files while waiting on alignment.

## Thanks!

* For trying this out
* For finding the inevitable bugs
* And for letting me know how everything goes and what I can do to make this more useful :)
* Correspondence to: cbrewster@stowers.org
