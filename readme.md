# Syrah
Syrah is an R package that provides read 1 error correction for [Slide-seqV2](https://www.nature.com/articles/s41587-020-0739-1) and [Curio Seeker](https://curiobioscience.com/seeker/) spatial transcriptomic data. 

**NOTE:** This is the R package version of Syrah that only corrects the read 1 FASTQ and is intended for use in conjunction with a standard [Slide-seqV2](https://broadinstitute.github.io/warp/docs/Pipelines/SlideSeq_Pipeline/README) or Curio Seeker analysis pipeline. If you are looking for the standalone analysis pipeline, it is [here](https://cran.r-project.org/).

## Tutorial

### Prerequisites

If you do not have [R](https://r-project.org/) installed, simply use the [R installer](https://rstudio-education.github.io/hopr/starting.html) for your operating system: [Linux](https://cran.r-project.org/bin/linux/), [macOS](https://cran.r-project.org/bin/macosx/), or [Windows](https://cran.r-project.org/bin/windows/).

Once R is installed, use the icon to open an R session and install the [devtools](https://cran.r-project.org/package=devtools) library by running

```install.packages("devtools")```

### Installing Syrah
In your R session, install Syrah with the command

```devtools::install_github("0x644BE25/Syrah",force=TRUE)```

and load the library with

```library("Syrah")```

### Downloading test data
We have [two test datasets available](#test-data). To download and unpack the Slide-seqV2 test data, run

```
download.file("https://github.com/0x644BE25/Syrah/raw/refs/heads/main/test_data/planarian_test_data.tar.gz","SlideSeq.tar.gz")
untar("SlideSeq.tar.gz")
```

### Running Syrah
To run Syrah on the test data, use the `syrah()` function and pass it the name of your coordinate/barcode and Read 1 FASTQ files:

```
syrah(coords_file="planarian_L43430_coordinates.txt",r1_fastq="planarian_L43430_H5F73AFX2_1M_R1.fastq.gz")
```

**Youre done!** Once Syrah has finished running, simply use the corrected Read 1 FASTQ along with the original Read 2 FASTQ and barcode/coordinates file in your preferred Slide-seq or Curio Seeker analysis pipeline.

&nbsp;
&nbsp;
***
&nbsp;
&nbsp;
&nbsp;

## Installation
You will need to have [R](https://www.r-project.org/) and [devtools](https://www.rdocumentation.org/packages/devtools/versions/2.4.5) installed, then run
```
library(devtools)
install_github("0x644BE25/Syrah")
```

## Inputs
Syrah takes as input the tab-delimited bead coordinates file and the read 1 FASTQ file.

## Basic Usage
```
library(Syrah)
syrah(coords_file="coodinates_file.tsv", r1_fastq="read1_file.fastq")
```

## Output
Syrah will output a bead deduplication map text file, a barcode whitelist file, and a corrected read 1 FASTQ with a filename ending in `.r1syrah`
This corrected read 1 FASTQ can now be input into your spatial transcriptomic pipeline of choice along with the original read 2 FASTQ and coordinates file.


<details>
  
<summary>Advanced Usage</summary>

## Advanced Usage
Syrah has three steps, and the `Syrah()` function is simply a wrapper for them. The steps can be run independently, if desired. This may be useful if you wish to use the same barcode whitelist for several read 1 FASTQs from the same puck or tile, such as if you have multiple lanes of the same library on a flowcell.

#### Step 1: Barcode deduplication
This step uses the barcode file to find beads with impossibly close x,y distance and barcodes one nucleotide apart. These groups of beads are virtual duplications, so Syrah reroutes all barcodes in the group to a single one.
```
make_bead_dedup_map(coords_file="coodinates_file.tsv")
```
This will output a text file of the deduplication mapping with the same name as the coordinates file but with `_dedup_map.txt` appended.

#### Step 2: Barcode whitelist generation
This step uses the barcode file and deduplication map (from the previous step) to generate a whitelist of all acceptable barcode matches. This whitelist automatically redirects duplicated beads when matching.
```
make_barcode_whitelist(dedup_map="coordinates_file.tsv_dedup_map.txt", coods_file="coordinates_file.tsv")
```
This will output a text file of the barcode matching whitelist with the same name as the coordinates file but with `_whitelist.txt` appended.

#### Step 3: Barcode correction
This step uses the barcode whitelist (from the previous step) and the read 1 FASTQ to generate a corrected read 1 FASTQ.
```
correct_barcodes(whitelist="coordinates_file.tsv_whitelist.txt", r1_fastq="read1_file.fastq")
```
This will output a corrected read 1 FASTQ with the same name as the original read 1 FASTQ but with `.r1syrah` appended. This FASTQ has the same reads in the same order as the original read 1 FASTQ, such that it is still a proper pair with the original read 2 FASTQ. It can now be used as input to your analysis pipeline of choice, such as the [WARP Slide-seq pipeline](https://broadinstitute.github.io/warp/docs/Pipelines/SlideSeq_Pipeline/README) which is also available for use in the cloud via [Terra](https://app.terra.bio/) at the [Slide-seq public workspace](https://app.terra.bio/#workspaces/warp-pipelines/Slide-seq). 

</details>

<details>
  
  <summary> Parameters </summary>
  
## Parameters

| Parameter name      | Description                                                                                                                                                                                                                                                                    | Example                                    |
|------------------|----------------------------|--------------------------|
| `coords_file`         | Path to a tab-delimited file containing the barcodes and coordinates for the puck | `"/path/to/myWorkingDir/coordinates.txt"` |
| `r1_fastq`          | Path to the read 1 FASTQ file, or a comma-delimited list of read 1 FASTQ files | `"path/to/myWorkingDir/read_1.fastq"` |
| `write_dir`        | (optional) Directory to write to. Defaults to current working directory. | `"/path/to/myWorkingDir/"`  |
| `n_cores`        | (optional) Number of CPU cores to use. Defaults to 1.  | `1`         |
| `max_slide_dist`          | (optional) Maxium allowable slide distance between beads to consider them duplicated. You are very unlikely to need to change this. Defaults to 10. | `10`  |
| `max_linker_dels`          | (optional) Maximum allowable number of deletions for an acceptable linker match. You are unlikely to need to change this. Defaults to 5. | `5`  |
| `batch_size`          | (optional) Number of reads to process at once.  You are very unlikely to need to change this. Defaults to 10^5. | `10^5`  | 

</details>

<details>
<summary> Test data</summary>
  
## [Test data](#test-data)

The [Curio Seeker test data](https://github.com/0x644BE25/Syrah/raw/refs/heads/main/test_data/example_input_mouse_spleen_1M.tar.gz) and [Slide-seqV2 test data](https://github.com/0x644BE25/Syrah/raw/refs/heads/main/test_data/planarian_test_data.tar.gz) archive files contain FASTQ files downsampled to 1 million reads and the corresponding barcode/coordinates file. Simply unzip and un-tar the files. [Test data README](https://github.com/0x644BE25/Syrah/blob/main/test_data/README.txt)

The Curio Seeker test data is a mirror of the official [Curio Seeker test data](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/TestDatasets/example_input_mouse_spleen_1M.tar.gz) provided for convenience and continuity.

</details>
