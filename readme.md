---
output:
  pdf_document: default
  html_document: default
---
# Syrah: an augmented Slide-seq & Curio pipeline

preprint info and link here, maybe one summary figure

### PURPOSE

This pipeline is intended to improve both the quantity and fidelity of usable data from Slide-seqV2 or Curio spatial transcriptomic datasets by correcting for two sources of error. We also hope that you find it easy to install and run.

## 1. Installation

### Method A: Singularity

Singularity is a container platform which allows you to build an encapsulated software environment ("container") that has all software components installed with the proper versions and configurations. If you do not have Singluarity installed, first [follow the installation instructions](https://docs.sylabs.io/guides/main/user-guide/quick_start.html#quick-installation-steps)

Once Singularity is installed and working, download the pre-built container from [wherever it ends up]()

### Method B: Manual

Alternatively, you can install each of the necessary components yourself. You should install the following software (versions in parenthesis are those used during pipeline development, other versions are likely to work as well):

-   [R](https://www.r-project.org/) (v4.3.1) and the R package [dbscan](https://github.com/mhahsler/dbscan)
-   [Samtools](https://www.htslib.org/download (v1.2)
-   [UMI-tools](https://umi-tools.readthedocs.io/en/latest/INSTALL.html) (v1.1.5)
-   [STAR aligner](https://github.com/alexdobin/STAR) (v2.7.10b)
-   [Subread](https://subread.sourceforge.net/) (v2.0.2)

## 2. Input data

We've tried to choose a minimal set of starting data that uses the most frequently published "raw" files. However, Syrah's barcode correction step does require a read 1 FASTQ that has not yet had the barcode extracted. If your read 1 FASTQ has reads of 42-44 nt long, you should be good to go. The total input data is:

- **Read 1 and read 2 FASTQ files:** These can be either zipped (.fastq.gz) or unzipped (.fastq). Read 1 and read 2 FASTQs must be in the same order, but this is generally how they are made. If you have multiple lanes in your sample, make sure to concatenate all the read 1 FASTQs together and all the read 2 FASTQs together so that you begin with two FASTQ files. If you are starting with BCL files instead, convert them to FASTQs with something like Illumina's [blc2fastq](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf).
- **Bead barcode + coordinates file:** This should have three tab-delimited columns (barcode, x, y) with no headers.
- **Reference genome:** You will need a path to both a STAR index (will be a directory) and a matching .gtf file.
- **Manifest file:** This contains the parameters for the pipeline. If you have multiple samples to process, create a manifest for each and then you can give all the paths as input to the script (i.e. `/path/to/myWorkingDir/*_manifest.txt`). The parameters in the manifest are as follows:

| Parameter name | Description | Example |
| -------- | ----------------- | --------------- |
| `resume` | If restarting the pipeline after a failure, should Syrah look for previous files and attempt to pick up where the last run left off? Must be `true` or `false`, defaults to `false` if no value provided. | `true` |
| `syrahDir` | Path to directory where Syrah code is located. Make sure the directory has the proper permissions for code to be executed from it! | `"/path/to/syrah/code/"` |
| `read1fastq` | Path to read 1 FASTQ file, zipped or not. Read 1 is the read that contains barcode and UMI information, probably 42-44 nucleotides long. | `"/path/to/myWorkingDir/r1.fastq"` |
| `read2fastq` | Path to read 2 FASTQ file, zipped or not. Read 2 contains the "biological" sequence data, probably 50-100 nucleotides long. | `"/path/to/myWorkingDir/r2.fastq"` |
| `puckFile` | Path to a tab-delimited file with **no header** containing three columns: bead barcodes, x coordinates, y coordinates. | `"/path/to/myWorkingDir/coordinates.txt"` |
| `STARindex` | Path to a STAR aligner index directory for the desired refrence genome or transcriptome. The [STAR manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf) has information on how to index a reference in the "Generating genome indexes" section. | `"/my/genomes/Gal_gal/STAR/"` |
| `gtf`  | Path to a GTF file corresponding to the `STARindex` | `"/my/genomes/Gal_gal/GRCg7b.Ens_110.gtf"`
| `writeDir` | Path to an **existing** directory where Syrah should write files. Make sure the directory has the proper permissions for Syrah to write to it! | | `"/path/to/myWorkDir/"` |
| `batchName` | A name to be used for this batch of data. | `"test_data"` |
| `nCores` | Maximum number of CPU cores to be used. If in doubt, use 1 (it'll be slow, though). | `1` |
| `minUMI` | Minimum UMI threshold for including beads in the final data. This only matters for the optional final "Seurat spatial" step. | `10` |
| `maxLinkerDistance` | Maximum allowd linker distance (comprised of deletions and 0 or 1 substitutions). Reads that fail to meet this criteria are discarded. A value of 5 is recommended based on [link to that part of the preprint] | `5` |
| `doNonSyrah` | Should the pipeline also use a generic process to output data comparable to that of the standard Curio Seeker or Slide-seqV2 pipeline? This will take additional time, but will allow for the comparison of results. This should be either `true` or `false`. | `true` |
| `nonSyrahDist2` | If doing the "non-Syrah" version, should we allow for distance=2 barcode matches (Curio-seeker does this)? Defaults to `false` if no value provided. | `false` |


Use the `manifest.txt` file as a template, pay attention to local vs global paths, and ensure that strings are quoted.

## 3. Running Syrah

As with installation, there are two ways to run Syrah: an automated way and a manual method. The manual method is helpful while troubleshooting, while the automated method allows for high-throughput processing of many datasets.
