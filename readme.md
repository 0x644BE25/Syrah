---
output:
  html_document: default
  pdf_document: default
---
# Syrah: an analysis pipeline to maximize Slide-seq & Curio Seeker data

preprint info and link here, maybe one summary figure

### PURPOSE

This pipeline is intended to improve both the quantity and fidelity of usable data from Slide-seqV2 or Curio spatial transcriptomic datasets by correcting for two sources of error. We also hope that you find it easy to install and run.

## 1. Installation

### Method A: Singularity

Singularity is a container platform which allows you to build an encapsulated software environment ("container") that has all software components installed with the proper versions and configurations. If you do not already have Singularity (if working on a shared compute resource you should check), you will need to install it first.

- [Linux installation instructions](https://apptainer.org/user-docs/master/quick_start.html#quick-installation-steps)
- [MacOS installation instructions](https://docs.sylabs.io/guides/3.5/admin-guide/installation.html#macOS)
- [Windows installation instructions](https://docs.sylabs.io/guides/3.5/admin-guide/installation.html#windows) are a bit more involved due to the need for a virtual machine.

Once Singularity is installed and working, download the pre-built container from [the Syrah v2.0.0-alpha release](https://github.com/0x644BE25/Syrah/archive/refs/tags/v2.0.0-alpha.sif)

**NOTE:** Scripts you run using Singularity will always have access to files in your user folder, but if you are reading or writing data in a different folder, you may need to [enable access to it using the  `--bind` command](https://apptainer.org/user-docs/master/quick_start.html#working-with-files).

### Method B: Manual

Other than the programming language R,  Syrah relies on several common bioinformatics tools. You will need the following software (if you are working on a shared compute resource, you should check if you already have some or all of them). Versions in parenthesis are those used during pipeline development, but most recent versions are likely to work as well:

-   [R](https://www.r-project.org/) (v4.4.1) and the R package [dbscan](https://github.com/mhahsler/dbscan)

    **NOTE:** if incorporating Syrah's barcode correction into your own pipeline, this is the only dependency. `Syrah_minimal.sh` will output a read 2 FASTQ with corrected bead barcodes and UMIs appended to the sequence IDs. If that's all you need, you're good just with R and dbscan.

-   [Samtools](https://www.htslib.org/download) (v1.20)
-   [STAR aligner](https://github.com/alexdobin/STAR) (v2.7.10b)
-   [UMI-tools](https://umi-tools.readthedocs.io/en/latest/INSTALL.html) (v1.1.5)
-   [Subread](https://subread.sourceforge.net/) (v2.0.2)

Once the pre-requisites are installed, you're basically done. Just [clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) the Syrah git repository to your desired location. Don't have `git` installed? No problem. Just copy the `.R` and `.sh` files manually to your desired location. Install complete!


## 2. Input files

We've tried to choose a minimal set of starting data that uses the most frequently published "raw" files. However, Syrah's barcode correction step does require a read 1 FASTQ that has not yet had the barcode extracted. If your read 1 FASTQ has reads of 42-44 nt long, you should be good to go. The total input data is:

- **Read 1 and read 2 FASTQ files:** These can be either zipped (.fastq.gz) or unzipped (.fastq). Read 1 and read 2 FASTQs must be in the same order, but this is generally how they are made. If you have multiple lanes in your sample, make sure to concatenate all the read 1 FASTQs together and all the read 2 FASTQs together so that you begin with two FASTQ files. If you are starting with BCL files instead, convert them to FASTQs with something like Illumina's [blc2fastq](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf).
- **Bead barcode + coordinates file:** This should have three tab-delimited columns (barcode, x, y) with no headers.
- **Reference genome:** You will need a path to both a STAR index (will be a directory) and a matching .gtf file (not required for `Syrah_minimal.sh`).
- **Manifest file:** This contains the parameters for the pipeline. If you have multiple samples to process, create a manifest for each and then you can give all the paths as input to the script (i.e. `/path/to/myWorkingDir/*_manifest.txt`). The parameters in the manifest are as follows:

| Parameter name | Description | Example |
| -------- | ----------------- | --------------- |
| `resume` | If restarting the pipeline after a failure, should Syrah look for previous files and attempt to pick up where the last run left off? Must be `true` or `false`, defaults to `false` if no value provided. | `true` |
| `syrahDir` | Path to directory where Syrah code is located. Make sure the directory has the proper permissions for code to be executed from it! | `"/path/to/syrah/code/"` |
| `read1fastq` | Path to read 1 FASTQ file, zipped or not. Read 1 is the read that contains barcode and UMI information, probably 42-44 nucleotides long. | `"/path/to/myWorkingDir/r1.fastq"` |
| `read2fastq` | Path to read 2 FASTQ file, zipped or not. Read 2 contains the "biological" sequence data, probably 50-100 nucleotides long. | `"/path/to/myWorkingDir/r2.fastq"` |
| `puckFile` | Path to a tab-delimited file with **no header** containing three columns: bead barcodes, x coordinates, y coordinates. | `"/path/to/myWorkingDir/coordinates.txt"` |
| `STARindex*` | Path to a STAR aligner index directory for the desired refrence genome or transcriptome. The [STAR manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf) has information on how to index a reference in the "Generating genome indexes" section. | `"/my/genomes/Gal_gal/STAR/"` |
| `gtf*`  | Path to a GTF file corresponding to the `STARindex` | `"/my/genomes/Gal_gal/GRCg7b.Ens_110.gtf"`
| `writeDir` | Path to an **existing** directory where Syrah should write files. Make sure the directory has the proper permissions for Syrah to write to it! | | `"/path/to/myWorkDir/"` |
| `batchName` | A name to be used for this batch of data. | `"test_data"` |
| `nCores` | Maximum number of CPU cores to be used. If in doubt, use 1 (it'll be slow, though). | `1` |
| `minUMI` | Minimum UMI threshold for including beads in the final data. This only matters for the optional final "Seurat spatial" step. | `10` |
| `maxLinkerDistance` | Maximum allowd linker distance (comprised of deletions and 0 or 1 substitutions). Reads that fail to meet this criteria are discarded. A value of 5 is recommended based on [link to that part of the preprint] | `5` |
| `doNonSyrah` | Should the pipeline also use a generic process to output data comparable to that of the standard Curio Seeker or Slide-seqV2 pipeline? This will take additional time, but will allow for the comparison of results. This should be either `true` or `false`. | `true` |

`*` = **Not** required for `Syrah_minimal.sh`

Use the `manifest.txt` file as a template, pay attention to local vs global paths, and ensure that strings are quoted.

## 3. Running Syrah

As with installation, there are two ways to run Syrah: an automated way and a manual method. The manual method is helpful for troubleshooting, while the automated method allows for high-throughput processing of many datasets.

### Method A: Automated

All you need to do is pass your manifest file(s) to `Syrah.sh`:

`bash /path/to/syrahDir/Syrah.sh /path/to/my_manifest.txt` for a single manifest, or possibly
`bash /path/to/syrahDir/Syrah.sh /path/to/*_manifest.txt` for many manifests.

**SINGULARITY USERS:** use `Singularity exec /path/to/Syrah-v2.0.0.sif /path/to/Syrah.sh /path/to/my_manifest.txt` instead.

### Method B: Manual

**SINGULARITY USERS:** The best way to run Syrah manually is to start a [Singluarity shell](https://docs.sylabs.io/guides/3.1/user-guide/cli/singularity_shell.html) and then run the commands in it, so start with `singularity shell /path/to/Syrah-v2.0.0.sif` using [`--bind`](https://apptainer.org/user-docs/master/quick_start.html#working-with-files) as needed to ensure access to files not in your home directory. To end the Singularity shell when you're all finished, use `exit`.

**NON-SINGULARITY USERS:** Just skip the above step, but make sure you've [added the install locations for the dependencies to your `$PATH`](https://linuxize.com/post/how-to-add-directory-to-path-in-linux/) so that they can be called from the command line (this should have occurred during installation).


The commands to run each step of the Syrah in order are

- `Rscript /path/to/syrahDir/determine_version.R /path/to/my_manifest.txt`
- `Rscript /path/to/syrahDir/create_bead_deduplication_map.R /path/to/my_manifest.txt`
- `Rscript /path/to/syrahDir/generate_bead_barcode_whitelist.R /path/to/my_manifest.txt`
- `Rscript /path/to/syrahDir/extract_bead_barcodes.R /path/to/my_manifest.txt`

    **NOTE:** This is the end of the minimal Syrah pipeline, when you reach a barcode corrected and UMI tagged read 2 FASTQ.

- `bash /path/to/syrahDir/STAR_alignment.sh /path/to/my_manifest.txt`
- `bash /path/to/syrahDir/quantify_counts.sh /path/to/my_manifest.txt`
- `Rscript /path/to/syrahDir/graphical_outputs.R /path/to/my_manifest.txt`

And that's -- you're done! You'll have a gene expression matrix called `batchName_counts.txt.gz` and a result summary called `batchName_Syrah_results_summary.pdf` (and ones for the non-Syrah version if `doNonSyrah=true`). There's a directory called `intermediate_files` which contains precisely that and can be deleted if you're sure everything went as planned.

## Questions? Problems? Reach out!

Thanks for trying out Syrah :) 

This is an alpha release and there **WILL** be bugs! Don't hesitate to contact Carolyn at cbrewster@stowers.org (or through GitHub) if you have any questions or issues.
