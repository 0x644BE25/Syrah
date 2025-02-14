Syrah: A standalone analysis pipeline to maximize Slide-seq & Curio Seeker data
==============

**Author:** [Carolyn Brewster](mailto:cbrewster@stowers.org)

preprint info and link here, maybe one summary figure

### PURPOSE

This pipeline is intended to improve both the quantity and fidelity of usable data from Slide-seqV2 or Curio spatial transcriptomic datasets by correcting for two sources of error. We also hope that you find it easy to install and run.

## 1. Installation

### Method A: Singularity

Singularity is a container platform which allows you to build an encapsulated software environment ("container") that has all software components installed with the proper versions and configurations. If you do not already have Singularity (if working on a shared compute resource you should check), you will need to install it first.

-   [Linux installation instructions](https://apptainer.org/user-docs/master/quick_start.html#quick-installation-steps)
-   [MacOS installation instructions](https://docs.sylabs.io/guides/3.5/admin-guide/installation.html#macOS)
-   [Windows installation instructions](https://docs.sylabs.io/guides/3.5/admin-guide/installation.html#windows) are a bit more involved due to the need for a virtual machine.

Once Singularity is installed and working, download the pre-built container from [the Syrah v2.0.0-alpha release](https://github.com/0x644BE25/Syrah/releases/download/v2.0.0-alpha/Syrah-v2.0.0.sif)

**NOTE:** Scripts you run using Singularity will always have access to files in your user folder, but if you are reading or writing data in a different folder, you may need to [enable access to it using the `--bind` command](https://apptainer.org/user-docs/master/quick_start.html#working-with-files).

### Method B: Manual

Other than the programming language R, Syrah relies on several common bioinformatics tools. You will need the following software (if you are working on a shared compute resource, you should check if you already have some or all of them). Versions in parenthesis are those used during pipeline development, but most recent versions are likely to work as well:

-   [R](https://www.r-project.org/) (v4.4.1) and the R package [dbscan](https://github.com/mhahsler/dbscan)

    **NOTE:** if incorporating Syrah's barcode correction into your own pipeline, this is the only dependency. `Syrah_minimal.sh` will output a read 2 FASTQ with corrected bead barcodes and UMIs appended to the sequence IDs. If that's all you need, you're good just with R and dbscan.

-   [Samtools](https://www.htslib.org/download) (v1.20)

-   [STAR aligner](https://github.com/alexdobin/STAR) (v2.7.10b)

-   [UMI-tools](https://umi-tools.readthedocs.io/en/latest/INSTALL.html) (v1.1.5)

-   [Subread](https://subread.sourceforge.net/) (v2.0.2)

Once the pre-requisites are installed, you're basically done. Just [clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) the Syrah git repository to your desired location. Don't have `git` installed? No problem. Just copy the `.R` and `.sh` files manually to your desired location. Install complete!

## 2. Input files

**NEED TEST DATA?** The Curio Seeker test dataset is small, formatted correctly, and available from <https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/TestDatasets/example_input_mouse_spleen_1M.tar.gz> You'll simply need a mouse genome and GTF for alignment (but not if you just care about `Syrah_minimal.sh`).

We've tried to choose a minimal set of starting data that uses the most frequently published "raw" files. However, Syrah's barcode correction step does require a read 1 FASTQ that has not yet had the barcode extracted. If your read 1 FASTQ has reads of 42-44 nt long, you should be good to go. The total input data is:

-   **Read 1 and read 2 FASTQ files:** These can be either zipped (.fastq.gz) or unzipped (.fastq). Read 1 and read 2 FASTQs must be in the same order, but this is generally how they are made. If you have multiple lanes in your sample, make sure to concatenate all the read 1 FASTQs together and all the read 2 FASTQs together so that you begin with two FASTQ files. If you are starting with BCL files instead, convert them to FASTQs with something like Illumina's [blc2fastq](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf).
-   **Bead barcode + coordinates file:** This should have three tab-delimited columns (barcode, x, y) with no headers.
-   **Reference genome:** You will need a path to both a STAR index (will be a directory) and a matching .gtf file (not required for `Syrah_minimal.sh`).
-   **Manifest file:** This contains the parameters for the pipeline. If you have multiple samples to process, create a manifest for each and then you can give all the paths as input to the script (i.e. `/path/to/myWorkingDir/*_manifest.txt`). The parameters in the manifest are as follows:

| Parameter name      | Description                                                                                                                                                                                                                                                                    | Example                                    |
|------------------|----------------------------|--------------------------|
| `resume`            | If restarting the pipeline after a failure, should Syrah look for previous files and attempt to pick up where the last run left off? Must be `true` or `false`, defaults to `false` if no value provided.                                                                      | `true`                                     |
| `syrahDir`          | Path to directory where Syrah code is located. Make sure the directory has the proper permissions for code to be executed from it!                                                                                                                                             | `"/path/to/syrah/code/"`                   |
| `read1fastq`        | Path to read 1 FASTQ file, zipped or not. Read 1 is the read that contains barcode and UMI information, probably 42-44 nucleotides long.                                                                                                                                       | `"/path/to/myWorkingDir/r1.fastq"`         |
| `read2fastq`        | Path to read 2 FASTQ file, zipped or not. Read 2 contains the "biological" sequence data, probably 50-100 nucleotides long.                                                                                                                                                    | `"/path/to/myWorkingDir/r2.fastq"`         |
| `puckFile`          | Path to a tab-delimited file with **no header** containing three columns: bead barcodes, x coordinates, y coordinates.                                                                                                                                                         | `"/path/to/myWorkingDir/coordinates.txt"`  |
| `STARindex*`        | Path to a STAR aligner index directory for the desired refrence genome or transcriptome. The [STAR manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf) has information on how to index a reference in the "Generating genome indexes" section. | `"/my/genomes/Gal_gal/STAR/"`              |
| `gtf*`              | Path to a GTF file corresponding to the `STARindex`                                                                                                                                                                                                                            | `"/my/genomes/Gal_gal/GRCg7b.Ens_110.gtf"` |
| `writeDir`          | Path to an **existing** directory where Syrah should write files. Make sure the directory has the proper permissions for Syrah to write to it!                                                                                                                                 |                                            |
| `batchName`         | A name to be used for this batch of data.                                                                                                                                                                                                                                      | `"test_data"`                              |
| `nCores`            | Maximum number of CPU cores to be used. If in doubt, use 1 (it'll be slow, though).                                                                                                                                                                                            | `1`                                        |
| `minUMI`            | Minimum UMI threshold for including beads in the final data. This only matters for the graphical outputs and optional "Seurat spatial" step.                                                                                                                                                   | `10`                                       |
| `maxLinkerDistance` | Maximum allowd linker distance (comprised of deletions and 0 or 1 substitutions). Reads that fail to meet this criteria are discarded. A value of 5 is recommended based on [link to that part of the preprint]                                                                | `5`                                        |
| `doNonSyrah`        | Should the pipeline also use a generic process to output data comparable to that of the standard Curio Seeker or Slide-seqV2 pipeline? This will take additional time, but will allow for the comparison of results. This should be either `true` or `false`.                  | `true`                                     |

`*` = **Not** required for `Syrah_minimal.sh`

Use the `manifest.txt` file as a template, pay attention to local vs global paths, and ensure that strings are quoted.

## 3. Running Syrah

As with installation, there are two ways to run Syrah: an automated way and a manual method. The manual method is helpful for troubleshooting, while the automated method allows for high-throughput processing of many datasets.

### Method A: Automated

All you need to do is pass your manifest file(s) to `Syrah.sh`:

`bash /path/to/syrahDir/Syrah.sh /path/to/my_manifest.txt` for a single manifest, or possibly `bash /path/to/syrahDir/Syrah.sh /path/to/*_manifest.txt` for many manifests.

**SINGULARITY USERS:** use `Singularity exec /path/to/Syrah-v2.0.0.sif /path/to/Syrah.sh /path/to/my_manifest.txt` instead.

### Method B: Manual

I promise it's not as daunting as this schematic makes it look! You just run each command in order, passing your manifest file as the only parameter each time.

![Syrah schematic.](https://github.com/0x644BE25/Syrah/blob/main/simple_Syrah_schematic.png?raw=true "Syrah schematic")

**SINGULARITY USERS:** The best way to run Syrah manually is to start a [Singluarity shell](https://docs.sylabs.io/guides/3.1/user-guide/cli/singularity_shell.html) and then run the commands in it, so start with `singularity shell /path/to/Syrah-v2.0.0.sif` using [`--bind`](https://apptainer.org/user-docs/master/quick_start.html#working-with-files) as needed to ensure access to files not in your home directory. To end the Singularity shell when you're all finished, use `exit`.

**NON-SINGULARITY USERS:** Just skip the above step, but make sure you've [added the install locations for the dependencies to your `$PATH`](https://linuxize.com/post/how-to-add-directory-to-path-in-linux/) so that they can be called from the command line (this should have occurred during installation).

The commands to run each step of the Syrah pipeline in order are

-   1 &ensp; **`Rscript /path/to/syrahDir/determine_version.R /path/to/my_manifest.txt`**

    <details>This step uses the first 100K read 1 sequences to estimate the nucleotide frequency at each position along read 1. Syrah determines which version of capture oligonculeotide was used on the beads on the puck in order to find the expected positions of the both parts of the bead barcode as well as the UMI. There will be a plot in the `intermediate_files` directory showing the nucleotide frequencies and barcode/UMI positions (`B`=bead barcode, `U`=UMI). A nucleotide frequency barplot with colors indicating A, C, G, or T. The x-axis indicates position along read 1, from 5' on the left to 3' on the right, with letters indicating the position of the bead barcode (B) and UMI (U).
       &nbsp; 
       
    **KEY OUTPUT FILES:** `estimated_read_1_nucleotide_frequencies.csv`, `r1_pattern.txt`, `r1_version.txt`</details>

-   2 &ensp; **`Rscript /path/to/syrahDir/create_bead_deduplication_map.R /path/to/my_manifest.txt`**

    <details>This step uses the bead coordinates file to find virtually duplicated beads. These "beads" are less than one bead diameter apart in x,y space and have barcodes that differ by one nucleotide substitution. The groups of virtually duplicated beads will be used by `generate_bead_whitelist.R` to "de-duplicate" the beads so that all reads from that group are assigned to a single bead.
    &nbsp; 
       
    **KEY OUTPUT FILES:** `deduplication_map.txt`</details>

-   3 &ensp; **`Rscript /path/to/syrahDir/generate_bead_barcode_whitelist.R /path/to/my_manifest.txt`**

    <details>This step uses the bead coordinates file and the deduplication map from the previous step to generate valid matches for each bead barcode. Valid matches may have one nucleotide deletion or one substitution. If a match has a deletion that occurs in the first part of the barcode, it will be 13 nucleotides long instead of the usual 14. If `doNonSyarh=true`, a second whitelist will be made without bead deduplication and only allowing one nucleotide substitution (no deletions). This will be used for barcode error correction during the barcode extraction step.
    &nbsp; 
       
    **KEY OUTPUT FILES:** `barcode_whitelist.txt`, `barcode_whitelist_nonSyrah.txt` (if `doNonSyrah=true`)</details>

-   4 &ensp; **`Rscript /path/to/syrahDir/extract_bead_barcodes.R /path/to/my_manifest.txt`**

    <details>This step takes barcode and UMI sequences from read 1 and appends them to the sequence ID in read 2. Syrah uses fuzzy matching to find the location of the invariant linker sequence between barcode parts 1 and 2 and uses this to determine the correct location of the barcode and UMI. The barcode sequence is then matched agains the whitelist from the previous step to correct and deduplicate the bead barcodes before appending the barcode and UMI to the sequence ID of the corresponding read 2. Reads lacking a high-confidnce linker position or valid bead barcode are discarded.
    &nbsp; 
    If `doNonSyrah=true` the non-Syrah version always takes the barcode and UMI from the canonically expected positions and matches agains the non-Syrah whitelist. It also appends the barcode and UMI to the read 2 sequence ID and discards reads without a valid barcode. 
    &nbsp; 
       
    **KEY OUTPUT FILES:** `r2_barcode_tagged.fastq`, `r2_barcode_tagged_nonSyrah.fastq` (if `doNonSyrah=true`)</details>

    **NOTE:** This is the end of the minimal Syrah pipeline, when you reach a barcode corrected and UMI tagged read 2 FASTQ.

-   5 &ensp; **`bash /path/to/syrahDir/STAR_alignment.sh /path/to/my_manifest.txt`**

    <details>This step uses the STAR aligner to align the barcode corrected and tagged read 2 FASTQ to the reference genome or transcriptome. You can modify the alignment parameters in the `STAR_alignment.sh` file.
    &nbsp; 
       
    **KEY OUTPUT FILES:** `Aligned.sortedByCoord.out.bam`, `nonSyrah_Aligned.sortedByCoord.out.bam` (if `doNonSyrah=true`)</details>

-   6 &ensp; **`bash /path/to/syrahDir/quantify_counts.sh /path/to/my_manifest.txt`**

    <details>This step first uses Subread's `featureCounts` function to assign reference-aligned reads to the gene features present in your GTF file. Then the BAM is sorted and indexed so that the UMI-tools function `count` can generate a digital gene expression matrix with genes/transcripts as the rows and beads as the columns.
    &nbsp; 
       
    **KEY OUTPUT FILES:** `counts.tzv.gz`, `nonSyrah_counts.tsv.gz` (if `doNonSyrah=true`)</details>

-   7 &ensp; **`Rscript /path/to/syrahDir/graphical_outputs.R /path/to/my_manifest.txt`**

    <details>This step first uses can generates a PDF with some summary and QA plots for the final data, showing info for eithat all beads in the dataset (>=1 UMIs) or those with at least 25 reads (>=25 UMIs).
       &nbsp; 
       
       **KEY OUTPUT FILES:** `Syrah_results_summary.pdf`, `nonSyrah_results_summary.pdf` (if `doNonSyrah=true`)</details>

And that's -- you're done! You'll have a gene expression matrix called `batchName_counts.txt.gz` and a result summary called `batchName_Syrah_results_summary.pdf` (and ones for the non-Syrah version if `doNonSyrah=true`). There's a directory called `intermediate_files` which contains precisely that and can be deleted if you're sure everything went as planned. Here's an example of what the results summary looks like:
&nbsp; 

![example Syrah results summary](https://github.com/0x644BE25/Syrah/blob/main/example_Syrah_results_summary.png?raw=true)
&nbsp; 


## Questions? Problems? Reach out!

Thanks for trying out Syrah :)

This is an alpha release and there **WILL** be bugs! Don't hesitate to contact [Carolyn](mailto:cbrewster@stowers.org) (or through GitHub) if you have any questions or issues.
