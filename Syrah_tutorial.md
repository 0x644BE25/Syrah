# Syrah Tutorial

This tutorial will walk you through installing Syrah, building a genome reference, and running the Syrah pipeline on a test dataset provided by Curio. It will require around 3hours (depending on your internet connection), 10 GB RAM, and 20GB disk space, including temporary files.

## 1. Installation (~10 min)

**Installing R:** If you don't have a version of R installed, you will need to install one. Fortunately there are wonderful installers that make it very easy: [find one here](https://www.r-project.org/).

**Installing everything else:** If you're on MacOS or linux, you can install the dependencies using `install_dependencies.sh`. This will require some minimal use of the terminal and will install software with default options. For alternative installation options, see the [README](https://github.com/0x644BE25/Syrah/blob/main/readme.md). Steps:

- 1 &ensp; Download Syrah's files to the directory where you want to install Syrah. <details>You can download a zip file of Syrah from [here](https://github.com/0x644BE25/Syrah/archive/master.zip), unzip it, and move the files from the unzipped `Syrah-main` folder (NOT the entire `Syrah-main` folder) to your install directory.</details>
- 2 &ensp; Open a terminal window and navigate to your chosen install directory <details>**MacOS:** `Applications > Utilities > Terminal.app` and navigate to the folder where you want to install software, OR right-click on the folder where you want to install software and choose `New Terminal at Folder` </details>
- 3 &ensp; Run `sudo bash install_dependencies.sh` **NOTE:** You will need to enter your user password and confirm some steps during installation, so watch the terminal. This process should take only a few mintues depending on your internet connection.

## 2. Download data files (~10 min)

Still in the terminal, we'll download and unzip the test data with the following commands
```
curl -O https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/TestDatasets/example_input_mouse_spleen_1M.tar.gz
tar -zxf example_input_mouse_spleen_1M.tar.gz
```

## 3. Prepare the reference genome (~30 min download, ~ 30 min build STAR index)
We'll need a genome reference to align to, so we will make a directory for our reference genome, move into it, download the FASTA and GTF files, and build our STAR index. Use your terminal to execute the following commands: 

(These downloads are big and building the STAR index takes a while and can be resource intensive. It's best not to do much else with your computer while STAR runs, so this is a great time to go get coffee. If you know you've got more cores available, you could bump up the `STAR` option `--runThreadN`, but if you do increase the threads and the process hangs for >1 hr, go back to 1 thread.)

```
mkdir GRCm39_mouse_genome
cd GRCm39_mouse_genome

curl --http1.0 -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
gunzip GCF_000001635.27_GRCm39_genomic.fna.gz

curl --http1.0 -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz
gunzip GCF_000001635.27_GRCm39_genomic.gtf.gz

# this will tell the terminal where to look for the stuff we installed (like STAR)
while IFS="" read -r p || [ -n "$p" ]; do PATH=$p:$PATH; done < ../paths.txt

STAR --runThreadN 2 \
     --runMode genomeGenerate \
     --genomeDir .\
     --genomeFastaFiles ./GCF_000001635.27_GRCm39_genomic.fna \
     --sjdbGTFfile ./GCF_000001635.27_GRCm39_genomic.gtf \
     --sjdbOverhang 75 \
     --genomeSAsparseD 3 \
     --genomeSAindexNbases 12 \
     --limitGenomeGenerateRAM 14000000000
     
cd ..
```
(In case you're curious, the `\` at the end of a line means that the command isn't finished -- this lets us split a command across multiple lines for convenience or clarity.)

**HELP, IT'S NOT WORKING!** Building a genome reference on a computer with limited hardware can be time consuming and prone to failure. Fortunately Zenodo has a pre-built [GRCm39 STAR 20.7.10b genome](https://zenodo.org/records/10655615) that you can download. Just unzip it and put the contents into the `GRCm39_mouse_genome` directory.

## 4. Fill out manifest

All the information that Syrah needs to run the pipeline is stored in the `manifest.txt` file. The main README has [details about specific manifest parameters](https://github.com/0x644BE25/Syrah/blob/main/readme.md#manifest-parameters) but for this tutorial the file `tutorial_manifest.txt` should have everything in the right place. 

**OPTIONAL:** The unedited `tutorial_manifest.txt` uses _relative_ paths instead of _global_ paths, but when possible it is recommended to use global paths to prevent errors. [Learn about relative versus absolute paths here](https://www.linuxbash.sh/post/understanding-absolute-and-relative-paths). Converting `tutorial_manifest.txt` to global paths is good practice to prepare for running Syrah on your own data! To find the global path to a given directory, navigate your terminal to (or open a terminal window at) the directory and use `pwd` 

## 5. Run Syrah (~60 min)

Using a terminal in your Syrah directory, pass your manifest to Syrah with

```
bash Syrah.sh tutorial_manifest.txt
```
or use this version to also save a copy of the terminal window output to a log file, which can help with troubleshooting
```
bash Syrah.sh tutorial_manifest.txt | tee -a syrah_output_log.txt
```
and that's it! Everything else should be handled by Syrah. If the pipeline is interrupted, you can use the prior command again to re-start after the last successful step. If that's not working right or you want to start from the beginning, delete the `intermediate_files` folder or set `resume=false` in `tutorial_manifest.txt`.

 Your Syrah folder should now have a PDF that will let you know how the analysis went and give you an idea of data quality.

## 5. What's next?

Once Syrah completes you'll have a counts file that can be imported into your favorite analysis pipeline, along with the coordinates file that stores bead location data. Syrah omits these steps to minimize installation requirements. Here are brief instructions for a couple of popular options:

### Seurat (R)

Here we'll create a Seurat object, add a spatial embedding, and do some basic processing with the same parameters used by the Curio-seeker pipeline. We will assume that you are working in an R environment, either from an [R session in your terminal](http://countbio.com/web_pages/left_object/R_for_biology/R_fundamentals/R_sessions.html), or an [IDE](https://www.codecademy.com/article/what-is-an-ide) such as [RStudio](https://posit.co/download/rstudio-desktop/). **NOTE** that all of the following code is **R** code, rather than bash terminal code like the previous sections.

In addition to R, you will need to have the [Seurat](https://satijalab.org/seurat/) library installed. Because it was originally developed for single-cell data, it refers to "cells" rather than "beads".


 - Setup: Once you've started R (RStudio will do this automatically) and use `setwd("/path/to/mySyrahDirectory/")` to make sure you're in the same directory where you ran Syrah. You can tell it's the right directory if it contains the `curio_test_counts.tsv.gz` file. If you haven't yet installed the [Seurat](https://satijalab.org/seurat/) library, do so with `install.packages("Seurat")`

 - Read in the expression matrix with genes = rows and beads = columns and put in in the variable `counts`, then filter to beads with at least 10 UMIs and use it to create a Seurat object.
```
counts <- read.delim("curio_test_counts.tsv.gz", row.names=1, header=2)
counts <- counts[,colSums(counts)>=10]
seu <- CreateSeuratObject(counts=counts,project="curio_test_data")
```
 - Read the bead coordinates into the matrix `coords` with the bead barcodes as rownames, filter to the same beads as in the Seurat object, and add the bead coordinates as a "SPATIAL" dimensional reduction.
```
coords <- as.matrix(read.delim("A0010_039_BeadBarcodes.txt",row.names=1,header=FALSE))
coords <- as.matrix(coords[Cells(seu),]
colnames(coords) <- c("SPATIAL_1","SPATIAL_2")
``` 
 - Stardard processing: normalization, PCA, clustering, and finding spatially variable features.
```
seu <- SCTransform(seu, assay="RNA", ncells=1000, verbose=TRUE, conserve.memory=TRUE)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims=1:30)
seu <- FindNeighbors(seu, dims=1:30)
seu <- FindClusters(seu, resolution=0.2)
seu <- FindSpatiallyVariableFeatures(seu, assay="SCT", slot="scale.data", features=VariableFeatures(seu)[1:200], selection.method="moransi", x.cuts=100, y.cuts=100, verbose=TRUE, nfeatures=200)
```
 - Now you can use `saveRDS(seu,"curio_test_seurat_object.rds")` to save your new Seurat object, and you're ready to proceed with analysis! Repeat this process with the `nonSyrah` version of the data if desired. **NOTE** We have used the same parameters as the Curio-seeker pipeline defaults for this dataset. It's _highly_ likely that you'll want to teak these when processing any other data.

   
### Scanpy (Python)

(!!!to do)

## References 
<details>
* https://zenodo.org/records/10655615
* all main README references
  </details>
