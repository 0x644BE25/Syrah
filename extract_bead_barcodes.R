######################################################
# EXTRACT BEAD BARCODES
#
# GOAL: Use fuzzy matching to identify linker regions,
# and if a high-confidence match is found use the
# linker start/end positions to determine where the
# barcode parts 1 and 2 and UMI are. If there's a
# match for the barcode, append barcode and UMI
# sequences to teh read 2 sequence ID and write out.
# Discard all other reads.
######################################################

# ================= IMPORTS ==========================

# ================= PARAMS ===========================

# get manifest params
manifestFile <- commandArgs(trailingOnly=TRUE)[1]
true <- TRUE
false <- FALSE
source(manifestFile)

if (!endsWith(writeDir,'/')) { writeDir <- paste0(writeDir,'/') }
if (!endsWith(syrahDir,'/')) { syrahDir <- paste0(syrahDir,'/') }

# computational parameters
batchSize <- 10^5

# linker matching parameters
# NOTE: this uses the vs1 oligo pattern
# but the relevant parts don't differ between versions,
# only the UMI position
pattern <- 'JJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJTCNNNNNNNNT'
UMIpos <- list('yesTC'=c(36,43),'noTC'=c(34,40),'curio'=c(34,39))
linker <- substr(pattern,9,26)
maxSubs <- 1
nts <- c('A','C','G','T')

# ================= CHECK INPUT/OUTPUT FASTQS ===============

# OUTPUT
r2out <- paste0(writeDir,'intermediate_files/',batchName,'_r2_barcode_tagged.fastq')
if (file.exists(r2out)) { if (file.rename(r2out,paste0(r2out,'.OLD'))) { cat(paste0('moved old barcode extracted read2 FASTQ (Syrah) to ',r2out,'.OLD\n')) } }

# INPUT
if (endsWith(read1fastq,'gz')) {
  nR1 <- as.integer(strsplit(system(paste('zcat',read1fastq,'| wc -l'),intern=TRUE),' ')[[1]][1])/4
  nR2 <- as.integer(strsplit(system(paste('zcat',read2fastq,'| wc -l'),intern=TRUE),' ')[[1]][1])/4
} else {
  nR1 <- as.integer(strsplit(system(paste0('wc -l ',read1fastq),intern=TRUE),' ')[[1]][1])/4
  nR2 <- as.integer(strsplit(system(paste0('wc -l ',read2fastq),intern=TRUE),' ')[[1]][1])/4
}
if (nR1!=nR2) { cat('Error: different numbers of reads in read1 and read2 fastqs'); exit() }

r1x10K <- sapply(strsplit(readLines(read1fastq,n=40000)[seq(1,40000,by=4)],' '),\(x){ x[1] })
r2x10K <- sapply(strsplit(readLines(read2fastq,n=40000)[seq(1,40000,by=4)],' '),\(x){ x[1] })
if (!all(r1x10K==r2x10K)) { cat('CAUTION!!!!! FASTQs may not in same order. Please check and sort if needed.') }

# ================= PREP =============================

# GET READ1 VERSION
vs <- readLines(paste0(writeDir,'intermediate_files/',batchName,'_r1_version.txt'))[1]

# POSITION INFORMATION
maxLinkerDels <- 5
basepos <- c(bc1end=8,linkstart=9,linkend=26,bc2start=27,bc2end=32,umistart=UMIpos[[vs]][1],umiend=UMIpos[[vs]][2])
coeff <- c(0,0,1,1,1,1,1)
r1coords <- do.call(rbind,lapply(0:maxLinkerDels,\(dels){
  rbind(basepos-(coeff*dels),
        (basepos-1)-(coeff*dels))
}))

# ================= LINKER AND BARCODE MATCHES =======

source(paste0(syrahDir,'string_operations.R'))
linkerMatches <- doMaxNdelMsub(linker,maxLinkerDels,maxSubs,nts=c('A','C','G','T','N'))

barcodes <- unlist(lapply(readLines(paste0(writeDir,'intermediate_files/',batchName,'_barcode_whitelist.txt')),\(x){
  x <- strsplit(x,'\t')[[1]]
  froms <- strsplit(x[2],',')[[1]]
  res <- setNames(rep(x[1],length(froms)),froms)
  return(res)
}))

# ================= CORRECT BARCODE + UMI ============

conR1 <- file(read1fastq,'r')
conR2 <- file(read2fastq,'r')
writeR2 <- file(r2out,'a')

totalReads <- 0
cat('Starting barcode extraction...')

# I'M SORRY THIS NEXT PART IS SO UGLY
# BUT THE VECTORIZATION MAKES IT SO FAST
# CONSIDER AVERTING YOUR EYES
repeat({
  r1s <- readLines(conR1,n=batchSize*4)
  if (length(r1s)==0) { break() }
  r2s <- readLines(conR2,n=batchSize*4)
  id_ind <- seq(1,length(r1s),by=4)
  r1seqs <- r1s[id_ind+1]
  inds <- 1:length(r1seqs)
  res <- setNames(rep(0,length(inds)),inds)
  umis <- bcs <- setNames(rep(NA,length(inds)),inds)
  i <- 1
  while (length(inds)>0 & i<=nrow(r1coords)) {
    incurr <- inds[substr(r1seqs[inds],r1coords[i,'linkstart'],r1coords[i,'linkend']) %in% linkerMatches]
    res[incurr] <- i
    inds <- setdiff(inds,incurr)
    i <- i+1
  }
  good <- res>0
  bcs[good] <- barcodes[paste0(substr(r1seqs[good],1,r1coords[res[good],'bc1end']),
                               substr(r1seqs[good],r1coords[res[good],'bc2start'],r1coords[res[good],'bc2end']))]
  good <- !is.na(bcs)
  umis[good] <- substr(r1seqs[good],r1coords[res[good],'umistart'],r1coords[res[good],'umiend'])
  bc_umi <- paste0('_',bcs[good],'_',umis[good])
  good_inds <- id_ind[!is.na(bcs)]
  allr2inds <- sort(c(good_inds,good_inds+1,good_inds+2,good_inds+3))
  r2s <- r2s[allr2inds]
  id_ind <- seq(1,length(r2s),by=4)
  spl <- strsplit(r2s[id_ind],' ')
  r2s[id_ind] <- paste0(sapply(spl,\(x){ x[1] }),bc_umi,' ',sapply(spl,\(x){ x[2] }))
  writeLines(r2s,writeR2)
  totalReads <- totalReads+length(r1seqs)
  if (totalReads%%(10^6)==0) { cat('  ',format(as.POSIXlt(Sys.time())),' ',totalReads,'reads processed\n') }
  if (length(r1seqs)<batchSize) { break() }
})
close(conR1)
close(conR2)
close(writeR2)
cat('Syrah barcode extraction finished\n')

if (doNonSyrah) {
  
  # ================= PREP =============================
  
  # OUTPUT
  r2out <- paste0(writeDir,'intermediate_files/',batchName,'_r2_barcode_tagged_nonSyrah.fastq')
  if (file.exists(r2out)) { if (file.rename(r2out,paste0(r2out,'.OLD'))) { cat(paste0('moved old barcode extracted read2 FASTQ (non-Syrah) to ',r2out,'.OLD\n')) } }

  
  # GET READ1 VERSION
  vs <- readLines(paste0(writeDir,'intermediate_files/',batchName,'_r1_version.txt'))[1]
  
  # POSITION INFORMATION
  bc1end <- 8
  bc2start <-27
  bc2end <- 32
  umistart <- UMIpos[[vs]][1]
  umiend <- UMIpos[[vs]][2]
  
  # ================= BARCODE WHITELIST ================
  
  barcodes <- unlist(lapply(readLines(paste0(writeDir,'intermediate_files/',batchName,'_barcode_whitelist_nonSyrah.txt')),\(x){
    x <- strsplit(x,'\t')[[1]]
    froms <- strsplit(x[2],',')[[1]]
    res <- setNames(rep(x[1],length(froms)),froms)
    return(res)
  }))
  
  
  # ================= CORRECT BARCODE + UMI ============
  
  conR1 <- file(read1fastq,'r')
  conR2 <- file(read2fastq,'r')
  writeR2 <- file(r2out,'a')
  
  totalReads <- 0
  cat('Starting barcode extraction...')
  
  # NOT AS UGLY AS THE SYRAH VERSION
  # AT LEAST
  repeat({
    r1s <- readLines(conR1,n=batchSize*4)
    if (length(r1s)==0) { break() }
    r2s <- readLines(conR2,n=batchSize*4)
    id_ind <- seq(1,length(r1s),by=4)
    r1seqs <- r1s[id_ind+1]
    inds <- 1:length(r1seqs)
    res <- setNames(rep(0,length(inds)),inds)
    umis <- bcs <- setNames(rep(NA,length(inds)),inds)
    bcs <- barcodes[paste0(substr(r1seqs,1,bc1end),
                           substr(r1seqs,bc2start,bc2end))]
    good <- !is.na(bcs)
    umis[good] <- substr(r1seqs[good],umistart,umiend)
    bc_umi <- paste0('_',bcs[good],'_',umis[good])
    good_inds <- id_ind[!is.na(bcs)]
    allr2inds <- sort(c(good_inds,good_inds+1,good_inds+2,good_inds+3))
    r2s <- r2s[allr2inds]
    id_ind <- seq(1,length(r2s),by=4)
    spl <- strsplit(r2s[id_ind],' ')
    r2s[id_ind] <- paste0(sapply(spl,\(x){ x[1] }),bc_umi,' ',sapply(spl,\(x){ x[2] }))
    writeLines(r2s,writeR2)
    totalReads <- totalReads+length(r1seqs)
    if (totalReads%%(10^6)==0) { cat('  ',format(as.POSIXlt(Sys.time())),' ',totalReads,'reads processed\n') }
    if (length(r1seqs)<batchSize) { break() }
  })
  close(conR1)
  close(conR2)
  close(writeR2)
  cat('Non-Syrah barcode extraction finished\n')
}
