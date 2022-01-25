#!/usr/bin/env Rscript --vanilla

######################################################
# CREATE BEAD BARCODE MAPPING
#
# GOAL: 
######################################################


# ================= IMPORTS ==========================

library(feather)
library(parallel)
library(doParallel)
library(foreach)
library(parallelDist)
#library(stringi)
library(stringdist)

# ================= PARAMS ===========================

# user
args <- commandArgs(trailingOnly=TRUE)
puckFile <- args[1]
writeDir <- args[2]
batchName <- args[3]
nCores <- as.integer(args[4])
batch <- args[2]
saveAs <- paste0(writeDir,batchName,'_barcode_mapping.fa')

# computational
chunkSize <- 100

# deforking
nBins <- 10
maxSlideDist <- 10
maxHammDist <- 1

# match generating
bcLength <- 14
bcPart1length <- 8
bcPart2start <- bcPart1length+1
nts <- c('A','T','G','C')

# ================= METHODS ==========================

# create distance = 1 variants
do1delOR1sub <- function(seq=linker,nts=c('A','T','G','C')) {
  # deletions in part 1, only 13 nt
  del1 <- sapply(1:bcPart1length,function(x) { 
    paste0(stringi::stri_sub_replace(seq,from=x,length=1,value=''))})
  # deletions in part 2, 14 nt
  del2 <- c(sapply((bcPart2start:bcLength),function(x) { 
    paste0(stringi::stri_sub_replace(seq,from=x,length=1,value=''),nts)}))
  # substitutions, 14nt
  sub <- c(sapply(1:bcLength,function(x) {
    stringi::stri_sub_replace(seq,from=x,length=1,value=nts)}))
  res <- unique(c(del1,del2,sub))
  res <- res[res!=seq]
  return(res)
}

doBatch <- function(pairs,nts=c('A','T','G','C')) {
  res <- NULL
  for (i in 1:nrow(pairs)) {
    bcs <- do1delOR1sub(pairs[i,'from'],nts=nts)
    res <- rbind(res,data.frame(from=bcs,to=pairs[i,'to']))
  }
  return(res)
}

# ================= DEFORKING ========================

puck <- read.csv(puckFile,sep='\t',header=FALSE,col.names=c('barcode','x','y'))

xbinSize <- (max(puck$x)-min(puck$x))/nBins
ybinSize <- (max(puck$y)-min(puck$y))/nBins

# find beads within distance threshold
# work with overlapping sub-areas of puck for speed
pairs <- NULL
for (i in 1:nBins) {
  xmin <- min(puck$x)+(i-1)*xbinSize
  xmax <- xmin+xbinSize+10
  for (j in 1:nBins) {
    ymin <- min(puck$y)+(j-1)*ybinSize
    ymax <- ymin+ybinSize+10
    
    curr <- puck[puck$x>=xmin & puck$x<xmax & puck$y>=ymin & puck$y<ymax,]
    #cat(i,' X:',xmin,xmax,'   Y:',ymin,ymax,'   ',nrow(curr),'\n')
    if (nrow(curr)>=2) {
      dists <- as.matrix(parallelDist::parDist(as.matrix(curr[,c('x','y')]),method="euclidean",threads=nCores))
      dists[!upper.tri(dists)] <- NA
      paircoords <- which(dists<maxSlideDist,arr.ind=TRUE)
      if (nrow(paircoords)>0) {
        pairNames <- apply(paircoords,1,function(coords){ sort(curr$barcode[coords]) })
        distHam <- stringdist::stringdist(pairNames[1,],pairNames[2,],method='hamming')
        pairs <- cbind(pairs,pairNames[,distHam==1])
      }
    }
  }
}
pairs <- dplyr::distinct(data.frame(t(pairs)))
colnames(pairs) <- c('from','to')

# deal with non-unique from barcodes
fromCounts <- table(pairs$from)
bad <- pairs[pairs$from %in% names(fromCounts[fromCounts>1]),]
pairs <- pairs[pairs$from %in% names(fromCounts[fromCounts==1]),]
for (b in unique(bad$from)) {
  curr <- bad[bad$from==b,][1,]
  pairs <- rbind(pairs,curr[1,])
  bads <- curr[2:nrow(curr),'to']
  pairs[pairs$to %in% bads,'to'] <- curr[1,'to']
}

# redirect secondary branches to root
if (nrow(pairs)>=2) {
  for (i in 1:nrow(pairs)) {
    while (pairs[i,'to'] %in% pairs[,'from']) {
      pairs[i,'to'] <- pairs$to[pairs$from==pairs[i,'to']][1]
    }
  }
}

# How did we do?
bcsInForks <- unique(c(pairs$from,pairs$to))

# ================= 1 DISTANCE MATCHES ===============

# add in non-forks
bcsNotForked <- puck$barcode[!puck$barcode %in% pairs$from]
pairs <- rbind(pairs,data.frame(to=bcsNotForked,from=bcsNotForked))

# build in 1 distance matches
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)

nPairs <- nrow(pairs)
nChunks <- ceiling(nPairs/chunkSize)

oneoff <- foreach(i=1:nChunks) %dopar% {
  start <- (i-1)*chunkSize+1
  end <- min(start+chunkSize-1,nPairs)
  curr <- pairs[start:end,]
  cat(start,end,nrow(curr),'\n')
  doBatch(pairs[start:end,])
}
parallel::stopCluster(cl)

oneoff <- data.table::rbindlist(oneoff)
good <- oneoff$from[!oneoff$from %in% pairs$from]
tab <- table(good)
good <- good[good %in% names(tab[tab==1])]
good <- good[!good %in% puck$barcode]
pairs <- rbind(pairs,oneoff[oneoff$from %in% good,])

# ================= SAVE DATA ========================

feather::write_feather(pairs,paste0(saveAs))
cat('map saved at',saveAs,'\n',
    paste0(round(100*length(bcsInForks)/nrow(puck),2),'% of puck barcodes are part of a fork!\n'),
    nrow(pairs),'total matches\n')
