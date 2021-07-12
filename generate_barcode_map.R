#!/usr/bin/env Rscript --vanilla

######################################################
# GENERATE BARCODE MAPPING
#
# GOAL: Find forked beads on the puck, create map so 
# that all branches of fork map to same barcode. Add
# non-forked barcodes mapped to self. Add in all 
# unique 1 hamming or 1 deletion matches. Save feather.
######################################################

# ================= DEPENDENCIES =====================

# I've referenced all these except for %dopar% 
# using :: so as long as  hey're all installed, no 
# need to load. This is really just for your reference
# as far as which script uses which packages.
library(foreach)

if (FALSE) {
  library(stringi)
  library(parallelDist)
  library(dplyr)
  library(stringdist)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(data.table)
  library(feather)
}

# ================= PARAMS ===========================

args <- commandArgs(trailingOnly=TRUE)
nCores <- as.integer(args[1])
puckFile <- args[2]
puckName <- args[3]
writeDir <- args[4]

chunkSize <- 10000

maxSlideDist <- 10  # max x,y distance on puck for forks
maxHammDist <- 1    # max Hamming distance for barcodes in forks

nts <- c('A','T','G','C')

# ================= METHODS ==========================

# get barcodes from coordinates in distance matrix
getpair <- function(coords){
  return(sort(puck$barcode[coords]))
}

# create distance = 1 variants
do1delOR1sub <- function(seq=linker,nts=c('A','T','G','C')) {
  positions <- 1:nchar(seq)
  del <- sapply(positions,function(x) { 
    paste0(stringi::stri_sub_replace(seq,from=x,length=1,value=''),nts)})
  sub <- sapply(positions,function(x) {
    stringi::stri_sub_replace(seq,from=x,length=1,value=nts)
  })
  res <- unique(c(del,sub))
  res <- res[res!=seq]
  return(res)
}

doBatch <- function(pairs,nts=c('A','T','G','C')) {
  res <- NULL
  for (i in 1:nrow(pairs)) {
    bcs <- do1delOR1sub(pairs[i,'from'],nts=nts)
    res <- rbind(res,data.frame(from=bcs,to=pairs[i,'to']))
  }
  good <- table(res$from)
  good <- names(good[good==1])
  good <- good[!good %in% uqpairs$from]
  return(res[res$from %in% good,])
}


# ================= FORKING DATA =====================

puck <- read.csv(puckFile,sep='\t',header=FALSE,col.names=c('barcode','x','y'))
rownames(puck) <- puck$barcode


# get bead-bead distances =======================

dists <- parallelDist::parDist(as.matrix(puck[,c('x','y')]),method="euclidean",threads=nCores)
mdists <- as.matrix(dists)
mdists[mdists==0] <- 999


# get forked pairs ==============================

paircoords <- which(mdists<maxSlideDist,arr.ind=TRUE)
pairs <- apply(paircoords,1,getpair)
uqpairs <- dplyr::distinct(data.frame(t(unname(pairs))))
colnames(uqpairs) <- c('to','from')
uqpairs$dist.ham <- stringdist::stringdist(uqpairs$from,uqpairs$to,method='hamming')
uqpairs <- uqpairs[uqpairs$dist.ham<=maxHammDist,]


# redirect secondary branches to root ===========

for (i in 2:nrow(uqpairs)) {
  while (uqpairs[i,'to'] %in% uqpairs[,'from']) {
    uqpairs[i,'to'] <- uqpairs$to[uqpairs$from==uqpairs[i,'to']][1]
  }
}
dups <- table(uqpairs$from)
dups <- names(dups[dups>1])
for (dup in dups) {
  uqpairs[uqpairs$from %in% dups,c('to','from')] <- uqpairs[uqpairs$from %in% dups,c('from','to')]
}


# ================= NON-FORKED BEADS =================

bcsNotForked <- puck$barcode[!puck$barcode %in% uqpairs$from]
uqpairs <- rbind(uqpairs,data.frame(to=bcsNotForked,from=bcsNotForked,dist.ham=0))

# ================= DISTANCE 1 MATCHES ===============

cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)

nPairs <- nrow(uqpairs)
nChunks <- ceiling(nPairs/chunkSize)

oneoff <- foreach(i=1:nChunks) %dopar% {
  start <- (i-1)*chunkSize+1
  end <- min(start+chunkSize-1,nPairs)
  curr <- uqpairs[start:end,]
  cat(start,end,nrow(curr),'\n')
  doBatch(uqpairs[start:end,])
}
parallel::stopCluster(cl)
oneoff <- data.table::rbindlist(oneoff)
good <- table(oneoff$from)
good <- names(good[good==1])
good <- good[!good %in% uqpairs$from]
oneoff <- oneoff[oneoff$from %in% good,]
uqpairs <- uqpairs[,c('from','to')]
uqpairs <- rbind(uqpairs,oneoff)

# ================= WRITE DATA =======================

feather::write_feather(uqpairs,paste0(writeDir,'puck_',puckName,'_map.fa'))
cat('\n## barcode mapping containing',nrow(uqpairs),'total entries saved at\n  ',paste0(writeDir,'puck_',puckName,'_map.fa'),'\n\n')
