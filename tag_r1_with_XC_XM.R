#!/usr/bin/env Rscript --vanilla

######################################################
# GENERATE CORRECTED ALIGNED SAM FROM BAMS
#
# GOAL: Use filtered, queryname-matched read 1 and 2
# BAMs to perform barcode and UMI correction and write
# out a corrected, aligned SAM file.
######################################################

# ================= IMPORTS ==========================


#! remove
source('~/rutils.R')

if (FALSE) {
  library(stringi)
  library(data.frame)
  library(feather)
  library(Rsamtools)
  library(parallel)
  library(readr)
}

# ================= PARAMS ===========================

args <- commandArgs(trailingOnly=TRUE)

# computational parameters ===========================

nCores <- as.integer(args[1])
batchSize <- 10^6
options(expressions=5e5) # needed for big old barcode matching hash


# file paths =========================================

writeDir <- args[2]
puckName <- args[3]
dataName <- args[4]

barcodesIn <- paste0(writeDir,'puck_',puckName,'_map.fa')
r1File <- paste0(writeDir,dataName,'_unaligned_r1.bam')
outFile <- paste0(writeDir,dataName,'_unaligned_r1_tagged.sam')


# linker matching parameters =========================

pattern <- 'JJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJTCNNNNNNNNT'
linker <- stringi::stri_sub(pattern,9,26)
maxDist <-5
maxSubs <- 1
nts <- c('A','T','G','C')

cat(nCores,barcodesIn,r1File,outFile,'\n')


# ================= METHODS ==========================

doNdelMsub <- function(seq=linker,n,m) {
  # generate acceptable linker matches with n dels and m subs
  len <- nchar(seq)
  positions <- 1:len
  if (m==0) {
    res <- c(combn(positions,n,function(x) {
      return(stringi::stri_sub_replace_all(seq,from=x,length=1,value=''))
    }))
  } else {
    res <- NULL
    for (i in combn(1:(len-1),m,list)) {
      # (len-1) = don't allow subs in last position, expect del instead
      curr <- unname(sapply(combn(rep(nts,m),m,list),function(x){ return(stringi::stri_sub_replace_all(seq,from=i,to=i,value=x)) }))
      curr <- combn(positions[-i],n,function(x) {
        return(stringi::stri_sub_replace_all(curr,from=list(x),length=1,value=''))
      })
      res <- c(res,c(curr))
    }
    if (m>1) {
      rm <- NULL
      for (i in combn(1:(len-1),(m-1),list)) {
        i <- c(i,len)
        curr <- unname(sapply(combn(rep(nts,m),m,list),function(x){ return(stringi::stri_sub_replace_all(seq,from=i,to=i,value=x)) }))
        curr <- combn(positions[-i],n,function(x) {
          return(stringi::stri_sub_replace_all(curr,from=x,length=1,value=''))
        })
        rm <- c(rm,c(curr))
      }
      res <- res[!res %in% unique(rm)]
    }
  }
  # res <- unique(res)
  return(unique(res))
}

doDdels <- function(seq=linker,d) {
  # generate acceptable linker matches with d deletions (including d-1 + 1 sub)
  return(unique(c((unlist(sapply(0:min(maxDist-d,maxSubs),function(x){ doNdelMsub(seq,n=d,m=x) }))))))
}

# ================= BUILD HASHMAPS ===================

options(expressions=5e5)

# linker match hashmap ================
nLinkHash <- 0
linkHash <- new.env(hash=TRUE)
linkHash[[linker]] <- TRUE
for (i in 0:maxDist) {
  lseqs <- doDdels(seq=linker,d=i)
  #cat(i,length(lseqs),'\n')
  for (lseq in lseqs) {
    linkHash[[lseq]] <- TRUE
  }
  nLinkHash <- nLinkHash + length(lseqs)
}

# barcode match/defork hashmap =======================
temp <- data.frame(feather::read_feather(barcodesIn))
bcHash <- new.env(hash=TRUE)
for (i in 1:nrow(temp)) {
  bcHash[[temp[i,'from']]] <- temp[i,'to']
}
rm(temp)

# ================= MORE METHODS =====================

# This version includes barcode mapping + deforking step
getXC_XM2 <- function(seq) {
  if (!is.null(linkHash[[stringi::stri_sub(seq,9,26)]])) {
    # 0 dels start 9
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,27,32))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,34,41))))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,25)]])) {
    # 0 dels start 8
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,26,31))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,33,40))))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,25)]])) {
    # 1 del start 9
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,26,31))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,33,40))))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,24)]])) {
    # 1 del start 8
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,25,30))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,32,39))))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,24)]])) {
    # 2 del start 9
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,25,30))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,32,39))))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,23)]])) {
    # 2 del start 8
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,24,29))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,31,38))))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,23)]])) {
    # 3 del start 9
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,24,29))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,31,38))))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,22)]])) {
    # 3 del start 8
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,23,28))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,30,37))))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,22)]])) {
    # 4 del start 9
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,23,28))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,30,37))))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,21)]])) {
    # 4 del start 8
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,22,27))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,29,36))))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,21)]])) {
    # 4 del start 9
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,22,27))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,29,36))))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,20)]])) {
    # 4 del start 8
    return(c(XC=paste0('XC:Z:',bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,21,26))]]),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,28,35))))
  } else { 
    return(c(XC='XC:Z:',XM='XM:Z:'))
  }
}

# ================= MAIN =============================

# setup BAM
r1Bam <- Rsamtools::BamFile(r1File)
Rsamtools::yieldSize(r1Bam) <- batchSize
open(r1Bam)

# clear old output file, if any
if (file.exists(outFile)) { cat('removing old XC + XM tagged SAM\n'); unlink(outFile) }

# setup parallelization
cl <- parallel::makeCluster(nCores,type="FORK")

nReads <- mappedReads <- 0
repeat {
  r1 <- cbind.data.frame(Rsamtools::scanBam(r1Bam,param=Rsamtools::ScanBamParam(what=c('qname','flag','rname','pos','mapq','cigar','mrnm','mpos','isize','seq','qual')))[[1]])
  r1$flag <- '4' # necessary to merge with aligned BAM
  r1$rname <- '*'
  r1$pos <- '0'
  r1$mapq <- '0'
  r1$cigar <- '*'
  r1$mrnm <- '*'
  r1$mpos <- '0'
  r1$isize <- '0'
  r1$XC <- ''
  r1$XM <- ''
  xcxm <- parallel::parSapplyLB(cl,r1$seq,getXC_XM2)
  r1$XC <- xcxm['XC',]
  r1$XM <- xcxm['XM',]
  good <- which(xcxm['XC',]>5)
  readr::write_tsv(x=r1[good,],file=outFile,append=TRUE,col_names=FALSE)
  readr::write_tsv(x=r1[-good,1:11],file=outFile,append=TRUE,col_names=FALSE)
  nReads <- nReads+nrow(r1)
  mappedReads <- mappedReads+length(good)
  if (nrow(r1)<batchSize) { break }
}
parallel::stopCluster(cl)
pct <- paste0('(',100*round((mappedReads/nReads),2),')%')
cat('## unaligned r1 with',nReads,'reads tagged with barcode and UMI.',mappedReads,'reads',pct,'mapped to barcodes. File at\n  ',outFile,'\n\n')
