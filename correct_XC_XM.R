######################################################
# CORRECT XC XM
#
# GOAL: Use read1 seqID + sequence data to align
# align linker sequnce, find BC position, map and
# defork BC, find UMI position, and replace both the
# BC (XC) and UMI (XM) in the aligned read2 data.
######################################################

# ================= IMPORTS ==========================

#library(stringi)
library(parallel)
library(data.table)
library(feather)

# ================= PARAMS ===========================

# file paths 
args <- commandArgs(trailingOnly=TRUE)
writeDir <- args[1]
batchName <- args[2]
vs <- args[3]
maxDist <- as.integer(args[4])
nCores <- as.integer(args[5])

samR1 <- paste0(writeDir,'intermediate_files/',batchName,'_r1_qname_filtered_qname_seq_only.txt')
samR2 <- paste0(writeDir,'intermediate_files/',batchName,'_r2_reads_only.sam')
bcMap <- paste0(writeDir,batchName,'_barcode_mapping.fa')
samOut <- paste0(writeDir,'intermediate_files/',batchName,'_XC_XM_corrected.txt')

# computational parameters
batchSize <- 10^6
options(expressions=5e5) # needed for big old barcode matching hash

# linker matching parameters
pattern <- 'JJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJTCNNNNNNNNT'
UMIpos <- list('vs1'=c(36,43),'vs2'=c(33,39))[[vs]]
linker <- stringi::stri_sub(pattern,9,26)
maxSubs <- 1
maxDistToCalc <- 8
nts <- c('A','T','G','C')

cat('read 1 in:',samR1,
    '\nread 2 in:',samR2,
    '\nUMI range:',UMIpos,
    '\nmax linker distance to calculate:',maxDistToCalc,
    '\nmax CPUs to use:',nCores,'\n')

# XC and XM replacement regex patterns
regexXC <- "XC\\:Z\\:[ATGC]{15}"
regexXM <- "XM\\:Z\\:[ATGC]{8}"

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
    for (i in combn(1:len,m,list)) {
      # (len-1) = don't allow subs in last position, expect del instead
      curr <- unname(sapply(combn(rep(nts,m),m,list),function(x){ return(stringi::stri_sub_replace_all(seq,from=i,to=i,value=x)) }))
      curr <- combn(positions[-i],n,function(x) {
        return(stringi::stri_sub_replace_all(curr,from=list(x),length=1,value=''))
      })
      res <- c(res,c(curr))
    }
    if (m>1) {
      rm <- NULL
      for (i in combn(1:len,(m-1),list)) {
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
  #return(unique(res))
  return(res)
}

doDdels <- function(seq=linker,d) {
  # generate acceptable linker matches with d deletions (including d-1 + 1 sub)
  #return(unique(c((unlist(sapply(0:min(maxDist-d,maxSubs),function(x){ doNdelMsub(seq,n=d,m=x) }))))))
  return(c((unlist(sapply(0:min(maxDistToCalc-d,maxSubs),function(x){ doNdelMsub(seq,n=d,m=x) })))))
}


# ================= BUILD HASHMAPS ===================
cat('starting\n')

# linker match hashmap
nLinkHash <- 0
linkHash <- new.env(hash=TRUE)
linkHash[[linker]] <- TRUE
for (i in 0:maxDistToCalc) {
  lseqs <- doDdels(seq=linker,d=i)
  #cat(i,length(lseqs),'\n')
  for (lseq in lseqs) {
    linkHash[[lseq]] <- TRUE
  }
  nLinkHash <- nLinkHash + length(lseqs)
}
cat('linker hash built!',nLinkHash,'entries\n')
  
# barcode match/defork hashmap
temp <- data.frame(feather::read_feather(bcMap))
bcHash <- new.env(hash=TRUE)
for (i in 1:nrow(temp)) {
  bcHash[[temp[i,'from']]] <- temp[i,'to']
}
cat('bc hash built,',nrow(temp),'entries\n')
rm(temp)

# ================= METHOD: GET_XC_XM ================

# This version includes barcode mapping + deforking step
getXC_XM_XX <- function(seq) {
  # 0 dels start 9
  if (!is.null(linkHash[[stringi::stri_sub(seq,9,26)]])) {
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,27,32))]]
    if (is.null(bc)) { fate <- '\tXX:Z:0d9s_noBC' } else { fate <- '\tXX:Z:0d9s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,33,39),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,25)]])) {
    # 0 dels start 8
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,26,31))]]
      if (is.null(bc)) { fate <- '\tXX:Z:0d8s_noBC' } else { fate <- '\tXX:Z:0d8s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,32,38),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,25)]])) {
    # 1 del start 9
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,26,31))]]
    if (is.null(bc)) { fate <- '\tXX:Z:1d9s_noBC' } else { fate <- '\tXX:Z:1d9s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,32,38),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,24)]])) {
    # 1 del start 8
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,25,30))]]
    if (is.null(bc)) { fate <- '\tXX:Z:1d8s_noBC' } else { fate <- '\tXX:Z:1d8s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,31,37),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,24)]])) {
    # 2 del start 9
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,25,30))]]
    if (is.null(bc)) { fate <- '\tXX:Z:2d9s_noBC' } else { fate <- '\tXX:Z:2d9s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,31,37),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,23)]])) {
    # 2 del start 8
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,24,29))]]
    if (is.null(bc)) { fate <- '\tXX:Z:2d8s_noBC' } else { fate <- '\tXX:Z:2d8s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,30,36),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,23)]])) {
    # 3 del start 9
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,24,29))]]
    if (is.null(bc)) { fate <- '\tXX:Z:3d9s_noBC' } else { fate <- '\tXX:Z:3d9s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,30,36),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,22)]])) {
    # 3 del start 8
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,23,28))]]
    if (is.null(bc)) { fate <- '\tXX:Z:3d8s_noBC' } else { fate <- '\tXX:Z:3d8s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,29,35),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,22)]])) {
    # 4 del start 9
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,23,28))]]
    if (is.null(bc)) { fate <- '\tXX:Z:4d9s_noBC' } else { fate <- '\tXX:Z:4d9s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,29,35),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,21)]])) {
    # 4 del start 8
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,22,27))]]
    if (is.null(bc)) { fate <- '\tXX:Z:4d8s_noBC' } else { fate <- '\tXX:Z:4d8s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,28,34),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,21)]])) {
    # 5 del start 9
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,22,27))]]
    if (is.null(bc)) { fate <- '\tXX:Z:5d9s_noBC' } else { fate <- '\tXX:Z:5d9s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,28,34),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,20)]])) {
    # 5 del start 8
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,21,26))]]
    if (is.null(bc)) { fate <- '\tXX:Z:5d8s_noBC' } else { fate <- '\tXX:Z:5d8s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,27,33),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,20)]])) {
    # 6 del start 9
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,21,26))]]
    if (is.null(bc)) { fate <- '\tXX:Z:6d9s_noBC' } else { fate <- '\tXX:Z:6d9s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,27,33),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,19)]])) {
    # 6 del start 8
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,20,25))]]
    if (is.null(bc)) { fate <- '\tXX:Z:6d8s_noBC' } else { fate <- '\tXX:Z:6d8s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,26,32),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,19)]])) {
    # 7 del start 9
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,20,25))]]
    if (is.null(bc)) { fate <- '\tXX:Z:7d9s_noBC' } else { fate <- '\tXX:Z:7d9s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,26,32),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,18)]])) {
    # 7 del start 8
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,19,24))]]
    if (is.null(bc)) { fate <- '\tXX:Z:7d8s_noBC' } else { fate <- '\tXX:Z:7d8s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,25,31),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,9,18)]])) {
    # 8 del start 9
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,8),stringi::stri_sub(seq,19,24))]]
    if (is.null(bc)) { fate <- '\tXX:Z:8d9s_noBC' } else { fate <- '\tXX:Z:8d9s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,25,31),fate)))
  } else if (!is.null(linkHash[[stringi::stri_sub(seq,8,17)]])) {
    # 8 del start 8
    bc <- bcHash[[paste0(stringi::stri_sub(seq,1,7),stringi::stri_sub(seq,18,23))]]
    if (is.null(bc)) { fate <- '\tXX:Z:8d8s_noBC' } else { fate <- '\tXX:Z:8d8s' }
    return(c(XC=paste0('XC:Z:',bc),
             XM=paste0('XM:Z:',stringi::stri_sub(seq,24,30),fate)))
  } else { 
    return(c(XC='XC:Z:',XM='XM:Z:\tXX:Z:0d0s'))
  }
}

# ================= MAIN =============================
# open file connections

conR1 <- file(samR1)
conR2 <- file(samR2)
open(conR1)
open(conR2)

cl <- parallel::makeCluster(nCores,type="FORK")
cat('cluster made\n')

if (file.exists(samOut)) { if (file.rename(samOut,paste0(samOut,'.OLD'))) { cat(paste0('moved old tagged SAM to ',samOut,'.OLD\n')) } }

if (TRUE) {
  nR1 <- as.numeric(stringi::stri_split_fixed(system(paste('wc -l',samR1),intern=TRUE),' ')[[1]][1])
  nR2 <- as.numeric(stringi::stri_split_fixed(system(paste('wc -l',samR2),intern=TRUE),' ')[[1]][1])
  if (nR1!=nR2) { stop('Number of lines in r1 and r2 files not equal. Ensure the same sequence IDs appear in both and that neither contain headers.')}
  nBatches <- ceiling(nR2/batchSize)
  cat(nBatches,'batches!\n')
  for (i in 1:nBatches) {
    hm <- ((i%%25)==1)
    if (hm) { cat('batch',i,'\n') }
    
    # GET READ 1
    r1 <- data.frame(do.call(rbind,stringi::stri_split_fixed(readLines(conR1,n=batchSize),'\t',n=2)),row.names=1)
    
    # GET READ 2
    r2 <- data.frame(do.call(rbind,stringi::stri_split_fixed(readLines(conR2,n=batchSize),'\t',n=2)))
    
    # REPLACE XC + XM
    xc_xm <- parallel::parSapplyLB(cl,r1[r2$X1,],getXC_XM_XX)
    corrReads <- list(paste(r2[,1],stringi::stri_replace_first_regex(stringi::stri_replace_first_regex(r2[,2],regexXC,xc_xm['XC',]),regexXM,xc_xm['XM',]),sep='\t'))
    data.table::fwrite(corrReads,file=samOut,append=TRUE,quote=FALSE)
  }
}

parallel::stopCluster(cl)
cat('ALL DONE!\n')
