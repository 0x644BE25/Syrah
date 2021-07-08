######################################################
# GENERATE CORRECTED ALIGNED SAM
#
# GOALS: 
# 1. Generate hash table of acceptable linker matches. 
# 2. Build hash table of barcode mappings from
#    previously generated feather file.
# 3. Use linker match alignments to find linker
#    position on read 1 (within given mapping 
#    parameters), pull barcode and UMI from read 1 
#    based on linker alignment start/stop. Look up 
#    barcode in mapping table, if barcode maps, write 
#    corrected read 2 to output SAM.
######################################################


# ================= IMPORTS ==========================

library(readr)
library(reader)
library(parallel)
library(doParallel)
library(data.table)

# ================= PARAMS ===========================

# computational parameters ===========================
nCores <- 100
batchSize <- 10^6
options(expressions=5e5) # needed for big old barcode matching hash
# The following parameter is the coefficient for how 
# large the read 1 chunk should be relative to the read 2 
# chunk. Increase this if you're getting "missing read"
# output, potentially decrease to improve total runtime.
# See Instructions for more info.
r1cushion <- 2

# linker matching parameters =========================
pattern <- 'JJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJNNNNNNNVVT'
linker <- stringi::stri_sub(pattern,9,26)
maxDist <-5
maxSubs <- 1
nts <- c('A','T','G','C')

# file paths =========================================
barcodesIn <- '<path to your barcode mapping feather>'
samR1 <- '<path to your headerless read 1 SAM>'
samR2 <- '<path to your headerless aligned read 2 SAM>'
samOut <- '<path to save the corrected aligned SAM>'

# SAM columns for unaligned read 1 ===================
colNamesR1 <- c('id',2:9,'seq',11:12)
colTypesR1 <- readr::cols(
  `id`=readr::col_character(),
  `2`=readr::col_double(),
  `3`=readr::col_character(),
  `4`=readr::col_double(),
  `5`=readr::col_double(),
  `6`=readr::col_character(),
  `7`=readr::col_character(),
  `8`=readr::col_double(),
  `9`=readr::col_double(),
  `seq`=readr::col_character(),
  `11`=readr::col_character(),
  `12`=readr::col_character()
)

# SAM columns for aligned read 2 =====================
colNamesR2 <- c('id',2:11,'xc',13:18,'xm',20:24)
colTypesR2 <- readr::cols(
  `id`=readr::col_character(),
  `2`=readr::col_double(),
  `3`=readr::col_character(),
  `4`=readr::col_double(),
  `5`=readr::col_double(),
  `6`=readr::col_character(),
  `7`=readr::col_character(),
  `8`=readr::col_double(),
  `9`=readr::col_double(),
  `10`=readr::col_character(),
  `11`=readr::col_character(),
  `xc`=readr::col_character(),
  `13`=readr::col_character(),
  `14`=readr::col_character(),
  `15`=readr::col_character(),
  `16`=readr::col_character(),
  `17`=readr::col_character(),
  `18`=readr::col_character(),
  `xm`=readr::col_character(),
  `20`=readr::col_character(),
  `21`=readr::col_character(),
  `22`=readr::col_character(),
  `23`=readr::col_character(),
  `24`=readr::col_character()
)

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
cat('starting\n')

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
cat('linker hash built!',nLinkHash,'entries\n')

# barcode match/defork hashmap =======================
temp <- data.frame(feather::read_feather(barcodesIn))
bcHash <- new.env(hash=TRUE)
for (i in 1:nrow(temp)) {
  bcHash[[temp[i,'from']]] <- temp[i,'to']
}
cat('bc hash built,',nrow(temp),'entries\n')
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

cl <- parallel::makeCluster(nCores,type="FORK")
cat('cluster made\n')

if (file.exists(samOut)) { cat('removing old tagged SAM'); unlink(samOut) }

if (TRUE) {
  totalReads <- 0
  leftovers <- NULL
  nR1 <- as.numeric(stringi::stri_split_fixed(system(paste('wc -l',samR1),intern=TRUE),' ')[[1]][1])
  nR2 <- as.numeric(stringi::stri_split_fixed(system(paste('wc -l',samR2),intern=TRUE),' ')[[1]][1])
  r1r2ratio <- nR1/nR2
  nBatches <- ceiling(nR2/batchSize)
  for (i in 1:nBatches) {
    
    # GET READ 1
    r1skip <- max(((i-(1+r1cushion/2))*batchSize*r1r2ratio),0)
    r1nReads <- r1cushion*batchSize*r1r2ratio
    r1 <- data.table::fread(file=samR1,sep='\t',col.names=colNamesR1,nrows=r1nReads,header=FALSE,skip=r1skip,data.table=FALSE)
    rownames(r1) <- r1$id
    
    # GET READ 2
    r2skip <- (i-1)*batchSize
    r2 <- data.table::fread(file=samR2,sep='\t',col.names=colNamesR2,nrows=batchSize,header=FALSE,skip=r2skip,data.table=FALSE)
    r2 <- rbind(leftovers,r2)
    
    # PULL OUT R2 NOT IN R1
    ind <- which(r2$id %in% r1$id)
    leftovers <- r2[-ind,]
    r2 <- r2[ind,]
    
    # CORRECT BC AND UMI
    xcxm <- parSapplyLB(cl,r1[r2$id,'seq'],getXC_XM2)
    r2$xc <- xcxm['XC',]
    r2$xm <- xcxm['XM',]
    r2 <- r2[nchar(r2$xc)>5,] # comment out to keep reads that don't map to puck
    
    readr::write_tsv(x=r2,file=samOut,append=TRUE,col_names=FALSE)
    totalReads <- totalReads+nrow(r2)
    if (nrow(leftovers)>0) { cat('batch',i,'written',totalReads,'total reads, missing',nrow(leftovers),'total reads\n') }
  }
}

parallel::stopCluster(cl)
cat('Done! Should be',totalReads,'total reads\n')

