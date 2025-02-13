######################################################
# GENERATE BEAD BARCODE WHITELIST
#
# GOAL: Make a whitelist of 1 deletion OR
# substiution bead barcode matches that correct to the
# "root" for forked beads. 

# If needed, make a max 1 hamming distance 
# whitelist for the "non-Syrah" pipeline.
######################################################

# ================= IMPORTS ==========================

library('parallel')

# ================= PARAMS ===========================

bcPart1length <- 8
bcPart2start <- bcPart1length+1
bcLength <- 14
nts <- c('A','C','G','T')

manifestFile <- commandArgs(trailingOnly=TRUE)[1]
true <- TRUE; false <- FALSE
source(manifestFile)

if (!endsWith(writeDir,'/')) { writeDir <- paste0(writeDir,'/') }

# ================= METHODS ==========================

# FOR SYRAH MATCHES
do1delOR1sub <- function(seq) {
  seq2 <- strsplit(seq,'')[[1]]
  
  del1 <- sapply(1:bcPart1length,function(x) { 
    paste0(seq2[-x],collapse='') })
  
  del2 <- c(sapply(sapply(bcPart2start:bcLength,function(x){
    paste0(seq2[-x],collapse='')}),\(y){ paste0(y,nts) }))
  
  sub <- c(sapply(1:bcLength,function(x) {
    sapply(nts,\(nt){ curr <- seq2; curr[x] <- nt; paste0(curr,collapse='') }) }))
  
  delOrSub <- unique(c(del1,del2,sub))
  delOrSub <- delOrSub[delOrSub!=seq]
  return(list(seq=seq,delOrSub=delOrSub))
}

# FOR NON-SYRAH MATCHES
do1sub <- function(seq) {
  seq2 <- strsplit(seq,'')[[1]]
  sub <- c(sapply(1:bcLength,function(x) {
    sapply(nts,\(nt){ curr <- seq2; curr[x] <- nt; paste0(curr,collapse='') }) }))
  
  return(sub)
}


# ================= INIT DATA ========================

puck <- read.delim(puckFile,header=FALSE)
colnames(puck) <- c('barcode','x','y')
dedup <- readLines(paste0(writeDir,'intermediate_files/',batchName,'_deduplication_map.txt'))


# ================= SYRAH VERSION ====================

# GENERATE MATCHES 
cl <- makeCluster(nCores,type='FORK')
matchBCs <- parLapply(cl,dedup,function(x){
  bcs <- strsplit(x,',')[[1]]
  d1 <- unique(unlist(lapply(bcs,do1delOR1sub)))
  d1 <- setdiff(d1,puck$barcode)
  return(data.frame(dedup=x,match=d1))
})
stopCluster(cl)
matchBCs <- do.call(rbind,matchBCs)

# DE-DUPLICATE
dups <- matchBCs$match[duplicated(matchBCs$match)]
matchBCs <- matchBCs[!matchBCs$match %in% dups,]

# COMBINE INTO CORRECT + dedup ROWS
cl <- makeCluster(nCores,type='FORK')
rows <- parLapply(cl,dedup,function(x){
  froms <- matchBCs[matchBCs$dedup==x,'match']
  root <- strsplit(x,',')[[1]][1]
  return(paste0(root,'\t',x,',',paste(froms,collapse=',')))
})
stopCluster(cl)

# WRITE FILE
fileName <- paste0(writeDir,'intermediate_files/',batchName,'_barcode_whitelist.txt')
writeLines(unlist(rows),fileName)
cat('Syrah barcode correction + deduping whitelist written to',fileName,'\n')


# ================= NON-SYRAH VERSION ================

if (doNonSyrah) {
  
  cl <- makeCluster(nCores,type='FORK')
  matchBCs <- parLapply(cl,puck$barcode,\(x){
    return(data.frame(puck=x,match=do1sub(x)))
  })
  stopCluster(cl)
  matchBCs <- do.call(rbind,matchBCs)
  matchBCs <- matchBCs[!matchBCs$match %in% puck$barcodes,]
  
  # DEDUPLICATE
  dups <- matchBCs$match[duplicated(matchBCs$match)]
  matchBCs <- matchBCs[!matchBCs$match %in% dups,]
  
  # COMBINE INTO WHITELIST ROWS
  cl <- makeCluster(nCores,type='FORK')
  rows <- parLapply(cl,puck$barcode,\(x){
    matches <- matchBCs[matchBCs$puck==x,'match']
    paste0(x,'\t',x,',',paste0(matches,collapse=','))
  })
  stopCluster(cl)
  
  # WRITE FILE
  fileName <- paste0(writeDir,'intermediate_files/',batchName,'_barcode_whitelist_nonSyrah.txt')
  writeLines(unlist(rows),fileName)
  cat('Non-Syrah barcode correction whitelist written to',fileName,'\n')
}
