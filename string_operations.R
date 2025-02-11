######################################################
# STRING OPERATIONS
#
# GOAL: Various string operations used by other R
# scripts to reduce dependencies on other libraries.
######################################################


del1 <- function(seq,pos=1:nchar(seq)) {
  seq2 <- strsplit(seq,'')[[1]]
  return(sapply(pos,\(i){
    paste0(seq2[-i],collapse='')
  }))
}

sub1 <- function(seq,pos=1:nchar(seq),nts=c('A','C','G','T')) {
  seq2 <- strsplit(seq,'')[[1]]
  return(c(sapply(pos,\(i){
    sapply(nts,\(nt){
      curr <- seq2; curr[i] <- nt; paste0(curr,collapse='')
    })
  })))
}

sub1notLast <- function(seq,pos=1:nchar(seq),nts=c('A','C','G','T')) {
  sub1(seq,pos=1:(nchar(seq)-1),nts=nts)
}

doMaxNdelMsub <- function(seq=linker,n,m,nts=c('A','C','G','T')) {
  # generate acceptable linker matches with n dels and m subs
  seqs <- res <- seq
  if (n>0) {
    for (x in 1:n) {
      seqs <- sapply(seqs,del1)
      res <- c(res,seqs)
    }
  }
  if (m>0) {
    seqs <- unique(c(res))
    for (x in 1:m) {
      seqs <- unlist(lapply(seqs,sub1,nts=nts))
      res <- c(res,seqs)
    }
  } 
  return(unique(unlist(res)))
}

do1hamm <- sub1

combn2hamm <- combn(1:14,2,list)
do2hamm <- function(seq) {
  seq2 <- strsplit(seq)
  sapply(combn2hamm)
}