#' Make barcode whitelist
#'
#' Use bead deduplication map to generate a whitelist of all acceptable + non-ambiguous barcode matches.
#' @param whitelist path to a bead dedupliction map, as created by make_barcode_whitelist()
#' @param r1_fastq path to a read 1 fastq file -- may be comma delimited list of fastqs from same puck/library
#' @param write_dir (optional) directory to write to, defaults to current directory
#' @param n_cores (optional) integer indicating number of cores for parallelization, defaults to 1
#' @param max_linker_dels (optional) integer for max acceptable deletions in linker, defaults to 5
#' @param batch_size (optional) integer determining number of reads per batch, defaults to 10^5
#' @return none
#' @examples
#' make_barcode_whtielist(dedup_map='coordinates.tsv_deduplication_map.txt')
#' make_barcode_whtielist(dedup_map='coordinates.tsv_deduplication_map.txt',write_dir='~/mydir/',n_cores=12)
#' @export
correct_barcodes <- function(whitelist,r1_fastq,write_dir='.',max_linker_dels=5,batch_size=10^5) {
  
  if (!endsWith(write_dir,'/')) { write_dir <- paste0(write_dir,'/') }
  
  r1_fastq <- unlist(strsplit(r1_fastq,',')[[1]])
  
 # PARAMS ==================================
  nts <- c('A','C','G','T')
  
  # computational parameters
  batchSize <- 10^5
  
  # linker matching parameters
  linker <- 'TCTTCAGCGTTCCCGAGA'
  maxSubs <- 1
  maxLinkerDels <- max_linker_dels

  # barcode matching
  bcPart1length <- 8
  bcPart2start <- bcPart1length+1
  bcLength <- 14
  r1bc2end <- bcLength+nchar(linker)
  
  # STRING METHODS =========================
  del1 <- function(seq,pos=1:nchar(seq)) {
    seq2 <- strsplit(seq,'')[[1]]
    return(sapply(pos,function(i){
      paste0(seq2[-i],collapse='')
    }))
  }
  
  sub1 <- function(seq,pos=1:nchar(seq),nts=c('A','C','G','T')) {
    seq2 <- strsplit(seq,'')[[1]]
    return(c(sapply(pos,function(i){
      sapply(nts,function(nt){
        curr <- seq2; curr[i] <- nt; paste0(curr,collapse='')
      })
    })))
  }
  
  sub1notLast <- function(seq,pos=1:nchar(seq),nts=c('A','C','G','T')) {
    sub1(seq,pos=1:(nchar(seq)-1),nts=nts)
  }
  
  doMaxNdelMsub <- function(seq=linker,n,m,nts=c('A','C','G','T')) {
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
  
  # LINKER MATCHES =========================
  
  linkerMatches <- doMaxNdelMsub(linker,maxLinkerDels,maxSubs,nts=c('A','C','G','T','N'))
  
  # BARCODE MATCHES ========================
  barcodes <- unlist(lapply(readLines(whitelist),function(x){
    x <- strsplit(x,'\t')[[1]]
    froms <- strsplit(x[2],',')[[1]]
    res <- setNames(rep(x[1],length(froms)),froms)
    return(res)
  }))
  
  # POSITIONS ==============================
  basepos <- c(bc1end=bcPart1length,
               linkstart=(bcPart1length+1),
               linkend=(bcPart1length+nchar(linker)),
               bc2start=(bcPart1length+nchar(linker)+1),
               bc2end=r1bc2end)
  coeff <- c(0,0,1,1,1)
  r1coords <- do.call(rbind,lapply(0:maxLinkerDels,function(dels){
    rbind(basepos-(coeff*dels),
          (basepos-1)-(coeff*dels))
  }))
  
  # UN-MATCHABLE SEQUENCE ==================
  r1_length <- nchar(readLines(r1_fastq[1],n=2)[2])
  unmatch <- paste0(rep('N',r1_length),collapse='')
  
  for (r1_fq in r1_fastq) {
  # FILEPATHS ==============================
  out_file <- strsplit(r1_fq,'/')[[1]]
  out_file <- out_file[length(out_file)]
  out_file <- paste0(write_dir,out_file,'.r1syrah')
  if (file.exists(out_file)) { file.rename(out_file,paste0(out_file,'.OLD')) }
  
  # barcode correction
  conR1 <- file(r1_fq,'r')
  writeR1 <- file(out_file,'a')
  totalReads <- 0
  cat(paste0('Starting barcode correction for ',r1_fq,' ...\n'))
  
  # I'M SORRY THIS NEXT PART IS SO UGLY
  # BUT THE VECTORIZATION MAKES IT SO FAST
  # CONSIDER AVERTING YOUR EYES
  repeat({
    r1s <- readLines(conR1,n=batchSize*4)
    if (length(r1s)==0) { break() }
    id_ind <- seq(1,length(r1s),by=4)
    r1seqs <- r1s[id_ind+1]
    inds <- 1:length(r1seqs)
    res <- setNames(rep(0,length(inds)),inds)
    bcs <- setNames(rep(NA,length(inds)),inds)
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
    new_r1_seqs <- paste0(substr(bcs[good],1,8),linker,substr(bcs[good],9,14),substr(r1seqs[good],basepos['bc2end']+1,r1_length))
    good_inds <- id_ind[!is.na(bcs)]+1
    bad_inds <- id_ind[is.na(bcs)]+1
    r1s[good_inds] <- new_r1_seqs
    r1s[bad_inds] <- unmatch
    writeLines(r1s,writeR1)
    
    totalReads <- totalReads+length(r1seqs)
    if (totalReads%%(10^6)==0) { cat('  ',format(as.POSIXlt(Sys.time())),' ',totalReads,'reads processed\n') }
    if (length(r1seqs)<batchSize) { break() }
  })
  closeAllConnections()
  cat('\nCorrected read 1 fastq written to',out_file,'\n')
  }
}
