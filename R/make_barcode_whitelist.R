#' Make barcode whitelist
#'
#' Use bead deduplication map to generate a whitelist of all acceptable + non-ambiguous barcode matches.
#' @param dedup_map path to a bead dedupliction map, as created by make_bead_dedup_map()
#' @param coords_file path to a tab-delimited file containing the barcodes and coordinates for the puck
#' @param write_dir (optional) directory to write to, defaults to current directory
#' @param n_cores (optional) integer indicating number of cores for parallelization, defaults to 1
#' @return none
#' @examples
#' make_barcode_whtielist(dedup_map='coordinates.tsv_deduplication_map.txt')
#' make_barcode_whtielist(dedup_map='coordinates.tsv_deduplication_map.txt',write_dir='~/mydir/',n_cores=12)
#' @export
make_barcode_whitelist <- function(dedup_map,coords_file,write_dir='.',n_cores=1) {
  
  if (!endsWith(write_dir,'/')) { write_dir <- paste0(write_dir,'/') }
  
  # PARAMS =================================
  nts <- c('A','C','G','T')
  maxHammDist <- 1
  
  bcPart1length <- 8
  bcPart2start <- bcPart1length+1
  bcLength <- 14
  
  # STRING METHODS =========================
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
  
  do1sub <- function(seq) {
    seq2 <- strsplit(seq,'')[[1]]
    sub <- c(sapply(1:bcLength,function(x) {
      sapply(nts,\(nt){ curr <- seq2; curr[x] <- nt; paste0(curr,collapse='') }) }))
    
    return(sub)
  }
  
  # PUCK INFO ==============================
  puck <- read.delim(coords_file,header=FALSE)
  if (ncol(puck)==4) { puck <- puck[,c(1,3,4)] }
  colnames(puck) <- c('barcode','x','y')
  
  # DEDUP MAP ==============================
  dedup <- readLines(dedup_map)
  # GENERATE MATCHES 
  cl <- parallel::makeCluster(n_cores,type='FORK')
  matchBCs <- parallel::parLapply(cl,dedup,function(x){
    bcs <- strsplit(x,',')[[1]]
    d1 <- unique(unlist(lapply(bcs,do1delOR1sub)))
    d1 <- setdiff(d1,puck$barcode)
    return(data.frame(dedup=x,match=d1))
  })
  parallel::stopCluster(cl)
  matchBCs <- do.call(rbind,matchBCs)
  
  # DE-DUPLICATE
  dups <- matchBCs$match[duplicated(matchBCs$match)]
  matchBCs <- matchBCs[!matchBCs$match %in% dups,]
  
  # COMBINE INTO CORRECT + dedup ROWS
  cl <- parallel::makeCluster(n_cores,type='FORK')
  rows <- parallel::parLapply(cl,dedup,function(x){
    froms <- matchBCs[matchBCs$dedup==x,'match']
    root <- strsplit(x,',')[[1]][1]
    return(paste0(root,'\t',x,',',paste(froms,collapse=',')))
  })
  parallel::stopCluster(cl)
  
  # WRITE FILE
  puck_name <- strsplit(dedup_map,'/')[[1]]
  puck_name <- gsub('dedup_map','whitelist',puck_name[length(puck_name)])
  file_out <- paste0(write_dir,puck_name)
  writeLines(unlist(rows),file_out)
  cat('Syrah barcode correction + deduping whitelist written to',file_out,'\n')
}