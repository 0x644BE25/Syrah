#' Make bead deduplication map
#'
#' Identify groups of virtually duplicted beads and re-direct each group to a single barcode, write out text file.
#' @param coords_file path to a tab-delimited file containing the barcodes and coordinates for the puck
#' @param write_dir (optional) directory to write to, defaults to current directory
#' @param max_slide_dist (optional) maximum slide distance to be considered a duplication, defaults to 10
#' @return none
#' @examples
#' make_bead_dedup_map(coords_file='coordinates.tsv')
#' make_bead_dedup_map(coords_file='coordinates.tsv',write_dir='~/mydir/',max_linker_dels=3)
#' @export
make_bead_dedup_map <- function(coords_file,write_dir='.',max_slide_dist=10) {
  
  if (!endsWith(write_dir,'/')) { write_dir <- paste0(write_dir,'/') }
  
  # PARAMS =======================
  nNeighbors <- 10
  maxSlideDist <- 10
  maxHammDist <- 1
  
  # hamming distance
  hamm <- function(a,b) {
    ab <- strsplit(c(a,b),'')
    sum(ab[[1]]!=ab[[2]])
  }
  
  # COORDINATES ============================
  puck <- read.delim(coords_file,header=FALSE)
  if (ncol(puck)==4) { puck <- puck[,c(1,3,4)] }
  colnames(puck) <- c('barcode','x','y')
  cat('Puck at',coords_file,'has',nrow(puck),'barcodes\n')
  
  # FILTER BY SLIDE + BC DISTANCE ==========
  nn <- dbscan::kNN(puck[,c('x','y')],k=nNeighbors)
  
  d <- data.frame(a=rep(puck$barcode,nNeighbors),b=puck$barcode[c(nn$id)],dist=c(nn$dist))
  d <- d[d$dist<=maxSlideDist,]
  
  #compare strings to remove rows where a,b are just switched order
  d <- d[d$a<d$b,]
  d$hamm <- apply(d,1,function(x){ hamm(as.character(x['a']),as.character(x['b'])) })
  d <- d[d$hamm==1,]
  
  # FIND DUPLICATES =========================
  dedup_groups <- unique(unlist(lapply(unique(d$a),function(x){
    bcs1 <- unique(c(x,d[d$a==x,'b']),d[d$b==x,'a'])
    curr <- d[(d$a %in% bcs1) | (d$b %in% bcs1),]
    bcs2 <- unique(c(curr$a,curr$b))
    curr <- d[(d$a %in% bcs2) | (d$b %in% bcs2),]
    bcs3 <- unique(c(curr$a,curr$b))
    if (length(bcs2)==length(bcs3)) {
      bcs3 <- sort(bcs3)
      return(paste(bcs3,collapse=','))
    }
  })))
  if (length(dedup_groups)==0) {
    dedup_groups <- character(0); dup_beads <- character(0)
  } else { 
    dup_beads <- unique(unlist(strsplit(dedup_groups,',')))
  }
  cat('...and',round((100*length(dup_beads)/nrow(puck)),1),'% are part of a virtual duplication\n')
  
  # WRITE DEDUP MAP ========================
  # add in un-duplicated barcodes so we don't lose them!
  dedup_groups <- c(dedup_groups,setdiff(puck$barcode,dup_beads))
  puck_name <- strsplit(coords_file,'/')[[1]]; puck_name <- puck_name[length(puck_name)]
  file_out <- paste0(write_dir,puck_name,'_dedup_map.txt')
  writeLines(dedup_groups,file_out)
  
  cat('Bead deduplication map is at',file_out,'\n')
}