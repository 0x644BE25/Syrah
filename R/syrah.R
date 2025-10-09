#' Syrah
#'
#' Wrapper for make_bead_dedup_map() + make_barcode_whitelist() + correct_barcodes()
#' @param coords_file path to a tab-delimited file containing the barcodes and coordinates for the puck
#' @param r1_fastq path to a read 1 fastq file -- may be comma delimited list of fastqs from same puck/library
#' @param write_dir (optional) directory to write to, defaults to current directory
#' @param n_cores (optional) integer indicating number of cores for parallelization, defaults to 1
#' @param max_slide_dist (optional) maximum slide distance to be considered a duplication, defaults to 10
#' @param max_linker_dels (optional) integer for max acceptable deletions in linker, defaults to 5
#' @param batch_size (optional) integer determining number of reads per batch, defaults to 10^5
#' @return none
#' @examples
#' syrah('coordinates.txt','read1.fastq')
#' syrah(coords_file='coordinates.txt',r1_fastq='read1.fastq,other_read1.fastq',write_dir='~/mydir/',n_cores=12)
#' @export
syrah <- function(coords_file,r1_fastq,write_dir='.',n_cores=1,max_slide_dist=1,max_linker_dels=5,batch_size=10^5) {
  
  if (!endsWith(write_dir,'/')) { write_dir <- paste0(write_dir,'/') }
  
  make_bead_dedup_map(coords_file=coords_file,write_dir=write_dir,max_slide_dist=max_slide_dist)
  puck_name <- strsplit(coords_file,'/')[[1]]; puck_name <- puck_name[length(puck_name)]
  dedup_map <- paste0(write_dir,puck_name,'_dedup_map.txt')
  
  make_barcode_whitelist(dedup_map=dedup_map,coords_file=coords_file,write_dir=write_dir,n_cores=n_cores)
  whitelist <- paste0(write_dir,puck_name,'_whitelist.txt')
  
  correct_barcodes(whitelist=whitelist,r1_fastq=r1_fastq,write_dir=write_dir,max_linker_dels=max_linker_dels,batch_size=batch_size)
}