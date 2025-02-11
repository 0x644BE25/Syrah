######################################################
# CREATE DEFORK MAP
#
# GOAL: Use x,y distance and barcode hamming distance
# to determine forks, then route all barcodes in each
# fork to a single barcode. This will be used during
# the generate_whitelist step.
######################################################

# ================= IMPORTS ==========================

# ================= PARAMS ===========================

nNeighbors <- 10
maxSlideDist <- 10
maxHammDist <- 1

manifestFile <- commandArgs(trailingOnly=TRUE)[1]
true <- TRUE; false <- FALSE
source(manifestFile)

if (!endsWith(writeDir,'/')) { writeDir <- paste0(writeDir,'/') }

# ================= METHODS ==========================

hamm <- function(a,b) {
  ab <- strsplit(c(a,b),'')
  sum(ab[[1]]!=ab[[2]])
}

# ================= INIT DATA ========================

puck <- read.delim(puckFile,header=FALSE)
colnames(puck) <- c('barcode','x','y')
cat('Puck at',puckFile,'has',nrow(puck),'barcodes\n')

# ================= ANALYSIS =========================

nn <- dbscan::kNN(puck[,c('x','y')],k=nNeighbors)

d <- data.frame(a=rep(puck$barcode,nNeighbors),b=puck$barcode[c(nn$id)],dist=c(nn$dist))
d <- d[d$dist<=maxSlideDist,]

#compare strings to remove rows where a,b are just switched order
d <- d[d$a<d$b,]
d$hamm <- apply(d,1,function(x){ hamm(as.character(x['a']),as.character(x['b'])) })
d <- d[d$hamm==1,]

defork_groups <- unique(unlist(lapply(unique(d$a),function(x){
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

fork_beads <- unique(unlist(strsplit(defork_groups,',')))
cat('...and',round((100*length(fork_beads)/nrow(puck)),1),'% are part of a fork\n')

# ================= SAVE DATA ========================

# add in non-forked barcodes so we don't lose them!
defork_groups <- c(defork_groups,setdiff(puck$barcode,fork_beads))
fileName <- paste0(writeDir,batchName,'_defork_groups.txt')
writeLines(defork_groups,fileName)

cat('Deforking map is at',fileName,'\n')
