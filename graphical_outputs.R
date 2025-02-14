######################################################
# GRAPHICAL OUTPUTS
#
# GOAL: Plot estimated read 1 nucleotide frequency as
# well as min=1 and min=25 nUMI versions of nFeature
# and nUMI violin plots and nUMI slide embedding.
######################################################

# ================= IMPORTS ==========================

# NONE, Y'ALL. Not even ggplot.

# ================= PARAMS ===========================

manifestFile <- commandArgs(trailingOnly=TRUE)[1]
true <- TRUE
false <- FALSE
source(manifestFile)
if (!endsWith(writeDir,'/')) { writeDir <- paste0(writeDir,'/') }

minUMI <- 10

spaces <- list(yesTC=c(rep(0,8),.25,rep(0,17),.25,rep(0,5),.25,rep(0,2),.25,rep(0,7),.25,rep(0,99)),
               noTC=c(rep(0,8),.25,rep(0,17),.25,rep(0,5),.25,.25,rep(0,7),.25,rep(0,99)),
               curio=c(rep(0,8),.25,rep(0,17),.25,rep(0,5),.25,.25,rep(0,5),.25,rep(0,1),.25,rep(0,99)))

# ================= METHODS ==========================


blendColors <- function(A,B,ratio) {
  ratio <- min(1,max(ratio,0))
  if (!startsWith(A,'#')) { A <- gplots::col2hex(A) }
  if (!startsWith(B,'#')) { B <- gplots::col2hex(B) }
  # ratio = b/total
  all <- NULL
  for (col in A) {
    #print(col)
    res <- '#'
    for (i in 1:3) {
      a <- strtoi(paste0('0x',substr(col,(2*i),2*i+1)))
      b <- strtoi(paste0('0x',substr(B,(2*i),2*i+1)))
      ab <- min(255,ceiling(((1-ratio)*a)+(ratio*b)))
      ab <- as.character(as.hexmode(ab))
      if (nchar(ab)==1) { ab <- paste0('0',ab)}
      res <- c(res,ab)
    }
    res <- paste0(res,collapse='')
    all <- c(all,res)
    #print(res)
  }
  return(all)
}


blendNcolors <- function(colors,ratio) {
  if (ratio==0) { return(colors[1]) }
  if (ratio==1) { return(colors[length(colors)]) }
  n <- length(colors)
  colors <- c(colors[1],colors)
  breakpts <- (0:n)/n
  thisn <- which(breakpts>ratio)[1]
  ratio <- n*(ratio-breakpts[thisn-1])
  return(blendColors(colors[thisn-1],colors[thisn],ratio))
}

# ================= INIT DATA ========================

vs <- readLines(paste0(writeDir,'intermediate_files/',batchName,'_r1_version.txt'))
patt <- readLines(paste0(writeDir,'intermediate_files/',batchName,'_r1_pattern.txt'))
freqs <- read.csv(paste0(writeDir,'intermediate_files/',batchName,'_estimated_read_1_nucleotide_frequencies.csv'),row.names=1)
puck <- read.delim(puckFile,row.names=1,header=FALSE)
syrah <- read.delim(paste0(writeDir,batchName,"_counts.tsv.gz"),row.names=1)

if (doNonSyrah) { std <-  read.delim(paste0(writeDir,batchName,"_nonSyrah_counts.tsv.gz"),row.names=1) }

# ================= PREP DATA ========================

chars <- c('X'='','C'='B','N'='U')[strsplit(patt,'')[[1]]]
colnames(freqs) <- paste0(chars,'\n',1:ncol(freqs))

# ================= SYRAH VERSION ========================

df <- data.frame(row.names=colnames(syrah),
                 nGenes=colSums(syrah>0),
                 nUMI=colSums(syrah),
                 x=puck[colnames(syrah),1],
                 y=puck[colnames(syrah),2])
df$colNum <- (log10(df$nUMI)-min(log10(df$nUMI)))/(max(log10(df$nUMI))-min(log10(df$nUMI)))
df$col <- sapply(df$colNum,\(x){ blendNcolors(colors=c('#ffffcc','#ffcc00','#ff8800','#ff0000','#660022','#000000'),x) })

df.overMin <- df[df$nUMI>=minUMI,]
df.overMin$colNum <- (log10(df.overMin$nUMI)-min(log10(df.overMin$nUMI)))/(max(log10(df.overMin$nUMI))-min(log10(df.overMin$nUMI)))
df.overMin$col <- sapply(df.overMin$colNum,\(x){ blendNcolors(colors=c('#ccffff','#00ccff','#0088ff','#0000ff','#220066','#000000'),x) })

nGenes <- sum(rowSums(syrah)>0)
nGenes.overMin <- sum(rowSums(syrah[,colSums(syrah>=minUMI)])>0)

try({ shh <- dev.off() },silent=TRUE)

pdf(paste0(writeDir,batchName,'_Syrah_results_summary.pdf'),width=8.5,height=11)
nf <- layout( matrix(c(1,1,2,4,3,4,5,7,6,7), ncol=2, byrow=TRUE))

par(mar=c(4,2,2,0))
barplot(as.matrix(freqs[c('T','G','C','A','N'),]), 
        col=c('#ff585f','#ffcb05','#73ab1c','#496efd','#CCCCCC'),
        border=NA,
        space=spaces[[vs]][1:ncol(freqs)],
        legend.text=c('T','G','C','A','N'),
        yaxt='n',
        cex.names=.45,
        args.legend=list(y=1,x=-.1,box.col='transparent'),
        offset=0,
        main=paste(batchName,'estimated read 1 nucleotide frequency by position'),
        xlab="5' ----------- nt position ---------> 3'")

par(mar=c(5,4,4,2)+0.1)
# MIN 1 NUMI
den <- density(log10(df$nUMI))
plot(den,
     col='brown',lwd=2,
     xlab='log10 UMI count',ylab='',yaxt='n',
     main=paste0('nUMI distribution (>=1 UMIs, ',nrow(df),' beads)'),
     sub=paste0('median: ',median(df$nUMI),'   mean: ',round(mean(df$nUMI),2)))
polygon(den,col='#ff000044')
abline(v=median(log10(df$nUMI)))
abline(v=mean(log10(df$nUMI)),lty='dashed')

den <- density(log10(df$nUMI))
plot(den,
     col='brown',lwd=2,
     xlab='log10 unique feature count',ylab='',yaxt='n',
     main=paste0('nFeature distribution (>=1 UMIs, ',nGenes,' features)'),
     sub=paste0('median: ',median(df$nGene),'   mean: ',round(mean(df$nGene),2)))
polygon(den,col='#ff000044')
abline(v=median(log10(df$nGene)))
abline(v=mean(log10(df$nGene)),lty='dashed')

plot(df$x~df$y,col=df$col,
     pch=20,cex=.2,xlab='',ylab='',main='nUMI (>=1 UMIs)',asp=1)

# MIN 25 NUMI
den <- density(log10(df.overMin$nUMI))
plot(den,
     col='#0000ff',lwd=2,
     xlab='log10 UMI count',ylab='',yaxt='n',
     main=paste0('nUMI distribution (>=',minUMI,' UMIs, ',nrow(df.overMin),' beads)'),
     sub=paste0('median: ',median(df.overMin$nUMI),'   mean: ',round(mean(df.overMin$nUMI),2)))
polygon(den,col='#0088ff44')
abline(v=median(log10(df.overMin$nUMI)))
abline(v=mean(log10(df.overMin$nUMI)),lty='dashed')

den <- density(log10(df.overMin$nGene))
plot(den,
     col='blue',lwd=2,
     xlab='log10 feature count',ylab='',yaxt='n',
     main=paste0('nFeature distribution (>=',minUMI,' UMIs, ',nGenes.overMin,' unique features)'),
     sub=paste0('median: ',median(df.overMin$nGene),'   mean: ',round(mean(df.overMin$nGene),2)))
polygon(den,col='#0088ff44')
abline(v=median(log10(df.overMin$nGene)))
abline(v=mean(log10(df.overMin$nGene)),lty='dashed')

plot(df.overMin$x~df.overMin$y,col=df.overMin$col,
     pch=20,cex=.2,xlab='',ylab='',main=paste0('nUMI (>=',minUMI,' UMIs)'),asp=1)

shh <- dev.off()

# ================= NON-SYRAH VERSION ================


if (doNonSyrah) {
  df <- data.frame(row.names=colnames(std),
                   nGenes=colSums(std>0),
                   nUMI=colSums(std),
                   x=puck[colnames(std),1],
                   y=puck[colnames(std),2])
  df$colNum <- (log10(df$nUMI)-min(log10(df$nUMI)))/(max(log10(df$nUMI))-min(log10(df$nUMI)))
  df$col <- sapply(df$colNum,\(x){ blendNcolors(colors=c('#ffffcc','#ffcc00','#ff8800','#ff0000','#660022','#000000'),x) })
  
  df.overMin <- df[df$nUMI>=minUMI,]
  df.overMin$colNum <- (log10(df.overMin$nUMI)-min(log10(df.overMin$nUMI)))/(max(log10(df.overMin$nUMI))-min(log10(df.overMin$nUMI)))
  df.overMin$col <- sapply(df.overMin$colNum,\(x){ blendNcolors(colors=c('#ccffff','#00ccff','#0088ff','#0000ff','#220066','#000000'),x) })
  
  nGenes <- sum(rowSums(std)>0)
  nGenes.overMin <- sum(rowSums(std[,colSums(std>=minUMI)])>0)
  
  try({ shh <- dev.off() },silent=TRUE)
  
  pdf(paste0(writeDir,batchName,'_nonSyrah_results_summary.pdf'),width=8.5,height=11)
  nf <- layout( matrix(c(1,1,2,4,3,4,5,7,6,7), ncol=2, byrow=TRUE))
  
  par(mar=c(4,2,2,0))
  barplot(as.matrix(freqs[c('T','G','C','A','N'),]), 
          col=c('#ff585f','#ffcb05','#73ab1c','#496efd','#CCCCCC'),
          border=NA,
          space=spaces[[vs]][1:ncol(freqs)],
          legend.text=c('T','G','C','A','N'),
          yaxt='n',
          cex.names=.45,
          args.legend=list(y=1,x=-.1,box.col='transparent'),
          offset=0,
          main=paste(batchName,'estimated read 1 nucleotide frequency by position'),
          xlab="5' ----------- nt position ---------> 3'")
  
  # MIN 1 NUMI
  den <- density(log10(df$nUMI))
  plot(den,
       col='brown',lwd=2,
       xlab='log10 UMI count',ylab='',yaxt='n',
       main=paste0('nUMI distribution (>=1 UMIs, ',nrow(df),' beads)'),
       sub=paste0('median: ',median(df$nUMI),'   mean: ',round(mean(df$nUMI),2)))
  polygon(den,col='#ff000044')
  abline(v=median(log10(df$nUMI)))
  abline(v=mean(log10(df$nUMI)),lty='dashed')
  
  den <- density(log10(df$nUMI))
  plot(den,
       col='brown',lwd=2,
       xlab='log10 unique feature count',ylab='',yaxt='n',
       main=paste0('nFeature distribution (>=1 UMIs, ',nGenes,' features)'),
       sub=paste0('median: ',median(df$nGene),'   mean: ',round(mean(df$nGene),2)))
  polygon(den,col='#ff000044')
  abline(v=median(log10(df$nGene)))
  abline(v=mean(log10(df$nGene)),lty='dashed')
  
  
  plot(df$x~df$y,col=df$col,
       pch=20,cex=.2,xlab='',ylab='',main='nUMI (>=1 UMIs)',asp=1)
  
  # MIN 25 NUMI
  den <- density(log10(df.overMin$nUMI))
  plot(den,
       col='#0000ff',lwd=2,
       xlab='log10 UMI count',ylab='',yaxt='n',
       main=paste0('nUMI distribution (>=',minUMI,' UMIs, ',nrow(df.overMin),' beads)'),
       sub=paste0('median: ',median(df.overMin$nUMI),'   mean: ',round(mean(df.overMin$nUMI),2)))
  polygon(den,col='#0088ff44')
  abline(v=median(log10(df.overMin$nUMI)))
  abline(v=mean(log10(df.overMin$nUMI)),lty='dashed')
  
  den <- density(log10(df.overMin$nGene))
  plot(den,
       col='blue',lwd=2,
       xlab='log10 feature count',ylab='',yaxt='n',
       main=paste0('nFeature distribution (>=',minUMI,' UMIs, ',nGenes.overMin,' unique features)'),
       sub=paste0('median: ',median(df.overMin$nGene),'   mean: ',round(mean(df.overMin$nGene),2)))
  polygon(den,col='#0088ff44')
  abline(v=median(log10(df.overMin$nGene)))
  abline(v=mean(log10(df.overMin$nGene)),lty='dashed')
  
  
  plot(df.overMin$x~df.overMin$y,col=df.overMin$col,
       pch=20,cex=.2,xlab='',ylab='',main=paste0('nUMI (>=',minUMI,' UMIs)'),asp=1)
  
  shh <- dev.off()
  
}
