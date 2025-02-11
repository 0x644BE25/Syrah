######################################################
# DETERMINE VERSION
#
# GOAL: Use the first 100K read 1 sequences to
# calculate nucleotide frequency by position and write
# to file. Use frequencies to determine the capture 
# oligo version and appropriate barcode extraction 
# parameters and write those out.
# Output nucleotide frequency by position plot with
# expected barcode/UMI positions annotated.
######################################################

# ================= IMPORTS ==========================

# ================= PARAMS ===========================

manifestFile <- commandArgs(trailingOnly=TRUE)[1]
true <- TRUE
false <- FALSE
source(manifestFile)
if (endsWith(writeDir,'/')) { writeDir <- paste0(writeDir,'/') }

nts <- c('A','C','G','T','N')


# ================= INIT DATA ========================

r1s <- readLines(read1fastq,n=400000)

# ================= ANALYSIS =========================

# READ 1 NUCLEOTIDE FREQUENCY
r1s <- r1s[seq(2,400000,by=4)]
r1s <- do.call(rbind,strsplit(r1s,''))
freqs <- apply(r1s,2,function(x){ table(factor(x,levels=nts))/length(x) })
colnames(freqs) <- 1:ncol(freqs)
write.csv(freqs,paste0(writeDir,batchName,'_estimated_read_1_nucleotide_frequencies.csv'),row.names=TRUE)

# VERSION INFO
if (freqs['T',34]>.3 & freqs['C',35]>.3) {
  vs <- 'yesTC' 
} else if (freqs['T',40]<.25 & freqs['T',41]<.25) {
  vs <- 'curio'
} else {
  vs <- 'noTC' 
}
writeLines(vs,paste0(writeDir,batchName,'_r1_version.txt'))

# NON-SYRAH BARCODE/UMI EXTRACTION PATTERN
link2len <- c('yesTC'=3,'noTC'=1,'curio'=1)[vs]
umilen <- c('yesTC'=8,'noTC'=8,'curio'=6)[vs]
patt <- c(rep('C',8),rep('X',18),rep('C',6),rep('X',link2len),rep('N',umilen),rep('X',99))
patt <- paste0(patt[1:ncol(r1s)],collapse='')
writeLines(patt,paste0(writeDir,batchName,'_r1_pattern.txt'))

# ================= PLOTTING =========================

# FOR READ1 STRUCTURE PLOT
spaces <- list(yesTC=c(rep(0,8),.25,rep(0,17),.25,rep(0,5),.25,rep(0,2),.25,rep(0,7),.25,rep(0,99)),
               noTC=c(rep(0,8),.25,rep(0,17),.25,rep(0,5),.25,.25,rep(0,7),.25,rep(0,99)),
               curio=c(rep(0,8),.25,rep(0,17),.25,rep(0,5),.25,.25,rep(0,5),.25,rep(0,1),.25,rep(0,99)))
spc <- spaces[[vs]][1:ncol(freqs)]

chars <- c('X'='','C'='B','N'='U')[strsplit(patt,'')[[1]]]
colnames(freqs) <- paste0(chars,'\n',colnames(freqs))

png(filename=paste0(writeDir,batchName,'_estimated_read_1_nucleotide_frequencies.png'),width=(2.5+(.3*ncol(freqs))),height=6,units='in',res=100)
barplot(as.matrix(freqs[c('T','G','C','A','N'),]), 
        col=c('#ff585f','#ffcb05','#73ab1c','#496efd','#CCCCCC'),
        border=NA,
        space=spaces[[vs]][1:ncol(freqs)],
        legend.text=c('T','G','C','A','N'),
        yaxt='n',
        cex.names=.75,
        args.legend=list(y=1,x=-.1,box.col='transparent'),
        offset=0,
        main=paste(batchName,'estimated read 1 nucleotide frequency by position'),
        xlab="5' ----------- nt position ---------> 3'")

dev.off()
