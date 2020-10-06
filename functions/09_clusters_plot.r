# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
	help="Input Excel file for clustering (custom file) [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
	help="output graph of clusters [default= %default]", metavar="character"),
  make_option(c("-K", "--keepN"), type="logical", default=FALSE, 
	help="Sholud N (the control) be kept? [default= %default]", metavar="character"),
  make_option(c("-T", "--outtable"), type="character", default="", 
	help="output table of hard cluster assignment [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-G' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-G' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

if (is.null(opt$keepN)) {
  stop("WARNING: No keepN specified with '-K' flag.")
} else {  cat ("keepN ", opt$keepN, "\n")
  keepN <- opt$keepN  
  }

if (is.null(opt$outtable)) {
  stop("WARNING: No outtable specified with '-T' flag.")
} else {  cat ("outtable ", opt$outtable, "\n")
  outtable <- opt$outtable  
  }

cfuzz<-function(infile,outfile,outtable,keepN)
{
library(openxlsx)
library(Mfuzz)
library(edgeR)
library(SummarizedExperiment)
library(data.table)
pp<-read.xlsx(infile)
#Try to use Mfuzz
#Tutorial here:
#https://2-bitbio.com/post/fuzzy-cmeans-clustering-of-rnaseq-data-using-mfuzz/
setnames(pp,c("RPKM","X6","X7","X8"),c("N","TA","TC","TE"))

ifelse(keepN,toc<-pp[,c("gene.ID","N","TA","TC","TE")],toc<-pp[,c("gene.ID","TA","TC","TE")])
toc<-toc[2:nrow(toc),]
row.names(toc)<-toc[,1]
toc[,1]<-NULL
if(keepN) toc$N<-as.numeric(toc$N)
toc$TA<-as.numeric(toc$TA)
toc$TC<-as.numeric(toc$TC)
toc$TE<-as.numeric(toc$TE)
#Remove rows with zero variance, because they will result in an error when clustering
#NA/NaN/Inf in foreign function call (arg 1)
myvar<-apply(toc,1,var)
toc<-toc[myvar>0,]

timepoint <- seq(1,ncol(toc))
toc<-rbind(timepoint,toc)
row.names(toc)[1]<-"time"

#save it to a temp file so ti doesnt clutter up my blog directory
tmp <- tempfile()
write.table(toc,file=tmp, sep='\t', quote = F, col.names=NA)

#read it back in as an expression set
data <- table2eset(file=tmp)
#Standardise data 
data.s <- standardise(data)
#Estimate the "fuzzifier", m1
m1 <- mestimate(data.s)

#Select optimal number of clusters (we use the elbow rule, I set it as when the decrease in centroid means is less that 0.1)
set.seed(1)
knumb<-Dmin(data.s, m=m1, crange=seq(2,15,1), repeats=20, visu=FALSE)
dk<-round(abs(diff(knumb)),3)
clust<-1+max(which(dk>0.1))
maremma <- mfuzz(data.s,c=clust,m=m1)
if(clust==2) 
{
myrow<-1
mycol<-2
mywidth<-14
myheight<-7

}
if(clust>2&clust<=4) 
{
myrow<-2
mycol<-2
mywidth<-7
myheight<-7
}
if(clust>4) 
{
myrow<-3
mycol<-3
mywidth<-7
myheight<-7
}
pdf(outfile,width=mywidth,height=myheight)
mfuzz.plot2(data.s,cl=maremma,mfrow=c(myrow,mycol),xlab="Treatment",time.labels=names(toc),x11=FALSE,centre=TRUE)
dev.off()
pdf(gsub(".pdf","_fancy.pdf",outfile),width=mywidth,height=myheight)
mfuzz.plot2(data.s,cl=maremma,mfrow=c(myrow,mycol),xlab="Treatment",time.labels=names(toc),colo="fancy",x11=FALSE,centre=TRUE)
dev.off()
write.table(maremma$cluster,outtable,col.names=FALSE,quote=F,sep="\t")
}
cfuzz(infile=infile,outfile=outfile,outtable=outtable,keepN=keepN)
