# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="",
  help="Parent input directory for hisat2 results [default= %default]", metavar="character"), 
  make_option(c("-M", "--metadata"), type="character", default="",
  help="Metadata [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
  help="output file (FPKM counts per replicates) [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$indir)) {
  stop("WARNING: No indir specified with '-I' flag.")
} else {  cat ("indir ", opt$indir, "\n")
  indir <- opt$indir  
  }

if (is.null(opt$metadata)) {
  stop("WARNING: No metadata specified with '-M' flag.")
} else {  cat ("metadata ", opt$metadata, "\n")
  metadata <- opt$metadata  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-G' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

FPKM.table<-function(indir,metadata,outfile)
{
library(data.table)
meta<-fread(metadata,data.table=F,header=F)
myfiles<-paste(indir,"/",meta[,1],"/",meta[,1],".gtf",sep="")
#Trick to skip "#" and "exon" at the same time, since fread doesn't have comment.char
myfpkm<-fread(cmd=paste("grep FPKM", myfiles[1]),data.table=F)
myfpkm$gene_id<-unlist(lapply(strsplit(myfpkm$V9,"\""),"[",2))
myfpkm$gene_id<-gsub("gene:","",myfpkm$gene_id)
myfpkm$FPKM<-as.numeric(as.character(unlist(lapply(strsplit(myfpkm$V9,"\""),"[",8))))
myfpkm<-myfpkm[,c("gene_id","FPKM")]
myfpkm<-aggregate(myfpkm$FPKM,by=list(myfpkm$gene_id),FUN="sum")
setnames(myfpkm,c("gene_id",meta[1,1]))
finalfpkm<-myfpkm
for(aaa in 2:length(myfiles))
{
cat(aaa,"out of",length(myfiles),"\n")
myfpkm<-fread(cmd=paste("grep FPKM", myfiles[aaa]),data.table=F)
myfpkm$gene_id<-unlist(lapply(strsplit(myfpkm$V9,"\""),"[",2))
myfpkm$gene_id<-gsub("gene:","",myfpkm$gene_id)
myfpkm$FPKM<-as.numeric(as.character(unlist(lapply(strsplit(myfpkm$V9,"\""),"[",8))))
myfpkm<-myfpkm[,c("gene_id","FPKM")]
myfpkm<-aggregate(myfpkm$FPKM,by=list(myfpkm$gene_id),FUN="sum")
setnames(myfpkm,c("gene_id",meta[aaa,1]))
finalfpkm<-merge(finalfpkm,myfpkm,by="gene_id",all=T)
cat("size=",nrow(finalfpkm),"\n")
}

write.table(finalfpkm,outfile,sep="\t",row.names=F,quote=F)

}
FPKM.table(indir=indir,metadata=metadata,outfile=outfile)
