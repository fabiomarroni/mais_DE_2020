# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
	help="Formatted GO file with entrez codes (output of 07_biomaRt_convert.r) [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
	help="output KEGG enrichment file [default= %default]", metavar="character")
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

kegg.enrich<-function(organism="zma",infile,outfile,cutoff=0.05)
{
library(data.table)
library(clusterProfiler)
library(DOSE)
mydat<-fread(infile,data.table=F)
#Only keep significant genes (because only using that we will comput enrichment).
# The universe is determined based on all the genes.
mydat<-mydat[mydat$significant=="yes",]
mykegg<-enrichKEGG(mydat$entrezgene,organism=organism,pvalueCutoff=cutoff,qvalueCutoff=cutoff)
newkegg<-as.data.frame(mykegg)
if(nrow(newkegg)>0)
{
	newkegg$V4_ID<-NA
	#For the sake of completeness, I also attach the V4 gene names to the results
	for(aaa in 1:nrow(newkegg))
	{
	myentrez<-data.frame(unlist(strsplit(newkegg$geneID[aaa],"/")))
	names(myentrez)<-"entrezgene_id"
	tdf<-merge(myentrez,mydat,by="entrezgene_id",all=F)
	myzm<-paste(unique(tdf$test_id),sep="/",collapse="/")
	newkegg$V4_ID[aaa]<-myzm
	}
	write.table(newkegg,outfile,sep="\t",quote=F,row.names=F)
}
}
kegg.enrich(infile=infile,outfile=outfile)
