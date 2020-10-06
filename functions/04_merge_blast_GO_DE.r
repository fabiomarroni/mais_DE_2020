# Run with --help or -h flag for help.
# Written 12/31/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-G", "--gofile"), type="character", default="",
	help="Formatted GO file (output of function 03_reformat_GO.r) [default= %default]", metavar="character"), 
  make_option(c("-B", "--blastfile"), type="character", default="",
	help="Output of blastx of genes read from DE and associated with GO [default= %default]", metavar="character"), 
  make_option(c("-N", "--name2tag"), type="character", default="",
	help="File in which both the gene name and gene tag are reported. Fasta file has both, blast only report gene name and we need the gene tag [default= %default]", metavar="character"), 
  make_option(c("-D", "--defile"), type="character", default="",
	help="DE result file (cuffdiff) [default= %default]", metavar="character"), 
  make_option(c("-O", "--out"), type="character", default="", 
	help="output file (contains DE with GO) [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$gofile)) {
  stop("WARNING: No GO file specified with '-G' flag.")
} else {  cat ("DE file ", opt$gofile, "\n")
  gofile <- opt$gofile  
  }

if (is.null(opt$blastfile)) {
  stop("WARNING: No blast file specified with '-B' flag.")
} else {  cat ("Blast file ", opt$blastfile, "\n")
  blastfile <- opt$blastfile  
  }

if (is.null(opt$name2tag)) {
  stop("WARNING: No name2tag file specified with '-N' flag.")
} else {  cat ("name2tag file ", opt$name2tag, "\n")
  name2tag <- opt$name2tag  
  }

if (is.null(opt$defile)) {
  stop("WARNING: No DE file specified with '-D' flag.")
} else {  cat ("DE file ", opt$defile, "\n")
  defile <- opt$defile  
  }

if (is.null(opt$out)) {
  stop("WARNING: No output file '-O' flag.")
} else {  cat ("Output file is ", opt$out, "\n")
  outfile <- opt$out  
  #setwd(wd_location)  
  }

reformat.go<-function(gofile,blastfile,name2tag,defile,outfile,minfrac=0.01)
{
library(data.table)
go<-fread(gofile,data.table=F)
blast<-fread(blastfile,data.table=F)
n2tag<-fread(name2tag,data.table=F,sep=" ",header=F)
n2tag<-n2tag[!duplicated(n2tag[,2]),]
#Only keep gene names and gene descritpions in the blast results.
blast<-blast[,c(1,2,13)]
blast<-blast[!duplicated(blast),]
blast$Gene_id<-unlist(lapply(strsplit(unlist(lapply(strsplit(blast$V2,"\\|"),"[",4)),"\\."),"[",1))
blast$V2<-NULL
blast<-merge(blast,n2tag,by="V1",sort=F)
blast$gene_name<-unlist(lapply(strsplit(unlist(lapply(strsplit(blast$V13,"RecName: Full="),"[",2)),";"),"[",1))
setnames(blast,"V2","Gene_tag")
blastgo<-merge(go,blast,by="Gene_id",all=T)
blastgo$Gene_tag<-gsub("ZEAMMB73_","",blastgo$Gene_tag)
DE<-fread(defile,data.table=F)
DE$test_id<-gsub("gene:","",DE$test_id)
dego<-merge(DE,blastgo,by.x="test_id",by.y="Gene_tag",all.x=T,sort=F)
write.table(dego,outfile,quote=F,row.names=F,sep="\t")
}
reformat.go(gofile=gofile,blastfile=blastfile,name2tag=name2tag,defile=defile,outfile=outfile)
