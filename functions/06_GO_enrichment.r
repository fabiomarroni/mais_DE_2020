
# Run with --help or -h flag for help.
# Written 12/31/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-G", "--GOfile"), type="character", default="",
	help="Formatted GO file (one transcript pe rline, with GO terms separated by semi-colons) [default= %default]", metavar="character"), 
  make_option(c("-I", "--infile"), type="character", default="",
	help="DE results, already annotated with GO terms [default= %default]", metavar="character"), 
  make_option(c("-N", "--ontology"), type="character", default="BP", 
	help="GO ontology [default= %default]", metavar="character"),
  make_option(c("-A", "--algorithm"), type="character", default="weight01", 
	help="Algorithm for testing enrichment [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", 
	help="output file containing results of the enrichment test [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$GOfile)) {
  stop("WARNING: No GOfile specified with '-G' flag.")
} else {  cat ("GOfile ", opt$GOfile, "\n")
  GOfile <- opt$GOfile  
  }

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-I' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$ontology)) {
  stop("WARNING: No ontology specified with '-N' flag.")
} else {  cat ("ontology ", opt$ontology, "\n")
  ontology <- opt$ontology  
  }

if (is.null(opt$algorithm)) {
  stop("WARNING: No algorithm specified with '-A' flag.")
} else {  cat ("algorithm ", opt$algorithm, "\n")
  algorithm <- opt$algorithm  
  }

  
  
if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

topgosel<-function(x)
{
y<-x>0
return(y)
}

#This function is prepared for reading cuffdiff inpunt
GO.enrich<-function(GOfile,infile,outfile,algorithm,ontology,top.res=50)
{
library(topGO)
library(data.table)
myres<-fread(infile,data.table=F)
#List of all the genes obtained in the De results (specific for fungi or plants)
allgenes<-unique(unlist(strsplit(myres$Gene_id,";")))
#Significantly DE genes
myInterestingGenes<-unique(unlist(strsplit(myres$Gene_id[myres$significant=="yes"],";")))
GO<-fread(GOfile,data.table=F)
#From the universe of GO, only extract genes to which our transcripts match
GO<-GO[GO$Gene_id%in%allgenes,]
fullset <- as.integer(allgenes %in% myInterestingGenes)
fullset<-as.numeric(fullset)
names(fullset) <- allgenes
#Create the gene2GO object
mygene2GO<-as.list(strsplit(GO$GO,";"))
mygene2GO<-setNames(mygene2GO,GO$Gene_id)
#Finally create this topGO object
mytopgo<-new("topGOdata",ontology=ontology,allGenes=fullset,annot = annFUN.gene2GO, geneSel= topgosel , gene2GO = mygene2GO)
sampleGOdata<-mytopgo
#default algorithm of topGO, weight01
#resultFisher <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
#New data: algorithm is now specified as a parameter
resultFisher <- runTest(sampleGOdata, algorithm = algorithm, statistic = "fisher")
#Extract Fisher results
allRes <- GenTable(sampleGOdata, pvalue = resultFisher, orderBy = "resultFisher", topNodes = length(resultFisher@score))
#Compute FDR
allRes$qvalue<-p.adjust(allRes$pvalue,method="BH")
allRes<-allRes[order(allRes$qvalue),]
#Only select top 50 results
allRes<-allRes[1:top.res,]
#Remove non-significant results
allRes<-allRes[allRes$qvalue<=0.05,]
if(nrow(allRes>0)) write.table(allRes,outfile,row.names=F,quote=F,sep="\t")
}
GO.enrich(GOfile=GOfile,infile=infile,algorithm=algorithm,ontology=ontology,outfile=outfile)



