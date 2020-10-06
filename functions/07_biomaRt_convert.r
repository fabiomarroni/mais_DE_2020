# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
	help="Input file containing DE and GO file (output of 04_merge_blast_GO_DE.r) [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
	help="output file with entrez gene names and eventually uniprot names [default= %default]", metavar="character")
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

retrieve.biomaRt<-function(host="plants.ensembl.org",org="mays",infile,outfile)
{
library(data.table)
library(biomaRt)
mylist<-listMarts(host=host)
mysets<-listDatasets(useMart("plants_mart", host = host))
myorg<-grep(org,mysets$description,ignore.case = T,value=T)
#Grep organism name in the description and select the corresponding dataset name
mydataset<-mysets$dataset[grep(org,mysets$description,ignore.case = T,value=F)]
allattributes<-listAttributes(useDataset(dataset = as.character(mydataset), mart    = useMart("plants_mart",host = "plants.ensembl.org")))
allfilters<-listFilters(useDataset(dataset = as.character(mydataset), mart    = useMart("plants_mart",host = "plants.ensembl.org")))
#Use a given dataset for analysis
myusemart <- useDataset(as.character(mydataset), mart = useMart("plants_mart", host = host))
#Get all genes and their annotations
cat("Retrieving genes and annotations...\n")
#resultTable <- getBM(attributes=c("ensembl_gene_id","entrezgene","start_position","end_position","go_id","goslim_goa_accession","goslim_goa_description","kegg_enzyme","uniprotswissprot"), mart = myusemart)
#Since at present I don't really need all kegg entries and GOs, I just download gene name conversion to avoid timeout errors.
resultTable <- getBM(attributes=c("ensembl_gene_id","entrezgene_id","start_position","end_position","uniprotswissprot"), mart = myusemart)
myindata<-fread(infile,data.table=F)
#Myoutdata may have some more lines than the input data because some ensembl_gene_id are associated to more than one entrez_gene_id
#It's not my fault!!!
myoutdata<-merge(myindata,resultTable,by.x="test_id",by.y="ensembl_gene_id",all.x=T,sort=F)
write.table(myoutdata,outfile,sep="\t",quote=F,row.names=F)
}
retrieve.biomaRt(infile=infile,outfile=outfile)
