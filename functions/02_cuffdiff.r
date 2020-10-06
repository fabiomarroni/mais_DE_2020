# Run with --help flag for help.
# Modified 12/30/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input_file"), type="character", default="",
              help="Input file (cuffdiff output) [default= %default]", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$input_file)) {
  stop("WARNING: No input file specified with '-I' flag.")
} else {  cat ("Input file is ", opt$input_file, "\n")
  input_file <- opt$input_file  
  }

  if (is.null(opt$out)) {
  stop("WARNING: No input directory specified with '-I' flag.")
} else {  cat ("Output file is ", opt$out, "\n")
  outfile <- opt$out  
  }


  
readcuffdiff<-function(input_file, outfile) 
{
	library(data.table)
	diff<-fread(input_file,data.table=F)	
	diff<-diff[diff$status!="NOTEST",]	
	genestosel<-paste("gene",c("Zm00001d037242","Zm00001d027557","Zm00001d027582","Zm00001d008587","Zm00001d048979"),sep=":")	
	seldiff<-diff[diff$gene_id%in%genestosel,]
	write.table(seldiff,outfile,row.names=F,quote=F,sep="\t")	
}
readcuffdiff(input_file=input_file,outfile=outfile)
