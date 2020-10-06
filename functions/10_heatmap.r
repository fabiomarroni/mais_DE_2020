# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
  help="Input excel file with significant genes (custom format provided by Laura) [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
  help="output graph with heatmap [default= %default]", metavar="character")
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

heat.plot<-function(infile,outfile)
{
library(openxlsx)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
pp<-read.xlsx(infile)
pp<-pp[,c("gene_id","N","TA","TC","TE")]
row.names(pp)<-pp$gene_id
pp$gene_id<-NULL
pdf(outfile)
par(mar=c(1,1,1,1))
cc<-Heatmap(pp, 
	show_row_names=FALSE,
	column_names_gp = gpar(fontsize = 8),
	name="FPKM",
	width=unit(9,"cm"),
	height=unit(13,"cm"),
	show_heatmap_legend = TRUE)
draw(cc,adjust_annotation_extension=TRUE,padding=unit(c(0.1,0.1,0.1,0.1),"mm"))
#draw(pd,x=unit(0.92,"npc"),y=unit(0.5,"npc"))
dev.off()

pp<-log10(0.1+pp)
pdf(gsub("FPKM","log10FPKM",outfile))
par(mar=c(1,1,1,1))
cc<-Heatmap(pp, 
	show_row_names=FALSE,
	column_names_gp = gpar(fontsize = 8),
	name="Log10\nFPKM",
	width=unit(9,"cm"),
	height=unit(13,"cm"),
	show_heatmap_legend = TRUE)
draw(cc,adjust_annotation_extension=TRUE,padding=unit(c(0.1,0.1,0.1,0.1),"mm"))
#draw(pd,x=unit(0.92,"npc"),y=unit(0.5,"npc"))
dev.off()

}
heat.plot(infile=infile,outfile=outfile)
