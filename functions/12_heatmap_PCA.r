# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
  help="Input file containing FPKM counts per replicate (output of 11_FPKM_by_rep.r) [default= %default]", metavar="character"), 
  make_option(c("-M", "--metafile"), type="character", default="",
  help="Metafile [default= %default]", metavar="character"), 
  make_option(c("-F", "--FPKM"), type="numeric", default=5,
  help="Min FPKM to retain transcript [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="",
  help="output graph file (should contain the strings 'heatmap_' and 'FPKM' to allow name changing for additional graphs. [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-I' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-M' flag.")
} else {  cat ("metafile ", opt$metafile, "\n")
  metafile <- opt$metafile  
  }

if (is.null(opt$FPKM)) {
  stop("WARNING: No FPKM specified with '-F' flag.")
} else {  cat ("FPKM ", opt$FPKM, "\n")
  FPKM <- opt$FPKM  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-G' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

heat.plot<-function(infile,metafile,FPKM,outfile,top=500)
{
library(openxlsx)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(cluster)
library(ggfortify)

metadata<-fread(metafile,data.table=F,header=F)
row.names(metadata)<-metadata[,1]
metadata[,1]<-NULL
pp<-fread(infile,data.table=F)
row.names(pp)<-pp$gene_id
pp$gene_id<-NULL
pp$mean<-apply(pp,1,"mean")
pp<-pp[pp$mean>=FPKM,]
pp<-pp[order(pp$mean,decreasing=T),]
pp$mean<-NULL
mypca<-prcomp(t(pp))

scores <- as.data.frame(mypca$x)           
expvar1<-round(100*mypca$sdev[1]^2/sum(mypca$sdev^2),0)
expvar2<-round(100*mypca$sdev[2]^2/sum(mypca$sdev^2),0)

ggdata <- merge(scores, metadata,by="row.names")
setnames(ggdata,"V2","Treatment")
outgraph<-gsub("heatmap_","PCA_",outfile)
pdf(outgraph)
gpgo<-ggplot(ggdata) +
  geom_point(aes(x=PC1, y=PC2, color=factor(Treatment)), size=5, shape=20) +
  xlab( paste ("PC1 (",expvar1,"% Explained variance)",sep="")) +
  ylab( paste ("PC2 (",expvar2,"% Explained variance)",sep="")) +
  # stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
               # geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Treatment"),fill=guide_legend("Treatment"))
print(gpgo)
dev.off()
mycol<-list(Treatment=c("N"="red","T_A"="blue","T_C"="darkgreen","T_E"="purple"))
trt<-HeatmapAnnotation(Treatment=ggdata[,"Treatment"],col=mycol,show_legend=F,show_annotation_name=F)
legtrt<-Legend(labels=sort(unique(ggdata[,"Treatment"])),title="Treatment",legend_gp=gpar(fill=mycol$Treatment),
		grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
		labels_gp = gpar(fontsize = 6),title_gp = gpar(fontsize = 7, fontface = "bold"))
maincolful<-circlize::colorRamp2(c(min(pp),mean(unlist(pp)), max(pp)), c("blue", "white", "red"))
at<-100*round(seq(from=min(pp),to=max(pp),length.out=5)/100,0)
legheat<-Legend(labels=at,title="FPKM",legend_gp=gpar(fill=maincolful(at)),
		grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
		labels_gp = gpar(fontsize = 6),title_gp = gpar(fontsize = 7, fontface = "bold"))
pd = packLegend(legtrt, legheat)

names(pp)<-unlist(lapply(strsplit(names(pp),"_"),"[",4))
#pp<-pp[1:top,]
pdf(outfile)
par(mar=c(1,1,1,1))
cc<-Heatmap(pp, 
	show_row_names=FALSE,
	column_names_gp = gpar(fontsize = 8),
	name="FPKM",
	width=unit(9,"cm"),
	height=unit(13,"cm"),
	top_annotation=trt,
	cluster_rows = FALSE,
	show_heatmap_legend = FALSE)
draw(cc,adjust_annotation_extension=TRUE,padding=unit(c(0.1,0.1,0.1,0.1),"mm"))
draw(pd,x=unit(0.92,"npc"),y=unit(0.5,"npc"))
dev.off()

pp<-log10(0.1+pp)
maincolful<-circlize::colorRamp2(c(min(pp),mean(unlist(pp)), max(pp)), c("blue", "white", "red"))
at<-round(seq(from=min(pp),to=max(pp),length.out=5),2)
legheat<-Legend(labels=at,title=bquote(log[10]~"FPKM"),legend_gp=gpar(fill=maincolful(at)),
		grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
		labels_gp = gpar(fontsize = 6),title_gp = gpar(fontsize = 7, fontface = "bold"))
pd = packLegend(legtrt, legheat)

pdf(gsub("FPKM","log10FPKM",outfile))
par(mar=c(1,1,1,1))
cc<-Heatmap(pp, 
	show_row_names=FALSE,
	column_names_gp = gpar(fontsize = 8),
	name="Log10\nFPKM",
	width=unit(9,"cm"),
	height=unit(13,"cm"),
	top_annotation=trt,
	cluster_rows = FALSE,
	show_heatmap_legend = FALSE)
draw(cc,adjust_annotation_extension=TRUE,padding=unit(c(0.1,0.1,0.1,0.1),"mm"))
draw(pd,x=unit(0.92,"npc"),y=unit(0.5,"npc"))
dev.off()

}
heat.plot(infile=infile,outfile=outfile,metafile=metafile,FPKM=FPKM)
