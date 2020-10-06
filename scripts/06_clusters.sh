#!/bin/bash

#Variables you need to set
BASEDIR=/path/to/your/results #Exiting folder that you want to be the top folder of your results
FUNCDIR=/path/to/R/functions #Folder in which you will store the R function provided
##############################

#Plot clusters of genes with similar expression
#The plan is to use mfuzz (https://2-bitbio.com/post/fuzzy-cmeans-clustering-of-rnaseq-data-using-mfuzz/)
Rscript ${FUNCDIR}/09_clusters_plot.r \
-I ${BASEDIR}/11_clusters/Clustering_RNAseq2019.xlsx
-K TRUE \
-O ${BASEDIR}/11_clusters/fuzzy_clust_RNAseq2019.pdf \
-T ${BASEDIR}/11_clusters/fuzzy_clust_RNAseq2019.txt

Rscript ${FUNCDIR}/09_clusters_plot.r \
-I ${BASEDIR}/11_clusters/Clustering_RNAseq2019.xlsx \
-K FALSE \
-O ${BASEDIR}/11_clusters/fuzzy_clust_RNAseq2019_noN.pdf \
-T ${BASEDIR}/11_clusters/fuzzy_clust_RNAseq2019_noN.txt


Rscript ${FUNCDIR}/10_heatmap.r \
-I ${BASEDIR}/11_clusters/RNAseq-maize2019_only significant.xlsx \
-O ${BASEDIR}/11_clusters/heatmap_FPKM.pdf


#Create a table with FPKM of every single replicate (no filtering)
Rscript ${FUNCDIR}/11_FPKM_by_rep.r \
-I ${BASEDIR}/06_hisat2 \
-M ${BASEDIR}/06_hisat2/samplecondition.txt \
-O ${BASEDIR}/11_clusters/FPKM_rep_table.txt

#Plot PCA and some (presumably old) heatmaps. The "official" heatmap is produced in function 10 using complexHeatmaps
Rscript ${FUNCDIR}/12_heatmap_PCA.r \
-I ${BASEDIR}/11_clusters/FPKM_rep_table.txt \
-M ${BASEDIR}/06_hisat2/samplecondition.txt \
-O ${BASEDIR}/11_clusters/heatmap_per_rep_FPKM.pdf
























