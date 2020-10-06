#!/bin/bash

#Variables you need to set
BASEDIR=/path/to/your/results #Exiting folder that you want to be the top folder of your results
FUNCDIR=/path/to/R/functions #Folder in which you will store the R function provided
##############################



module load sw/bio/cufflinks/2.2.0
GTF=/genomes/zea_mays/genes/assembly_build_v4/v4.37/annotation/Zea_mays.AGPv4.37.gtf
BAMDIR=${BASEDIR}/06_hisat2/

RESDIR=${BASEDIR}/08_DE_cuffdiff
mkdir -p $RESDIR/logs
cuffdiff -L N,T_A,T_C,T_E --dispersion-method pooled --num-threads 8  --library-type fr-firststrand --library-norm-method geometric \
--output-dir ${RESDIR} $GTF ${BAMDIR}/45754_ID1462_1-1ZANIN_S1_sort.bam,${BAMDIR}/45755_ID1462_2-2ZANIN_S2_sort.bam,${BAMDIR}/45756_ID1462_3-3ZANIN_S3_sort.bam\
 ${BAMDIR}/45757_ID1462_4-4ZANIN_S4_sort.bam,${BAMDIR}/45758_ID1462_5-5ZANIN_S5_sort.bam,${BAMDIR}/45759_ID1462_6-6ZANIN_S6_sort.bam\
 ${BAMDIR}/45760_ID1462_7-7ZANIN_S7_sort.bam,${BAMDIR}/45761_ID1462_8-8ZANIN_S8_sort.bam,${BAMDIR}/45762_ID1462_9-9ZANIN_S9_sort.bam\
 ${BAMDIR}/45763_ID1462_10-10ZANIN_S10_sort.bam,${BAMDIR}/45764_ID1462_11-11ZANIN_S11_sort.bam,${BAMDIR}/45765_ID1462_12-12ZANIN_S12_sort.bam > ${RESDIR}/logs/log.out 2> ${RESDIR}/logs/log.err

#Individual comparisons

RESDIR=${BASEDIR}/08_DE_cuffdiff/N_vs_TA
mkdir -p $RESDIR
cuffdiff -L N,T_A --dispersion-method pooled --num-threads 8  --library-type fr-firststrand --library-norm-method geometric \
--output-dir ${RESDIR} $GTF ${BAMDIR}/45754_ID1462_1-1ZANIN_S1_sort.bam,${BAMDIR}/45755_ID1462_2-2ZANIN_S2_sort.bam,${BAMDIR}/45756_ID1462_3-3ZANIN_S3_sort.bam ${BAMDIR}/45757_ID1462_4-4ZANIN_S4_sort.bam,${BAMDIR}/45758_ID1462_5-5ZANIN_S5_sort.bam,${BAMDIR}/45759_ID1462_6-6ZANIN_S6_sort.bam

RESDIR=${BASEDIR}/08_DE_cuffdiff/N_vs_TC
mkdir -p $RESDIR
cuffdiff -L N,T_C --dispersion-method pooled --num-threads 8  --library-type fr-firststrand --library-norm-method geometric \
--output-dir ${RESDIR} $GTF ${BAMDIR}/45754_ID1462_1-1ZANIN_S1_sort.bam,${BAMDIR}/45755_ID1462_2-2ZANIN_S2_sort.bam,${BAMDIR}/45756_ID1462_3-3ZANIN_S3_sort.bam ${BAMDIR}/45760_ID1462_7-7ZANIN_S7_sort.bam,${BAMDIR}/45761_ID1462_8-8ZANIN_S8_sort.bam,${BAMDIR}/45762_ID1462_9-9ZANIN_S9_sort.bam

RESDIR=${BASEDIR}/08_DE_cuffdiff/N_vs_TE
mkdir -p $RESDIR
cuffdiff -L N,T_E --dispersion-method pooled --num-threads 8  --library-type fr-firststrand --library-norm-method geometric \
--output-dir ${RESDIR} $GTF ${BAMDIR}/45754_ID1462_1-1ZANIN_S1_sort.bam,${BAMDIR}/45755_ID1462_2-2ZANIN_S2_sort.bam,${BAMDIR}/45756_ID1462_3-3ZANIN_S3_sort.bam ${BAMDIR}/45763_ID1462_10-10ZANIN_S10_sort.bam,${BAMDIR}/45764_ID1462_11-11ZANIN_S11_sort.bam,${BAMDIR}/45765_ID1462_12-12ZANIN_S12_sort.bam

RESDIR=${BASEDIR}/08_DE_cuffdiff/TA_vs_TC
mkdir -p $RESDIR
cuffdiff -L T_A,T_C --dispersion-method pooled --num-threads 8  --library-type fr-firststrand --library-norm-method geometric \
--output-dir ${RESDIR} $GTF ${BAMDIR}/45757_ID1462_4-4ZANIN_S4_sort.bam,${BAMDIR}/45758_ID1462_5-5ZANIN_S5_sort.bam,${BAMDIR}/45759_ID1462_6-6ZANIN_S6_sort.bam ${BAMDIR}/45760_ID1462_7-7ZANIN_S7_sort.bam,${BAMDIR}/45761_ID1462_8-8ZANIN_S8_sort.bam,${BAMDIR}/45762_ID1462_9-9ZANIN_S9_sort.bam

RESDIR=${BASEDIR}/08_DE_cuffdiff/TA_vs_TE
mkdir -p $RESDIR
cuffdiff -L T_A,T_E --dispersion-method pooled --num-threads 8  --library-type fr-firststrand --library-norm-method geometric \
--output-dir ${RESDIR} $GTF ${BAMDIR}/45757_ID1462_4-4ZANIN_S4_sort.bam,${BAMDIR}/45758_ID1462_5-5ZANIN_S5_sort.bam,${BAMDIR}/45759_ID1462_6-6ZANIN_S6_sort.bam\
 ${BAMDIR}/45763_ID1462_10-10ZANIN_S10_sort.bam,${BAMDIR}/45764_ID1462_11-11ZANIN_S11_sort.bam,${BAMDIR}/45765_ID1462_12-12ZANIN_S12_sort.bam

RESDIR=${BASEDIR}/08_DE_cuffdiff/TC_vs_TE
mkdir -p $RESDIR
cuffdiff -L T_C,T_E --dispersion-method pooled --num-threads 8  --library-type fr-firststrand --library-norm-method geometric \
--output-dir ${RESDIR} $GTF ${BAMDIR}/45760_ID1462_7-7ZANIN_S7_sort.bam,${BAMDIR}/45761_ID1462_8-8ZANIN_S8_sort.bam,${BAMDIR}/45762_ID1462_9-9ZANIN_S9_sort.bam ${BAMDIR}/45763_ID1462_10-10ZANIN_S10_sort.bam,${BAMDIR}/45764_ID1462_11-11ZANIN_S11_sort.bam,${BAMDIR}/45765_ID1462_12-12ZANIN_S12_sort.bam


#Write only the output for the selected genes on which we focused (gene names are written in the function)
for aaa in N_vs_TA N_vs_TC N_vs_TE TA_vs_TC TA_vs_TE TC_vs_TE
do
Rscript ${FUNCDIR}/02_cuffdiff.r -I ${BASEDIR}/08_DE_cuffdiff/${aaa}/gene_exp.diff -O ${BASEDIR}/08_DE_cuffdiff/${aaa}/${aaa}_selected.txt
done
