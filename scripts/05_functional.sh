#!/bin/bash
#Variables you need to set
BASEDIR=/path/to/your/results #Exiting folder that you want to be the top folder of your results
FUNCDIR=/path/to/R/functions #Folder in which you will store the R function provided
BLASTDIR=/path/to/blast/database #Path to the blast local database
LOGDIR=/path/to/log/files #Path to folder in which you want to write log files 
UNIPROT_GAF=/path/to/gaf #Path to gaf file from uniprot (name is goa_uniprot_all.gaf, retrieved from ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz)
##############################


# I was not able to identify GO annotation online. So I am going to create mine.

# 0. get fasta sequences of genic regions (we don't have them)
cd ${BASEDIR}/10_annotation/gene_fasta
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_cds_from_genomic.fna.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_protein.faa.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_rna_from_genomic.fna.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_translated_cds.faa.gz

# I decided to use the rna from genomics file (the largest, so I assume it has more features) to annotate.
# I checked that it has the same name convention of the gtf (obviously!!)
# 1. blastx against uniprot
cd ${LOGDIR}
BLASTDB=${BLASTDIR}/latest/swissprot
MYQ=${BASEDIR}/10_annotation/gene_fasta/GCA_000005005.6_B73_RefGen_v4_rna_from_genomic.fna
MYB=${BASEDIR}/10_annotation/genes_GO/rna_from_genomic.out
MYGO=${UNIPROT_GAF}
blast1=`echo "module load aligners/blast/latest; blastx -db $BLASTDB -query $MYQ -num_threads 12\
 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -evalue 1e-10 -out $MYB" | \
qsub -N blast1 -l vmem=20G,walltime=168:00:00,nodes=1:ppn=12`

##########################################################################
#Extract the Swissprot ID from blast output and write to file
##########################################################################
BSWISS=${BASEDIR}/10_annotation/genes_GO/blasted_swissprot_ID.txt
cut -f2 $MYB | cut -d"|" -f4| cut -d"." -f1 | sort | uniq > $BSWISS

##########################################################################
#Read the gene_id in $BSWISS and return them if they have a correspondence in the second column of $MYGO
##########################################################################

MYMAPPING=${BASEDIR}/10_annotation/genes_GO/transcripts_GO.txt
#This might be very slow, and we might rather run the command below on the split files.
awk 'FNR==NR{k[$1]=1;next;} k[$2]' $BSWISS $MYGO > $MYMAPPING 


#Reformat GO file, so that each gene has only one row and all GO terms are separated by ";"
FORMGO=${BASEDIR}/10_annotation/formatted_transcripts_GO.txt
Rscript ${FUNCDIR}/03_reformat_GO.r \
-I $MYMAPPING \
-O $FORMGO

#Extract the correct names from fasta files and save them as table. We will use it later 
NAMETOTAG=${BASEDIR}/10_annotation/name_to_tag_converter.txt
grep ">" $MYQ | cut -d" " -f1,2 | sed -e 's/\[locus_tag=//g' | sed -e 's/\]//g' | sed -e 's/>//g' > $NAMETOTAG



#Merge GO terms (in $MYMAPPING) to the DE file and produce a summary table aggregating data by phylum
for COMP in N_vs_TA N_vs_TC N_vs_TE TA_vs_TC TA_vs_TE TC_vs_TE
do
RESP=${BASEDIR}/08_DE_cuffdiff/${COMP}/gene_exp.diff
DEGOFILE=${BASEDIR}/10_annotation/GO_results/${COMP}_DE_GO.txt
Rscript ${FUNCDIR}/04_merge_blast_GO_DE.r \
-G $FORMGO \
-B $MYB \
-N $NAMETOTAG \
-D $RESP \
-O $DEGOFILE
#Run GO enrichment test using topGO
for ONTO in BP CC MF
do
for ALGORITHM in weight01 classic elim
do
OUTFILE=${BASEDIR}/10_annotation/${COMP}_GO_enrich_${ALGORITHM}_${ONTO}.txt
Rscript ${FUNCDIR}/06_GO_enrichment.r -G $FORMGO -I $DEGOFILE -A $ALGORITHM -N $ONTO -O $OUTFILE
done
done
done

#ADD KEGG TERMS (takes a lot of times and might be useless... Still I keep the output)
#It is possible that I never used this function in the paper, I think I used the strategy of 07_biomaRt_convert
for COMP in N_vs_TA N_vs_TC N_vs_TE TA_vs_TC TA_vs_TE TC_vs_TE
do
INFILE=${BASEDIR}/10_annotation/GO_results/${COMP}_DE_GO.txt
OUTFILE=${BASEDIR}/10_annotation/GO_results/${COMP}_DE_GO_KEGG.txt
if [ -e $OUTFILE ]
then
echo "File" $OUTFILE "exists and I will skip this operation" 
else
echo "File" $OUTFILE "doesn't exist and I will run on it"
Rscript ${FUNCDIR}/05_assign_kegg.r \
-I $INFILE \
-O $OUTFILE
fi
done



#ADD ENTREZ TERMS (will be needed to use enrichKEGG in clusterprofiler)
for COMP in N_vs_TA N_vs_TC N_vs_TE TA_vs_TC TA_vs_TE TC_vs_TE
do
INFILE=${BASEDIR}/10_annotation/GO_results/${COMP}_DE_GO.txt
OUTFILE=${BASEDIR}/10_annotation/GO_results/${COMP}_DE_GO_entrez.txt
if [ -e $OUTFILE ]
then
echo "File" $OUTFILE "exists and I will skip this operation" 
else
echo "File" $OUTFILE "doesn't exist and I will run on it"
Rscript ${FUNCDIR}/07_biomaRt_convert.r \
-I $INFILE \
-O $OUTFILE
fi
done


#PERFORM KEGG ENRICHMENT
for COMP in N_vs_TA N_vs_TC N_vs_TE TA_vs_TC TA_vs_TE TC_vs_TE
do
INFILE=${BASEDIR}/10_annotation/GO_results/${COMP}_DE_GO_entrez.txt
OUTFILE=${BASEDIR}/10_annotation/KEGG_results/${COMP}_KEGG_enrich.txt
if [ -e $OUTFILE ]
then
echo "File" $OUTFILE "exists and I will skip this operation" 
else
echo "File" $OUTFILE "doesn't exist and I will run on it"
Rscript ${FUNCDIR}/08_KEGG_enrichment.r \
-I $INFILE \
-O $OUTFILE
fi
done






####################################à
#OLD STUFF BELOW!!!!
####################################à



module load sw/bio/cufflinks/2.2.0
GTF=/genomes/zea_mays/genes/assembly_build_v4/v4.37/annotation/Zea_mays.AGPv4.37.gtf

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


#Write only the output for the selected genes
for aaa in N_vs_TA N_vs_TC N_vs_TE TA_vs_TC TA_vs_TE TC_vs_TE
do
Rscript ${FUNCDIR}/02_cuffdiff.r -I ${BASEDIR}/08_DE_cuffdiff/${aaa}/gene_exp.diff -O ${BASEDIR}/08_DE_cuffdiff/${aaa}/${aaa}_selected.txt
done
