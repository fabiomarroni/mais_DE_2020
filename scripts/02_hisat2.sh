#!/bin/bash 


#Variables you need to set
BASEDIR=/path/to/your/results #Exiting folder that you want to be the top folder of your results
FUNCDIR=/path/to/R/functions #Folder in which you will store the R function provided
SOFTWAREDIR=/path/to/software #Path to software you need (if it is in your path, you don't need to specify it)
##############################


module load sw/aligners/hisat2/2.0.4
#Build index (only once)

REF=/genomes/zea_mays/assembly/reference/v4/Zea_mays.AGPv4.dna.toplevel.fa
INDEXDIR=${BASEDIR}/05_Indices
mkdir -p $INDEXDIR
hisat2-build $REF $INDEXDIR/Zea_mays.AGPv4.dna.toplevel.fa -p 16

READDIR=${BASEDIR}/00_Runs_collector
ALDIR=${BASEDIR}/06_hisat2
GTF=/genomes/zea_mays/genes/assembly_build_v4/v4.37/annotation/Zea_mays.AGPv4.37.gtf
NCORES=8
mkdir -p $ALDIR

#Align to reference
for INFILE in ${READDIR}/*R1*fastq.gz
do
echo $INFILE
R1=$INFILE
OUTFILE=${INFILE/$READDIR/$ALDIR}
OUTFILE=${OUTFILE/_L001_R1_001.fastq.gz/.sam}
echo $OUTFILE
echo "module load sw/aligners/hisat2/2.0.4;hisat2 $INDEXDIR/Zea_mays.AGPv4.dna.toplevel.fa --no-unal -U $R1 -p $NCORES -S ${OUTFILE}.sam" | qsub -N h2_align -l vmem=32G,walltime=24:00:00,nodes=1:ppn=$NCORES
done

#Convert sam to bam
module load it/tools/samtools/1.7
for INFILE in ${ALDIR}/*sam
do
echo "module load it/tools/samtools/1.7; samtools view -b $INFILE --threads 7 -o ${INFILE/.sam/.bam}" | qsub -N stv -l vmem=16G,walltime=12:00:00,nodes=1:ppn=8
done

#Sort bam files
for INFILE in ${ALDIR}/*bam
do
echo "module load it/tools/samtools/1.7; samtools sort $INFILE -o ${INFILE/.bam/_sort.bam}" | qsub -N stv -l vmem=16G,walltime=12:00:00,nodes=1:ppn=8
done


#Remove sam files
rm ${ALDIR}/*sam

