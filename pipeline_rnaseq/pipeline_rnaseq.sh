#!/usr/bin/env bash
SELF=${BASH_SOURCE[0]}

##### Input arguments
read1=$1 ####  e.g. test_R1_.fastq
read2=$2 ####  e.g. test_R2_.fastq
NCORE=${3:-6} #### number of threads, default to 6


DIR=$PWD
ALI1=$(bname $read1)
ALI2=$(bname $read2)
ALI=${ALI1%_R1_*}


#####  $GFF,$GTF,$FA_HISAT,$GSIZE should have been defined before proceeds
##### Checking $ENVDIR/bin/activate
GSIZE=${ENVDIR}/ref/genome.sizes

echo $ALI1
echo $ALI2
echo $ALI

echo $GTF
echo $IDX_HISAT



T00=`datefloat`
echo '== Proposed alias ${ALI} =='
echo 'command,time(s)'>$ALI.time

#INPUT="$INDIR/Exp0024-ZT8-Bdphyc_S2_R1_raw.fastq $INDIR/Exp0024-ZT8-Bdphyc_S2_R2_raw.fastq"

##### Simply 'head' a .fastq file not accepted by trimmomatic
#INPUT=" $INDIR/testR1.fastq $INDIR/testR2.fastq"
#INPUT=" $INDIR/t1.fq $INDIR/t1.fq"
#######################



########### Starting pipeline

pipeline_trim.sh $read1 $read2 $NCORE 
assert "$? -eq 0" $LINENO "Trimmomatic/fastqc failed"

pipeline_hisat.sh  $ALI1.fastq  $ALI2.fastq  $IDX_HISAT $NCORE
assert "$? -eq 0" $LINENO "HISAT2 failed"

pipeline_samtools.sh ${ALI}.sam $NCORE
assert "$? -eq 0" $LINENO "SAM2BAM/SORT failed"

pipeline_picard.sh ${ALI}.sorted.bam $NCORE
assert "$? -eq 0" $LINENO "Picard Deduplication failed"

pipeline_bamqc.sh ${ALI}.bam $GSIZE 
assert "$? -eq 0" $LINENO "BAMQC/conversion failed"

CMD="stringtie -p $NCORE --rf ${ALI}.bam -G $GTF -o ${ALI}.stringtie.gtf -A ${ALI}.stringtie.count &> ${ALI}.stringtie.log"

{   
    T0=`datefloat`
    echo $CMD
    time `eval $CMD`
    DUR=$(echo $(datefloat) - $T0 | bc)
    ARR=($CMD); echo $(which ${ARR[0]}),$DUR >>$ALI.time 
}
assert "$? -eq 0" $LINENO "Stringtie failed"


CMD="htseq-count -s reverse -f bam ${ALI}.bam $GTF -r pos -o $ALI.htseq.sam >$ALI.htseq.count 2>${ALI}.htseq.log"

##### trying out -s parameter
#CMD="htseq-count -s yes -f bam ${ALI}.bam $GTF -r pos -o $ALI.htseq.sam >$ALI.htseq.count 2>${ALI}.htseq.log"
{   
    T0=`datefloat`
    echo $CMD
    time `eval $CMD`
    DUR=$(echo $(datefloat) - $T0 | bc)
    ARR=($CMD); echo $(which ${ARR[0]}),$DUR >>$ALI.time 
}
assert "$? -eq 0" $LINENO "HTSeq-count failed"

DUR=$(echo $(datefloat) - $T00 | bc)
echo $SELF,$DUR >>$ALI.time 
