#!/usr/bin/env bash
local SELF
SELF=`readlink -f ${BASH_SOURCE[0]}`
SELFALI=$(bname $SELF)

### Kind of weird here...
source $(dirname $SELF)/activate

######################################
#### Hand-coded environment variables
#### $GFF,$GTF,$IDX_HISAT,$GSIZE should have been defined before proceeds

####### Adapter FASTA
export ADADIR="/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters"
export FA_ADAPTER="$ENVDIR/adapters/TruSeq3-PE-all.fa"

###### Genome annotation .gtf and .gff3 (optional)
export GTF=$(echo "$ENVDIR/ref/annotation/*.gtf")
export GFF=$(echo "$ENVDIR/ref/annotation/*.gene_exons.gff3")
export GSIZE="${ENVDIR}/ref/genome.sizes"

###### HISAT2 index
A=$(ls -1 $ENVDIR/ref/HISAT2Index/* | head -1)
export IDX_HISAT=${A%%.*}

#### Hand-coded environment variables
######################################


######################################
echo ==== Parse Input arguments
{
    read1=$1 ####  e.g. test_R1_.fastq
    read2=$2 ####  e.g. test_R2_.fastq
    NCORE=${3:-6} #### number of threads, default to 6
    DIR=$PWD
    ALI1=$(bname $read1)
    ALI2=$(bname $read2)
    ALI=${ALI1%_R1_*}
#         echo $ALI1; echo $ALI2; echo $ALI
    LOGFILE=$ALI.pipeline_rnaseq.log
    echo []Proposed alias ${ALI} ==
    echo []Logfile $LOGFILE
}

######################################
echo ==== main program
{
    T00=`datefloat`
    

    
    echo 'command,time(s)'>$ALI.time
    #######################



    ########### Starting pipeline

    pipeline_trim_pe.sh $read1 $read2 $NCORE 
    assert "$? -eq 0" $LINENO "Trimmomatic/fastqc failed"

    pipeline_hisat.sh  $ALI1.fastq  $ALI2.fastq $IDX_HISAT $NCORE
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



    cleanup() {
        set -e
        ALI=$1
        mkdir -p output
        ln -f *.log *.time output/
        ln -f *.count *.gtf output/
        ln -f *.bw *.bdg output/
        ln -f $ALI.bam $ALI.bai output/
        ln -f fastqc/*.html output/
    }

    cleanup ${ALI} &>${ALI}.cleanup.log
    assert "$? -eq 0" $LINENO "Output/cleanup failed"

    DUR=$(echo $(datefloat) - $T00 | bc)
    echo $SELF,$DUR >>$ALI.time 

    echo "[$SELFALI]:OUTPUT has been deposited into $PWD/output"
    echo "[$SELFALI]: ...Ending..."
    echo ---- Main program
}  &> $LOGFILE