#!/usr/bin/env bash
main(){
    local SELF
    SELF=`readlink -f ${BASH_SOURCE[0]}`
    SELFALI=$(bname $SELF)

    # set -e
    ### Kind of weird here...
    source $(dirname $SELF)/activate

    ######################################
    #### Hand-coded environment variables
    #### $GFF,$GTF,$IDX_HISAT,$GSIZE should have been defined before proceeds

#     ####### Adapter FASTA
#     export ADADIR="/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters"
#     export FA_ADAPTER="$ENVDIR/adapters/TruSeq3-SE.fa"

#     ###### Genome annotation .gtf and .gff3 (optional)
#     export GTF=$(echo "$ENVDIR/ref/annotation/*.gtf")
#     export GFF=$(echo "$ENVDIR/ref/annotation/*.gene_exons.gff3")
#     export GSIZE="${ENVDIR}/ref/genome.sizes"

#     ###### BOWTIE2 index
#     A=$(ls -1 $ENVDIR/ref/Bowtie2Index/* | head -1)
#     export IDX_BOWTIE2=${A%%.*}
    #### Hand-coded environment variables
    ######################################
    
    ######################################
    echo ==== Parse Input arguments
    {
        read1=$1 ####  e.g. test_R1_.fastq
    #     read2=$2 ####  e.g. test_R2_.fastq
        NCORE=${2:-6} #### number of threads, default to 6
        DIR=$PWD
        ALI1=$(bname $read1)
    #     ALI2=$(bname $read2)
        ALI=${ALI1%_R1_*}
        LOGFILE=${ALI}.${SELFALI}.log
        echo []Proposed alias ${ALI} ==
        echo []Logfile $LOGFILE
        rm -f $LOGFILE
    }
    {
    
    checkVars GTF GSIZE FA_ADAPTER REF IDX_BOWTIE2 # GFF 
    echo "===== IMPORTANT VARS =====
### Call
read1=$1
SELF=$SELF
ALI=$ALI   ### alias
LOGFILE=$LOGFILE

### Genome
FA_ADAPTER=$FA_ADAPTER
GSIZE=$GSIZE
IDX_BOWTIE2=$IDX_BOWTIE2

### Not Used
GTF=$GTF
GFF=$GFF

" | tee -a $LOGFILE
    } 
    
    
assert ()
{
    E_PARAM_ERR=98;
    E_ASSERT_FAILED=99;
    if [ -z "$2" ]; then
        return $E_PARAM_ERR;
    fi;
    lineno=$2;
    if [ ! $1 ]; then
        echo "Assertion failed:  \"$1\".MSG:\"$3\"";
        echo "File \"$0\", line $lineno";
#         exit $E_ASSERT_FAILED;
    fi
}



    ######################################
    echo ==== main program
    {
        T00=`datefloat`

        echo 'command,time(s)'>$ALI.time
        #######################

        ########### Starting pipeline
        mkdir -p fastqc; cd fastqc; headqc ../$read1 100k; cd ..

        pipeline_trim_se.sh $read1 $NCORE 
        assert "$? -eq 0" $LINENO "Trimmomatic/fastqc failed"

    #     pipeline_hisat.sh  $ALI1.fastq  $ALI2.fastq $IDX_HISAT $NCORE
        pipeline_bowtie2_se.sh  $ALI1.fastq $IDX_BOWTIE2 $NCORE
        assert "$? -eq 0" $LINENO "BOWTIE2 failed"

        pipeline_samtools.sh ${ALI}.sam $NCORE
        assert "$? -eq 0" $LINENO "SAM2BAM/SORT failed"

        pipeline_picard.sh ${ALI}.sorted.bam $NCORE
        assert "$? -eq 0" $LINENO "Picard Deduplication failed"

        pipeline_bamqc.sh ${ALI}.bam $GSIZE 
        assert "$? -eq 0" $LINENO "BAMQC/conversion failed"

        cleanup() {
            set -e
            ALI=$1
            ln -f *.log *.time output/
    #         ln -f *.count *.gtf output/
            ln -f *.bw *.bdg output/
            ln -f $ALI.bam $ALI.bai output/
            ln -f fastqc/*.html output/
        }

        mkdir -p output
        touch output/DESCRIPTION
        cleanup ${ALI} &>${ALI}.cleanup.log
        assert "$? -eq 0" $LINENO "Output/cleanup failed"

        DUR=$(echo $(datefloat) - $T00 | bc)
        echo $SELF,$DUR >>$ALI.time 

        echo "[$SELFALI]:OUTPUT has been deposited into $PWD/output"
        echo "[$SELFALI]: ...Ending..."
        echo ---- Main program
    }  &> $LOGFILE
}
echo [MAIN]main "$@"
main "$@"