#!/usr/bin/env bash
main(){
    local SELF
    SELF=`readlink -f ${BASH_SOURCE[0]}`
    SELFALI=$(bname $SELF)

    # set -e
    ### Kind of weird here...
    source $(dirname $SELF)/activate

    
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
        rm -f $LOGFILE; touch $LOGFILE
    }
    
    {
    echo "===== IMPORTANT VARS ====="
    checkVars read1 SELF LOGFILE 
    echo "===== Genome Vars ====="
    checkVars GSIZE FA_ADAPTER_SE REF IDX_BOWTIE2 
    #checkVars GTF  GFF 
    } | tee -a $LOGFILE

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

        pipeline_bamqc.sh ${ALI}.bam $GSIZE $NCORE
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