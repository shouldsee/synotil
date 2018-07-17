#!/usr/bin/env bash
main()
{
    local SELF
    SELF=`readlink -f ${BASH_SOURCE[0]}`
    SELFALI=$(bname $SELF)
    set -e ###exit on error
    local T0=$(datefloat)

    #### $FA_ADAPTER must be set

    #ADADIR="/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters"
    #ADADIR="$ENVDIR/adapters"
    echo 
    echo "===== Trimmomatic ====="

    read1=`readlink -f $1`
    read2=`readlink -f $2`
    NCORE=${3:-4}

    read ALI1 ALI2 PHRED <<< $(check_PE $read1 $read2)
    echo $ALI1
    echo $ALI2
    echo Using $NCORE threads
    ALI=${ALI1%_R1_*}

    ####!! overiding
    #PHRED=phred33
    ####
    echo Phred quality version: $PHRED
    echo $FA_ADAPTER


    #### Paired-end routine
    OUTDIR=$PWD/trimmed
    mkdir -p $OUTDIR
    cd $OUTDIR


    SETTING="ILLUMINACLIP:${FA_ADAPTER_PE}:6:30:10"
    SETTING="$SETTING LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15"
    CMD="trimmomatic PE -threads $NCORE -$PHRED $read1 $read2 ${ALI1}_pass.fastq ${ALI1}_fail.fastq ${ALI2}_pass.fastq ${ALI2}_fail.fastq $SETTING" 

    echo $CMD
    time `$CMD &> ${ALI}.runlog`

    ##### Link files for later
    ln -f $PWD/${ALI1}_pass.fastq ../${ALI1}.fastq
    ln -f $PWD/${ALI2}_pass.fastq ../${ALI2}.fastq
    cd ..

    mkdir -p fastqc 
    cd fastqc
    routine_fastqc ../trimmed/${ALI1}_pass.fastq &
    routine_fastqc ../trimmed/${ALI2}_pass.fastq &
    routine_fastqc ../trimmed/${ALI1}_fail.fastq &
    routine_fastqc ../trimmed/${ALI2}_fail.fastq &
    cd ..

    # T0=$(datefloat)
    DUR=$(echo $(datefloat) - $T0 | bc)
    echo ${BASH_SOURCE[0]},$DUR >>$ALI.time 
    #$CMD test.out
}
main "$@"