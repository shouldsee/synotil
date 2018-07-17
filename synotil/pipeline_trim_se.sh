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
    echo "===== Trimmomatic ====="
    echo $1
    echo $2

    read1=`readlink -f $1`
    NCORE=${2:-4}

    read ALI1 ALI1 PHRED <<< $(check_PE $read1 $read1)
    ALI=${ALI1%_R1_*}

    echo Using $NCORE threads
    echo Phred quality version: $PHRED
    echo $FA_ADAPTER


    #### Paired-end routine
    OUTDIR=$PWD/trimmed
    mkdir -p $OUTDIR
    cd $OUTDIR


    SETTING="ILLUMINACLIP:${FA_ADAPTER_SE}:6:30:10"
    SETTING="$SETTING LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15"
    CMD="trimmomatic SE -threads $NCORE -$PHRED $read1 ${ALI1}_pass.fastq $SETTING" 

    echo $CMD
    time `$CMD &> ${ALI}.runlog`

    ##### Link files for later
    ln -f $PWD/${ALI1}_pass.fastq ../${ALI1}.fastq
    # ln -f $PWD/${ALI2}_pass.fastq ../${ALI2}.fastq
    cd ..

    mkdir -p fastqc 
    cd fastqc
    routine_fastqc ../trimmed/${ALI1}_pass.fastq &
    routine_fastqc ../trimmed/${ALI1}_fail.fastq &
    # routine_fastqc ../trimmed/${ALI2}_pass.fastq &
    # routine_fastqc ../trimmed/${ALI2}_fail.fastq &
    cd ..

    # T0=$(datefloat)
    DUR=$(echo $(datefloat) - $T0 | bc)
    echo ${BASH_SOURCE[0]},$DUR >>$ALI.time 
    #$CMD test.out
}
main "$@"