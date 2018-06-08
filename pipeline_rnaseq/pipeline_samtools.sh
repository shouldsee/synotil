#!/usr/bin/env bash
main(){
    local SELF
    SELF=`readlink -f ${BASH_SOURCE[0]}`
    SELFALI=$(bname $SELF)

    set -e ###exit on error

    #ADADIR="/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters"

    PROG=samtools
    read1=$(readlink -f $1)
    NCORE=${2:-4}

    #### Prepare Directory
    OUTDIR=$PWD/$PROG
    mkdir -p $OUTDIR
    cd $OUTDIR

    T0=`datefloat`
    MSG=".sam to sorted .bam"
    echo "===== Starting $MSG ====="
    echo "CALLED: $@"

    # echo "===== Ending $MSG  ====="


    ##### Checking input format
    echo $read1
    read ali1 ali1 PHRED <<< $(check_PE $read1 $read1)
    echo $ali1
    echo $ali2
    echo "Using $NCORE threads "
    # echo "(sorting with 10 threads)"
    echo Phred quality version: $PHRED
    echo "(Expect unkown encoding for a .sam input)"

    #### Shared alias
    ali=${ali1%_R1_*}
    echo $ali
    ALI=$ali

    OPT="--threads $NCORE"
    CMD="$PROG view $read1 -b $OPT -o ${ali}.bam"
    echo $CMD
    time `$CMD &> ${ali}.runlog`

    OPT="--threads $NCORE -m 4G"
    CMD="$PROG sort ${ali}.bam $OPT -o ${ali}.sorted.bam"
    echo $CMD

    time `$CMD &>> ${ali}.runlog `

    ln -f $PWD/${ali}.sorted.bam ../${ali}.sorted.bam
    cd ..

    echo "===== Ending $MSG  ====="
    DUR=$(echo $(datefloat) - $T0 | bc)
    echo ${BASH_SOURCE[0]},$DUR >> $ALI.time 
    # echo ${BASH_SOURCE[0]},$DUR | tee $ALI.time 
}
main "$@"