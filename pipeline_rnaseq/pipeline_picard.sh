#!/usr/bin/env bash
main(){

    local SELF
    SELF=`readlink -f ${BASH_SOURCE[0]}`
    SELFALI=$(bname $SELF)

    set -e ###exit on error


    PROG=picard

    read1=$(readlink -f $1)
    #NCORE=4 ### Hard-coded, maybe try benchmarking 
    NCORE=${2:-4}


    #### Prepare Directory
    OUTDIR=$PWD/$PROG
    mkdir -p $OUTDIR
    cd $OUTDIR

    T0=$(datefloat)
    MSG="Deduplication with Picard"
    echo "===== Starting $MSG ====="
    # echo "===== Ending $MSG  ====="


    ##### Checking input format
    echo $read1
    ali1=$(basename $read1)
    ali1=${ali1%%.*}

    #### Shared alias
    ali=${ali1%_R1_*}
    echo $ali
    ALI=$ali

    CMD="java -XX:ParallelGCThreads=$NCORE -jar"
    CMD="$CMD $JARLIB/MarkDuplicates.jar"
    CMD="$CMD I=$read1 O=${ali}_dedup.bam M=${ali}.dupstat.log"
    CMD="$CMD REMOVE_DUPLICATES=true "
    echo $CMD
    # CMD="$PROG view $read1 -b $OPT"
    time `$CMD &> ${ali}.$PROG.log`
    ln -f $PWD/${ali}_dedup.bam ../${ali}.bam
    ln -f $PWD/*.log ../
    cd ..

    echo "===== Ending $MSG  ====="
    DUR=$(echo $(datefloat) - $T0 | bc)
    echo ${BASH_SOURCE[0]},$DUR >>$ALI.time 
}
main "$@"