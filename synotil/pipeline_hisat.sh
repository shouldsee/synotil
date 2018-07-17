#!/usr/bin/env bash
main(){

    local SELF
    SELF=`readlink -f ${BASH_SOURCE[0]}`
    SELFALI=$(bname $SELF)

    set -e ###exit on error

    #ADADIR="/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters"

    PROG=hisat2


    T0=$(datefloat)
    MSG="spliced alignment with HISAT2"
    echo "===== Starting $MSG ====="
    # echo "===== Ending $MSG  ====="

    read1=$1
    read2=$2
    GIDX=$3
    NCORE=${4:-4}

    ##### Checking input format
    echo $read1
    read ali1 ali2 PHRED <<< $(check_PE $read1 $read2)
    echo $ali1
    echo $ali2
    echo Using $NCORE threads

    ####!! overiding
    #PHRED=phred33
    echo Phred quality version: $PHRED

    #### Shared alias
    ali=${ali1%_R1_*}
    echo $ali
    ALI=$(echo $ali)



    OPT="--threads $NCORE --no-mixed --rna-strandness RF --dta --fr"
    CMD="$PROG -x $GIDX -1 $read1 -2 $read2 -S ${ali}.sam $OPT" 

    echo $CMD
    time `$CMD &> ${ali}.$PROG.log`

    #ln -f $PWD/${ali}.sam ../${ali}.sam
    #cd ..

    echo "===== Ending $MSG  ====="
    DUR=$(echo $(datefloat) - $T0 | bc)
    echo ${BASH_SOURCE[0]},$DUR >>$ALI.time

}
main "$@"
#$CMD test.out
