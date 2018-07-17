#!/usr/bin/env bash
main(){
    local SELF
    SELF=`readlink -f ${BASH_SOURCE[0]}`
    SELFALI=$(bname $SELF)

    set -e ###exit on error

    #### Take a bam file, index and flagstat, followed by conversion to .bdg and .bw
    #### check util.sh::bamqc()  util.sh::bam2bigwig()

    INPUT=$1
    GSIZE=$2
    NCORE=${3:-1}

    T0=`datefloat`
    MSG="QC .dedup.bam and convert to .bdg .bw"
    echo "===== Starting $MSG ====="


    ALI=$(bname $INPUT)
    CMD="bamqc $INPUT &>> ${ALI}.runlog"
    echo $CMD
    time `eval $CMD` &  pid[0]=$!

    CMD="bam2bigwig $INPUT $GSIZE &>> ${ALI}.runlog"
    time `eval $CMD` &  pid[1]=$!

    for p in ${pid[*]}
    do
        wait $p
    done

    echo "===== Ending $MSG  ====="
    DUR=$(echo $(datefloat) - $T0 | bc)
    echo ${BASH_SOURCE[0]},$DUR >>$ALI.time 
}
main "$@"