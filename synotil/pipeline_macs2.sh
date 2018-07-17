#!/bin/bash
main()
{
    local PROG=macs2
    local SELF=${BASH_SOURCE}
    local SELF_ALI=`basename ${SELF%.*}`
    
    local BAM=$1  ### input bam file
#    GCOUNT=${2:-`size2sum $GSIZE`} ### lenght of your genome
    GCOUNT=`size2sum $GSIZE` ### lenght of your genome
    OPT="${@:2}"
    ALI=`basename ${BAM%.*}`
    CMD="macs2 callpeak $OPT -t $BAM --keep-dup 1 -n $ALI -g $GCOUNT -p 0.1"    
    runWithTimeLog "$CMD" | tail -1 >> ${ALI}.time
}
main "$@"
