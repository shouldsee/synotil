#!/bin/bash
### pipeline_ame.sh <in.fasta> <db.meme>
routine_mast()
{
    PROG=mast
    IN1=$1
    IN2=$2    
    ALI1=`basename ${IN1%.*}`
    ALI2=`basename ${IN2%.*}`
    OPT="${@:3}"
    OPT="-w -remcor $OPT"
    ODIR=`cmd2dir $PROG $OPT $ALI1 $ALI2`; mkdir -p $ODIR
#     echo []PROG=$PROG
#     echo []OPT=$OPT
#     echo []ODIR=$ODIR
#     type cmd2dir
    CMD="$PROG -oc $ODIR $OPT $IN1 $IN2"
    echo $CMD
    [[ $DRY -eq 1 ]] || $CMD 2>&1 | tee $ODIR/${PROG}.log
}   

main()
{
    set -e
    local PROG=ame
    local IN=$1
    local DB_MEME=${@:2}
    DB_MEME=${DB_MEME:-$DB_MOTIF}
    checkVars IN DB_MEME

    IN=`readlink -f $IN`
    ODIR="PROG=${PROG}_`bname $IN`_`basename -a $DB_MEME | tr '\n' '_' | tr '.' '-'`"
    
    mkdir -p $ODIR; cd $ODIR; 
    cp -l $IN . || cp -s $IN . || cp $IN .

    OPT="--verbose 2 --kmer 2 --control --shuffle--"
    OPT="$OPT --hit-lo-fraction 0.25 --evalue-report-threshold 10.0"
    CMD="$PROG $OPT --oc . $IN $DB_MEME"
    
    echo $CMD
    [[ $DRY -eq 1 ]] || $CMD 2>&1 | tee $PROG.log
    set +e

    ame2meme.sh ame.tsv 0.05 > MAST-INPUT
    routine_mast MAST-INPUT $IN
    cd $OLDPWD
}
main "$@"