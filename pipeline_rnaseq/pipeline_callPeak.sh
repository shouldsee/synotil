#!/bin/bash
# %%bash
main(){
    set -e
    local SELF=`readlink -f ${BASH_SOURCE[0]}`        
    local SELFALI=`basename ${SELF%.*}`
    
    show_help(){
    echo " No help for you yet :X
Usage:
    pipeline_callPeak.sh [OPTIONS] [-c control.bam] <IN.bam> 
Options:
    -d dry-run
    -r [RADIUS]  peaks within this radius will be reported
    -q [Q_cutoff] peaks with Q-value lower than the cutoff will be kept
    
Required environment variables:
    GSIZE: *.sizes, size of each chromosome
    GFF: GFF annotation 
    " 
        
    checkVars GSIZE GFF
        
    }
    if [[ $# -eq 0 ]] ; then
        show_help
        exit 0
    fi
    
    ###############
    #### ==== Parsing argument
    local CTRL=""
    local GSUM=2.7E9
    local PAD=1000
    local QLIM=5E-2
    local TARGET=""
         
    OPTIND=1         # Reset in case getopts has been used previously in the shell.
    while getopts "h?c:j:g:r:q:dt" opt; do
        case "$opt" in
        h|\?)
            show_help
            exit 0
            ;;
        c)  CTRL=$OPTARG
            ;;
        j)  NCORE=$OPTARG
            ;;
        g)  GSUM=$OPTARG
            ;;
        r)  PAD=$OPTARG
            ;;
        q)  QLIM=$OPTARG
            ;;
        d)  local DRY=1
            ;;
        t)  TARGET=1
            ;;
        esac
    done
    shift $((OPTIND-1))
    [ "${1:-}" = "--" ] && shift
    #### ---- Parsing argument
    ###############
    echo "[Positional]: $@"
    ## dbg
#     checkVars QLIM 
    IN=$1 ; 
    if [[ -d $IN ]]; then
        IN=`find $IN -type f -wholename "*.bam" | head -n1 `
    else
         echo [NOT A DIR?]
    fi         
    IN=`readlink -f $IN`
    ALI=`basename ${IN%.*}`     
    echo [ALI=$ALI]
    checkVars IN   
         
    #### Reference annotation
    grep gene $GFF | bedtools sort  >tmp.ref

    {
        echo "#### Step 1: MACS2 callpeak "
        PROG="macs2 callpeak"
        PROGALI=${PROG// /_}
        OPT="-q $QLIM -g $GSUM --keep-dup 1"

        ### Use a control if specified
        if [[ ! -z "$CTRL" ]]; then
            OPT="$OPT -c $CTRL"
            ALI="${ALI}_ctrl:`basename ${CTRL%.*}`"
        else
            ALI=${ALI}_ctrl:none
        fi

        CMD="$PROG $OPT -t $IN -n $ALI &>$ALI.$PROGALI.log"
        echo $CMD ; [[ ! -z $DRY ]] || eval $CMD    

        echo "#### Step 2:  Inner join annotation (GFF) to called peaks"
        INBED=${ALI}_peaks.narrowPeak       
        ### tabCount() used since The output of 'bedtools closest' is not dircetly accepted by bedtools ;
        CMD="
        bedtools slop -i tmp.ref -g $GSIZE -b $PAD >tmp.ref.s${PAD} ;
        bedtools closest -a tmp.ref.s${PAD} -b $INBED >tmp.out   ;

        tabCut tmp.out -9  | bedtools slop -g $GSIZE -b -$PAD -i - >tmp.out.-9 ; 
        tabCut tmp.out 10- | paste -d$'\t' tmp.out.-9 - >tmp.bed ;
        ln -f tmp.bed ${ALI}.bed ;
        "
        echo $CMD; [[ ! -z $DRY ]] || eval $CMD   

        echo "#### Step 3: Select the peak with highest"
        CMD="uniqMaxCol.py ${ALI}.bed ${ALI}_uniq.bed  2>${ALI}.uniq.log;" 
        echo $CMD; [[ ! -z $DRY ]] || eval $CMD   


        echo "#### Step 4: Output to dir"
        {
            OUTDIR="script:${SELFALI}_sample:${ALI}"
            CMD="mkdir -p $OUTDIR ; 
            cp *.log *.bed -t $OUTDIR"
            echo $CMD ; [[ ! -z $DRY ]] || eval $CMD
            
            REMDIR=`dirname $IN`                   
            CMD="cp -R `readlink -f $OUTDIR` -t $REMDIR"  
            echo "$CMD" | tee push.sh
            if [[ ! -z "$TARGET" ]]; then
                [[ ! -z $DRY ]] || . push.sh
            fi        
        }
        echo "[SUCC][DRY=$DRY] MACS2 pipeline finished for bam input:$IN"
    }
    set +e
}

main "$@"
         