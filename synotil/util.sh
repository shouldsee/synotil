#!/bin/bash
# assert.sh Source: http://tldp.org/LDP/abs/html/debugging.html
export UTIL=`readlink -f ${BASH_SOURCE}`
# export UTIL=$SELF

echo "[SOURCE]ing [UTILITIES] from ${UTIL}"

#######################################################################
assert ()                 #  If condition false,
{                         #+ exit from script
                          #+ with appropriate error message.
  E_PARAM_ERR=98
  E_ASSERT_FAILED=99


  if [ -z "$2" ]          #  Not enough parameters passed
  then                    #+ to assert() function.
    return $E_PARAM_ERR   #  No damage done.
  fi

  lineno=$2

  if [ ! $1 ] 
  then
#     echo "Assertion failed:  \"$1\"" 
    echo "Assertion failed:  \"$1\".MSG:\"$3\"" 
    echo "File \"$0\", line $lineno"    # Give name of file and line number.
    exit $E_ASSERT_FAILED
  # else
  #   return
  #   and continue executing the script.
  fi  
} # Insert a similar assert() function into a script you need to debug.    
#######################################################################
export -f assert
# assert "0 -eq 1" $LINENO "hi"

headqc ()
{
  FNAME=$1
  HEADLEN=${2:-100k}
  TMP=$(basename ${FNAME}).${HEADLEN}.tmp
  head -n $HEADLEN $FNAME > $TMP
  fastqc $TMP -o .
  rm $TMP
}
export -f headqc

routine_fastqc ()
{
  INDIR=${1:-$PWD}
  OUTDIR="$PWD/fastqc"
  mkdir -p $OUTDIR
  cd $OUTDIR

  ##### !!!WARN: This for-loop is dangerous and need to be replaced with PARALLEL
  for F in $INDIR/*.fastq
  do
    headqc $F 100k &>>run.log &
    fastqc $F -o . &>>run.log &
  done
}
export -f routine_fastqc

array_to_vars() {
  # Source: https://stackoverflow.com/a/45242139/8083313
  declare -n _arr=$1
  local var
  for var; do
    shift || return
    printf -v "$var" %s "$1"
  done
}
export -f array_to_vars

check_PE() {
  read1=$1
  read2=$2
  #echo read1:
  #ls -lh $read1; guessPhred $read1
  #echo read2:
  #ls -lh $read2; guessPhred $read2
  ali1=$(basename $read1)
  ali1=${ali1%.*}
  ali2=$(basename $read2)
  ali2=${ali2%.*}
  echo $ali1
  echo $ali2
  assert $(guessPhred $read1)=$(guessPhred $read2) $LINENO
  PHRED=$(guessPhred $read1)
  PHRED=${PHRED/+/''}
  PHRED=${PHRED,,}
  echo $PHRED
}
export -f check_PE

bamqc() {
  ### QC a sorted, dedup bam
  BAM=$1
  ALI=${BAM%.bam}
  echo $ALI
  samtools index $BAM $ALI.bai &
  pids[0]=$!
  samtools flagstat $BAM >$ALI.flagstat.log &  
  pids[1]=$!
  for pid in ${pids[*]}; do
    wait $pid
  done  
}
export -f bamqc

bam2bigwig() {
    local BAM=$1
    local GSIZE=$2
    local NORM=${3:-RPGC}
    ALI=$(bname $BAM)
    
    #### -split argument is essential !!!
    genomeCoverageBed -ibam $BAM -bg -split > $ALI.bdg 
    bedGraphToBigWig $ALI.bdg $GSIZE $ALI.bw 
    bamCoverage --normalizeUsing $NORM --effectiveGenomeSize `size2sum $GSIZE` \
      --smoothLength 10 --binSize 10 -p ${NCORE:-1} \
      -b $BAM -o ${ALI}_${NORM}.bw
}
export -f bam2bigwig

bam2normBW(){
    local BAM=$1
    local NORM=${2:-RPGC}
}
export -f bam2normBW

size2sum(){
    local GSIZE=$1
    cut -d$'\t' -f2 $GSIZE | paste -sd+ | bc
}
export -f size2sum

routine_fastqc() 
{ 
    F=$1
    headqc $F 100k &>> run.log &
    pid[0]=$!
    fastqc $F -o . &>> run.log &
    pid[1]=$!
    for p in ${pid[*]}
    do
        wait $p
    done
}
export -f routine_fastqc

bname() {
  basename ${1%%.*}
}
export -f bname

function datefloat() {
  date +%s.%N
}
export -f datefloat


headbam() {
    BAM=$1
    NLINE=${2:-100k}
    samtools view -h $BAM | head -n${NLINE} | samtools view -bS
}
export -f headbam

guessPhred() {
    inputfile=$1
    # Source: http://onetipperday.sterding.com/2012/10/code-snip-to-decide-phred-encoding-of.html
    # Update: The following part for checking the file extension can be simplified (Thanks to the comment from Unknown) 

    less $inputfile | head -n40 | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding!";}'
}
export -f guessPhred

headfastq() {
  FILE=$1  
  HEADLEN=${2:-100k}
  ALI=$(bname $FILE)
  if [[ $FILE == *.fastq.gz ]]; then
      gzip -dc $FILE | head -n $HEADLEN | gzip -c >$ALI.fastq.gz
  elif [[ $FILE == *.fastq ]]; then
      head -n $HEADLEN $FILE >$ALI.fastq
  fi
}
export -f headfastq

envPull() {
    chmod +x $ORIGIN/* 
    cp -f $ORIGIN/* $ENVDIR/bin
#     chmod +x $ENVDIR/bin/*.sh 
}
export -f envPull


dusort(){
 du -csh "$@" | sort -h 
}
export -f dusort

flatten(){
    #### Link the input directory to $PWD in a flattened manner
    BNAME=`basename $1`
    local ODIR=${2:-${BNAME}_flat}
    mkdir -p $ODIR
    find $1 -mindepth 1 -type f -exec ln -f -t $ODIR '{}' +
}
export -f flatten

flattenDir(){
    INDIR=`readlink -f $1`
    BNAME=`basename $1`
    mkdir -p ${BNAME}_flat
    cd ${BNAME}_flat
    for D in $INDIR/*
    do
        flatten $D
    done
    cd ..
}
export -f flattenDir


regGroup(){
    ##### Regroup files in the current directory into subfolders
    ##### depending on the output of a REGEX

    local REG=${1:-".*(_S[0-9]+_).*"}
    local ARR=`ls -1 | sed -E --expression="s/$REG/\1/g" |  uniq`
    for ID in ${ARR[@]}
    do
        DI=${ID//_/}
#         echo $ID $DI
        INF=`ls -1 | grep $ID`
        mkdir -p $DI; mv $INF $DI 
    done
}
export -f regGroup

tidyFastq(){
    #### Tidy-up a Project folder containing multiple runs
    INDIR=${1%/}
    ODIR=${2:-${INDIR}_tidy}
    echo $INDIR;    echo "[Output]ing to $ODIR"
    mkdir -p $ODIR
    flatten $INDIR $ODIR
    cd $ODIR; regGroup ".*(_S[0-9]+_).*"    
    cd ..
}
export -f tidyFastq

flattenFastq(){
    #### Take a .fastq directory from illumina download and flatten it.
    #### First level directory names are concatentaed as prefix
    INDIR=${1%/}
    local ODIR=${2:-${INDIR}_tidy}
    echo $INDIR;    echo "[Output]ing to $ODIR"
    mkdir -p $ODIR
    for D in $(ls -d1 $INDIR/*)
    do         
        echo $D, ${ODIR}/${D#*/}
        flatten $D ${ODIR}/${D#*/}
    done
    cd $ODIR
    flatten_keepname *
    cd ..
}
export -f flattenFastq

findLeaf(){
    ### Find leaf files in the input directories
    find "$@" -mindepth 1 -type f 
}
export -f findLeaf

flatten_keepname(){
    ### Find leaf files in the input directories
    set -e
    local DIR=(`find "$@" -maxdepth 0 -type d`)
    echo ${DIR[@]}
    for F in `find "${DIR[@]}" -mindepth 1 -type f`
    do
    mv -f $PWD/$F ${F//\//_}
    done    
    rm -rf "${DIR[@]}" 
}
export -f flatten_keepname

routine_indexGenome ()
{
#     F=$(echo *.fa);
    IN=${1}
    faidx -i chromsizes $IN > ${IN%.*}.sizes
}
export -f routine_indexGenome

checkVars ()
{
    for F in "$@";
    do
        FILE=${!F};
#         echo [FILE]=$FILE
#         [[ ! -z "$FILE" ]] || { echo $F variable not set ; exit 255 ; } 
        if [[ -z "$FILE" ]]; then echo $F variable not set; exit 255 ; fi
        echo [Test] $F=$FILE;
#         [[ -f "$FILE" ]] || [[ -d "$FILE" ]] || ls -l ${FILE}* || { echo "$F=$FILE does not exists!" ;  exit 255 ; }
        if [[ -f "$FILE" ]]; then : 
        else
            if [[ -d "$FILE" ]]; then
                :
            else 
                ls -l ${FILE}* || { echo "$F=$FILE does not exists!" ;  exit 255 ; }
            fi 
        fi
    done
}
export -f checkVars

tabCut(){
    cut -f"$2" -d$'\t' $1
}
export -f tabCut
addFunc () 
{ 
    echo
    local func=$1;
    type $func | tail -n +2;
    echo export -f $func
    echo
}
export -f addFunc
erun () 
{ 
    local CMD="$@";
    SHELOG=${SHELOG:-log.sh};
    echo "$CMD" >> $SHELOG && eval "$CMD"
}
export -f erun
bamHist () 
{ 
    NCORE=${NCORE:-4};
    BAM=$1;
    [[ -e ${1}.bai ]] || { echo $1 is not indexed, now indexing ; samtools index $1 ; } 
    ALI=`bname $BAM`;
    CMD="bamPEFragmentSize -b $BAM -o ${ALI}_bamHist.png -p $NCORE >${ALI}_bamHist.log";
    erun $CMD
}
export -f bamHist


GTF2CDSR(){
    IN=$1
    OFILE=${2:-$(basename $IN).cds}
    gtf2bed< $IN  | sed "s/\"//g" | grep CDS >tmp
    cat tmp | sort -k1,1 -k4,4 -k5,5 -k6,6  | bedtools groupby -g 1,4,5,6 -c 2,3,8 -o min,max,first \
    | awk -v OFS='\t' '{print $1, $5, $6, $2, $3, $4, $7, $8}' \
    > $OFILE
    echo "[GTF2CDSR] saved to $OFILE"
}
export -f GTF2CDSR

# |  awk -v OFS='\t' '{print $1, $2, $3, $4, $5}'>tmp.cds
# cat tmp | bedtools groupby -g 1,4,5,6 -c 2,3,8 -o min,max,first >tmp.cds
# head *.cds
# gtf2CDSR
autoinstallPython () 
{ 
    local BASE=$1;
    cd $BASE && ls -1 */*.py | entr ./setup.py install --user
}
export -f autoinstallPython


sortAll()
{
IN=$1
tr x '\000' < $IN | sort | tr '\000' x
}
export -f sortAll

fastq2seq () 
{ 
    IN=$1;
    paste - - - - < $IN | cut -f2 > $(basename $IN).raw
}
export -f fastq2seq

uniqCount () 
{ 
    IN=$1;
    OUT=$(basename $IN).count;
    wc -l $IN > $OUT;
    cat $IN | sort | uniq -c | sort -nr >> $OUT
}
export -f uniqCount
quickFasta()
{
    local IN=$1
    local FI=${2:-$FA_GENOME}
    checkVars IN FI
    local OUT=`basename ${IN%.*}`.fa
    bedtools getfasta -name+ -s -fi $FI -bed $IN -fo $OUT 
}
export -f quickFasta

wrapper_meme ()
{ 
    IN=$2;
    PROG=$1;
    ARG="${@:3}";
    ODIR="${ARG// /_}";
    ODIR="${ODIR//-/=}";
    CMD="$PROG $ARG -oc $ODIR $IN *.fa";
    echo $CMD;
    $CMD 2>&1 | tee ${ODIR}/${PROG}.log
}
export -f wrapper_meme
cmd2dir(){
    local PROG=$1
    ARG="${@:2}"
    ADIR="${ARG// /_}"
    ADIR="${ADIR//-/=}";
    ODIR="PROG=${PROG}_${ADIR}";
    echo $ODIR
}
export -f cmd2dir

routine_mast()
{
    local PROG=mast
    local IN1=$1
    local IN2=$2    
    local ALI1=`basename ${IN1%.*}`
    local ALI2=`basename ${IN2%.*}`
    local OPT="${@:3}"
    local OPT="-w -remcor $OPT"
    local ODIR=`cmd2dir $PROG $OPT $ALI1 $ALI2`; mkdir -p $ODIR
#     echo []PROG=$PROG
#     echo []OPT=$OPT
#     echo []ODIR=$ODIR
#     type cmd2dir
    local CMD="$PROG -oc $ODIR $OPT $IN1 $IN2"
    echo $CMD
    [[ $DRY -eq 1 ]] || $CMD 2>&1 | tee $ODIR/${PROG}.log
}   
export -f routine_mast

# routine_tomtom()
# {
#     local PROG=mast
#     local IN1=$1
#     local IN2=$2    
#     local ALI1=`basename ${IN1%.*}`
#     local ALI2=`basename ${IN2%.*}`
#     local OPT="${@:3}"
#     local OPT="-w -remcor $OPT"
#     local ODIR=`cmd2dir $PROG $OPT $ALI1 $ALI2`; mkdir -p $ODIR
# #     echo []PROG=$PROG
# #     echo []OPT=$OPT
# #     echo []ODIR=$ODIR
# #     type cmd2dir
#     local CMD="$PROG -oc $ODIR $OPT $IN1 $IN2"
#     echo $CMD
#     [[ $DRY -eq 1 ]] || $CMD 2>&1 | tee $ODIR/${PROG}.log
# }   
# export -f routine_tomtom


cpLink()
{
    cp -l "$@" || cp "$@"
}  
export -f cpLink
runWithTimeLog()
{
    local CMD="$@"
    local ALI=${ALI:-testALI}
    local PROG=${PROG:-testPROG}
    local T0 T1 Tdiff
    local SELF=${SELF:-testScript.sh}
    
#     echo $CMD
    [[ $DRY -eq 1 ]] || {
        T0=`datefloat`
        $CMD 2>&1 | tee ${ALI}.${PROG}.log 
        T1=`datefloat`
        Tdiff=`echo $T1 - $T0 | bc`
        echo $SELF,$Tdiff,\"$CMD\"
    }
}
export -f runWithTimeLog
    
routine_fimo () 
{ 
    local IN=$1;
    local ODIR="PROG\=fimo_IN=`basename ${IN%.*}`";
    local OPT="${@:2}";
    local PROG=fimo;
    local CMD="$PROG --oc $ODIR  $OPT $IN";
    echo $CMD;
    [[ $DRY -eq 1 ]] || $CMD
}
export -f routine_fimo

routine_findSummit () 
{ 
    local IN=$1;
    local ALI=`basename ${IN%.*}`;
    checkVars IN GSIZE;
    bedtools slop -i $IN -pct -r -0.5 -l -0.5 -g $GSIZE | tee ${ALI}_summit\=yes.bed | head
}
export -f routine_findSummit

quickTargz () 
{ 
    local IN=`readlink -f $1`;
    local ALI=`basename $IN`;
    rm -f ${ALI}.tar.gz;
    ls -1 "$IN" > .${ALI}.tmp;
    tar -cvzf ${ALI}.tar.gz -C $IN `cat .${ALI}.tmp`
}
export -f quickTargz
    
trimN () 
{ 
    local IN=$1;
    local N=${2:-'4-'};
    local CMD="cut -c${N} $IN";
    echo $CMD;
    [[ $DRY -eq 1 ]] || { 
        $CMD > "${IN}.tmp" && mv "${IN}.tmp" "$IN";
        echo "[Trimmed] $IN"
    }
}
export -f trimN
bed_sortMerge () 
{ 
    local BED=$1;
    local merge_arg="${@:2:-test}";
    [[ $merge_arg -eq "" ]] && { 
        merge_arg="-c 4 -o mean"
    };
    local CMD="bedtools sort -i $BED -g $GSIZE > ${BED}.sorted ;
    bedtools merge -i ${BED}.sorted $merge_arg";
    (>&2 echo $CMD);
    [[ $DRY -eq 1 ]] || eval $CMD
}
export -f bed_sortMerge


prompt_yn () 
{ 
    local MSG=$1;
    while true; do
        read -p "$MSG: (y/n)" yn;
        case $yn in 
            [Yy]*)
                break
            ;;
            [Nn]*)
                exit
            ;;
            *)
                echo "Please answer yes(y) or no(n)"
            ;;
        esac;
    done
}
export -f prompt_yn

