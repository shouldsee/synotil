#!/bin/bash
# assert.sh Source: http://tldp.org/LDP/abs/html/debugging.html

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
    BAM=$1
    GSIZE=$2
    ALI=$(bname $BAM)
    
    #### -split argument is essential !!!
    genomeCoverageBed -ibam $BAM -bg -split > $ALI.bdg 
    bedGraphToBigWig $ALI.bdg $GSIZE $ALI.bw
}
export -f bam2bigwig

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