set -e ###exit on error

#### Take a bam file, index and flagstat, followed by conversion to .bdg and .bw
#### check util.sh::bamqc()  util.sh::bam2bigwig()

INPUT=$1
GSIZE=$2

T0=`datefloat`
MSG="QC .dedup.bam and convert to .bdg .bw"
echo "===== Starting $MSG ====="


ALI=$(bname $INPUT)
time `bamqc $INPUT &>> ${ALI}.runlog` & 
pid[0]=$!
time `bam2bigwig $INPUT $GSIZE &>> ${ALI}.runlog` & 
pid[1]=$!

for p in ${pid[*]}
do
    wait $p
done

echo "===== Ending $MSG  ====="
DUR=$(echo $(datefloat) - $T0 | bc)
echo ${BASH_SOURCE[0]},$DUR >>$ALI.time 
