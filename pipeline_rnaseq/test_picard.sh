#./pipeline_rnaseq.sh 


INDIR="${ENVDIR}/out"
INPUT="$INDIR/Exp0024-ZT8-Bdphyc_S2.sorted.bam"

echo 
echo "==INPUT DIR:$INDIR"
##### Simply 'head' a .fastq file not accepted by trimmomatic
#INPUT=" $INDIR/testR1.fastq $INDIR/testR2.fastq"
#INPUT=" $INDIR/t1.fq $INDIR/t1.fq"

pipeline_picard.sh $INPUT 10


#routine_fastqc .
#head $ENVDIR/out/1.out -n100k | fastqc stdin -o $ENVDIR/out/.

