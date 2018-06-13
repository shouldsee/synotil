#./pipeline_rnaseq.sh 


INDIR="${ENVDIR}/out"
GIDX="/home/ref_genew/Brachypodium_Bd21_v3.1/HISAT2Index/Bdistachyon314_Bd"
INPUT="$INDIR/Exp0024-ZT8-Bdphyc_S2_R1_raw.fastq $INDIR/Exp0024-ZT8-Bdphyc_S2_R2_raw.fastq $GIDX"

echo 
echo "==INPUT DIR:$INDIR"
##### Simply 'head' a .fastq file not accepted by trimmomatic
#INPUT=" $INDIR/testR1.fastq $INDIR/testR2.fastq"
#INPUT=" $INDIR/t1.fq $INDIR/t1.fq"


# pipeline_hisat.sh $INPUT 10
pipeline_hisat_test.sh $INPUT 10

#routine_fastqc .
#head $ENVDIR/out/1.out -n100k | fastqc stdin -o $ENVDIR/out/.

