INDIR="${ENVDIR}/out"
INPUT="$INDIR/Exp0024-ZT8-Bdphyc_S2.bam"
GSIZE=${ENVDIR}/ref/genome.sizes
echo "==INPUT DIR:$INDIR"

pipeline_bamqc.sh $INPUT $GSIZE 
