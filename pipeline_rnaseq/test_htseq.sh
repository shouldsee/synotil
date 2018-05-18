# CMD="htseq-count -f bam $BAM $GTF >$ALI.htseq.gtf"
# time `eval $CMD` &>>log 



INPUT=$ENVDIR/out/Exp0024-ZT8-Bdphyc_S2.bam
BAM=$INPUT
ALI=${BAM%.bam}

CMD="htseq-count -s reverse -f bam $BAM $GTF  -o $ALI.htseq.sam >$ALI.htseq.count 2>htseq.log"
{
    echo $CMD
    time `eval $CMD` 
}&>>log 
