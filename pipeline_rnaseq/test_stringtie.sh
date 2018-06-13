# CMD="htseq-count -f bam $BAM $GTF >$ALI.htseq.gtf"
# time `eval $CMD` &>>log 

NCORE=8
CMD="stringtie -p $NCORE --rf $BAM -G $GTF -o test.gtf"
{
    echo $CMD
    time `$CMD`
} &>>log

