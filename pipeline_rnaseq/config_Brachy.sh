######################################
#### Hand-coded environment variables
#### $GFF,$GTF,$IDX_HISAT,$GSIZE should have been defined before proceeds

####### Adapter FASTA
export ADADIR="/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters"
export FA_ADAPTER="$ENVDIR/adapters/TruSeq3-SE.fa"

###### Genome annotation .gtf and .gff3 (optional)
export GTF=$(echo "$ENVDIR/ref/annotation/*.gtf")
export GFF=$(echo "$ENVDIR/ref/annotation/*.gene_exons.gff3")
export GSIZE="${ENVDIR}/ref/genome.sizes"

###### BOWTIE2 index
A=$(ls -1 $ENVDIR/ref/Bowtie2Index/* | head -1)
export IDX_BOWTIE2=${A%%.*}
#### Hand-coded environment variables
######################################