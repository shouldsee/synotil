######################################
#### Hand-coded environment variables
#### $GFF,$GTF,$IDX_HISAT,$GSIZE should have been defined before proceeds

####### Adapter FASTA
export FA_ADAPTER="/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters/TruSeq3-SE.fa"
\
export REF=/home/feng/ref/Brachypodium_Bd21_v3-1

###### Genome annotation .gtf and .gff3 (optional)
export GTF=$(echo $REF/annotation/*.gtf)
export GFF=$(echo $REF/annotation/*.gene_exons.gff3)
export GSIZE="$REF/genome.sizes"
A=$(ls -1 $REF/sequence/Bowtie2Index/* | head -1)
export IDX_BOWTIE2=${A%%.*}
echo A=$A

checkVars GTF GFF GSIZE FA_ADAPTER REF IDX_BOWTIE2

#### Hand-coded environment variables
######################################