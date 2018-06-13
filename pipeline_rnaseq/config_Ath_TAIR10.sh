######################################
#### Hand-coded environment variables
#### $GFF,$GTF,$IDX_HISAT,$GSIZE should have been defined before proceeds

####### Adapter FASTA
ADADIR="/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters"
export FA_ADAPTER_SE="$ADADIR/TruSeq3-SE.fa"
export FA_ADAPTER_PE="$ADADIR/TruSeq3-PE-2.fa"

export REF=/home/feng/ref/Arabidopsis_thaliana_TAIR10
###### Genome annotation .gtf and .gff3 (optional)
export GTF=$(echo $REF/annotation/*.gtf)
# export GFF=$(echo $REF/annotation/*.gene_exons.gff3)
export GSIZE="$REF/genome.sizes"
A=$(ls -1 $REF/sequence/Bowtie2Index/* | head -1)
export IDX_BOWTIE2=${A%%.*}

checkVars GTF GSIZE FA_ADAPTER REF IDX_BOWTIE2 # GFF

#### Hand-coded environment variables
######################################