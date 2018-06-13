######################################
#### Hand-coded environment variables
#### $GFF,$GTF,$IDX_HISAT,$GSIZE should have been defined before proceeds

####### Adapter FASTA
ADADIR="/home/feng/ref/adapters"
export FA_ADAPTER_SE="$ADADIR/TruSeq3-SE.fa"
export FA_ADAPTER_PE="$ADADIR/TruSeq3-PE-all.fa"

export REF=/home/feng/ref/Brachypodium_Bd21_v3-1
###### Genome annotation .gtf and .gff3 (optional)
export GTF=$(echo $REF/annotation/*.gtf)
export GFF=$(echo $REF/annotation/*.gene_exons.gff3)
export GSIZE="$REF/genome.sizes"
A=$(ls -1 $REF/sequence/Bowtie2Index/* | head -1)
export IDX_BOWTIE2=${A%%.*}

# checkVars GTF GFF GSIZE FA_ADAPTER REF IDX_BOWTIE2
checkVars GTF GSIZE FA_ADAPTER_SE FA_ADAPTER_PE REF IDX_BOWTIE2 # GFF IDX_HISAT

#### Hand-coded environment variables
######################################