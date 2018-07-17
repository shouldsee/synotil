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
export FA_GENOME="$REF/genome.fa"

A=$(ls -1 $REF/sequence/Bowtie2Index/* | head -1)
export IDX_BOWTIE2=${A%%.*}

###### HISAT2 index
A=$(ls -1 $REF/sequence/HISAT2Index/* | head -1)
export IDX_HISAT2=${A%%.*}

export DB_MOTIF=" \
/home/feng/ref/motif/CIS-BP/Brachypodium_distachyon.meme \
/home/feng/ref/motif/ARABD/ArabidopsisDAPv1.meme \
/home/feng/ref/motif/ARABD/ArabidopsisPBM_20140210.meme"

# checkVars GTF GFF GSIZE FA_ADAPTER REF IDX_BOWTIE2
checkVars GTF GSIZE FA_ADAPTER_SE FA_ADAPTER_PE REF IDX_BOWTIE2 IDX_HISAT2 \
    FA_GENOME DB_MOTIF
    
#### Hand-coded environment variables
######################################