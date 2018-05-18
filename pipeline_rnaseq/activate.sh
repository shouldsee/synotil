#!/usr/bin/env bash
#!/bin/bash

called=$_
[[ $called != $0 ]] && echo "Script is being sourced" || echo "Script is being run"
echo "\$BASH_SOURCE ${BASH_SOURCE[0]}"
# echo "\$BASH_SOURCE ${BASH_SOURCE[*]}"
DIR=${BASH_SOURCE[0]%/*}

echo $DIR


DIR=$(readlink -f $DIR)

source ${DIR}/util.sh


export ENVDIR=${DIR%/*}
export JARLIB=$ENVDIR/jar
mkdir -p $JARLIB
echo Adding $ENVDIR to PATH
export PATH="$PATH:$DIR:$ENVDIR"


export ADADIR="/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters"
export FA_ADAPTER="$ENVDIR/adapters/TruSeq3-PE-all.fa"
export GTF="$ENVDIR/ref/annotation/*.gtf"
export GTF=$(echo $GTF)
export GFF="$ENVDIR/ref/annotation/*.gene_exons.gff3"
export GFF=$(echo $GFF)



A=$(ls -1 $ENVDIR/ref/HISAT2Index/* | head -1)
export IDX_HISAT=${A%%.*}

