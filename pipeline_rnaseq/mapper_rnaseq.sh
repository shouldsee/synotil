#!/bin/bash
set -e
SELF=${BASH_SOURCE[0]}

# INDIR=/media/pw_synology3/PW_HiSeq_data/RNA-seq/Raw_data/testONLY
# OUTDIR=/media/pw_synology3/BrachyPhoton/Mapped_data
# ID=133R/BdPIFs-32747730/133E_23_DN-40235206


# INDIR=/media/pw_synology3/BrachyPhoton/raw

INDIR=$HOME/envs/bulkrna/raw
OUTDIR=/media/pw_synology3/BrachyPhoton/Mapped_data
PROG=pipeline_rnaseq.sh

{
    #### echo ==== Parsing arugments
    ID=$1
    # e.g.: ID=./150R/Doro_150R_Doro_1-43239982
    ID=${ID#./}
    echo $ID; echo PROGRAM:$PROG;
    
}

{
    #### echo ==== Downloading fastq files
    mkdir -p $OUTDIR
    TEMPDIR=`preprocessor.py $INDIR/$ID | tee -a bulk.log | tail -1 `
}

{
    #### echo ==== Running pipeline
    cd $TEMPDIR
    read1=(*_R1_raw.fastq)
    read2=(*_R2_raw.fastq)
    $PROG *_R1_raw.fastq *_R2_raw.fastq
}

{
    #### echo ==== Uploading outputs
    mkdir -p $OUTDIR/$ID
    cp output/* $OUTDIR/$ID
}  
cd ..

echo '[FINISH]: Outputed to $OUTDIR/$ID'
