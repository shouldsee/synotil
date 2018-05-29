#!/usr/bin/env bash
set -e
##### Initialising a working environment in current directory 
SELF=${BASH_SOURCE[0]}
ENVDIR=$PWD

##### install binaries
JARLIB=$PWD/jar
mkdir -p $PWD/bin
mkdir -p $JARLIB

# pip install --user pyfaidx
# fastqc is assumed installed

### Trimmomatic
# must supply absolute paths
TRIMDIR=/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32
printf '#!/usr/bin/env bash \njava -jar $TRIMDIR/trimmomatic-0.32.jar "$@"'>$ENVDIR/bin/trimmomatic

##### Prepare adapter seqs
mkdir -p $PWD/adapters
ln -sf $TRIMDIR/adapters/* adapters
cd adapters
if [ ! -e "TruSeq3-PE-all.fa" ]; then
#     :
# else
    cat TruSeq3-PE*.fa >TruSeq3-PE-all.fa
fi
cd ..


### HISAT2
ln -sf /home/Program_NGS_sl-pw-srv01/hisat2-2.1.0/* $ENVDIR/bin
### picard for dedup
ln -sf /home/Program_NGS_sl-pw-srv01/picard-tools-1.103/*.jar $JARLIB
### StringTie
ln -sf /home/Program_NGS_sl-pw-srv01/stringtie-1.3.3b.Linux_x86_64/stringtie $ENVDIR/bin


##### Index genome.fa
mkdir -p $PWD/ref
#ln -s /home/ref_genew/Brachypodium_Bd21_v3.1/assembly/Bdistachyon_314_v3.0.fa ref/genome.fa
ln -sf /home/ref_genew/Brachypodium_Bd21_v3.1/* ref
ln -sf $PWD/ref/assembly/*.fa ref/genome.fa
faidx ref/genome.fa -i chromsizes > ref/genome.sizes


cp -f ${SELF%/*}/activate.sh bin/activate
cp -f ${SELF%/*}/*.sh bin/
