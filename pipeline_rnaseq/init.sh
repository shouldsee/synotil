mkdir -p $JARLIB
#trimmomatic 

#### must supply absolute paths
TRIMDIR=/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32
printf '#!/usr/bin/env bash \njava -jar $TRIMDIR/trimmomatic-0.32.jar "$@"'>$ENVDIR/bin/trimmomatic
mkdir adapters
ln -s $TRIMDIR/adapters/* adapters
cd adapters
cat TruSeq3-PE*.fa >TruSeq3-PE-all.fa
cd ..


# pip install --user pyfaidx
# fastqc is assumed installed

### HISAT2
ln -s /home/Program_NGS_sl-pw-srv01/hisat2-2.1.0/* $ENVDIR/bin

### picard for dedup
ln -s /home/Program_NGS_sl-pw-srv01/picard-tools-1.103/*.jar $JARLIB

###
ln -s /home/Program_NGS_sl-pw-srv01/stringtie-1.3.3b.Linux_x86_64/stringtie $ENVDIR/bin


##### Index genome.fa
mkdir -p ref
#ln -s /home/ref_genew/Brachypodium_Bd21_v3.1/assembly/Bdistachyon_314_v3.0.fa ref/genome.fa
ln -s /home/ref_genew/Brachypodium_Bd21_v3.1/* ref
ln -s $PWD/ref/assembly/*.fa ref/genome.fa
faidx ref/genome.fa -i chromsizes > ref/genome.sizes
