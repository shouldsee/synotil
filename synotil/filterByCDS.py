#!/usr/bin/env python2

import synotil.util as sutil
import pymisca.util as pyutil
import sys,os

import argparse

parser= argparse.ArgumentParser()
parser.add_argument('-v','--verbose',default=0,type = int)

parser.add_argument('peakFile',)
parser.add_argument('-r','--peakRadius',default=1, type=int,)

parser.add_argument('-c','--cdsFile',default=pyutil.os.environ.get('GTF','none')+'.cds')
parser.add_argument('-g','--GSIZE',default=pyutil.os.environ.get('GSIZE',None))
parser.add_argument('-u','--upStream',default =1000, type=int, )
parser.add_argument('-d','--downStream',default=500, type=int,)

parser.add_argument('-f','--featFile',default=None)
parser.add_argument('-s','--centerSummit',default=0)

# stdout = sys.stdout
# stderr = sys.stderr

def main(peakFile, 
         peakRadius = parser.get_default('peakRadius'),
         cdsFile=parser.get_default('cdsFile'),
         upStream=parser.get_default('upStream'),
         downStream=parser.get_default('downStream'),
         GSIZE=parser.get_default('GSIZE'),
         featFile = parser.get_default('featFile'),
         verbose = parser.get_default('verbose'),
         centerSummit = parser.get_default('centerSummit'),
         center_summit = None
        ):
    if center_summit is None:
        center_summit = centerSummit
    else:
        pass
       
    if featFile is None:
        featFile = sutil.findPromoter( cdsFile,
                                   upStream=upStream, 
                                   downStream=downStream, 
                                   GSIZE=GSIZE,
                                   silent= ~verbose)
        msg = '[MSG]making feature file from cdsFile:\n %s\>\>\n%s \n'%(cdsFile,featFile)
        pyutil.sys.stderr.write(msg)
    '''

Filter a .narrowPeak file to discard peaks too far from a start codon.

Example:
    ./npk/186C/S12/1505-Zt12-27C_S12_peaks.narrowPeak 
    -g /home/feng/ref/Arabidopsis_thaliana_TAIR10/genome.sizes 
    -c /home/feng/ref/Arabidopsis_thaliana_TAIR10/annotation/genes.gtf.cds
'''
#     ! head -2 {qfile}; echo
#     ! head -2 {exon_file}

    peakTemp = sutil.npk_expandSummit(fname=peakFile,
                                      radius=peakRadius,
                                     center_summit=center_summit)
    LC = pyutil.lineCount(peakTemp)
    msg = '[tmp]:{peakTemp},lc={LC}\n'.format(**locals())
    sys.stderr.write(msg)
    
#     print ('[tmp]:',peakTemp,pyutil.lineCount(peakTemp))
    ofname = sutil.closestAnnotation(
        peakTemp,
        ANNOTATION_FILE=featFile,
        RANGE=0,
        GSIZE=GSIZE,silent= 1 - verbose)
    
    LC = pyutil.lineCount(ofname)
    msg = '[output]:{ofname},lc={LC}\n'.format(**locals())
    sys.stderr.write(msg)

    sys.stdout.write(ofname+'\n')
    
#     print ('[output]',ofname,pyutil.lineCount(ofname))
    return ofname

if __name__ =='__main__':
    args = parser.parse_args()

    main(**vars(args))