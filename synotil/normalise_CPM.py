#!/usr/bin/env python2
import pandas as pd
import sys, os
#import subprocess
#import pymisca.util as pyutil


### Remove extra rows
#CMD = "grep -v '$_' {INFILE} > {INFILE}.tmp".format(**locals())
#subprocess.check_output(CMD,shell=1)

#df = pd.read_csv( INFILE +'.tmp',header=None,sep='\t',index_col=None)

#CMD = "rm {INFILE}.tmp".format(**locals())
#subprocess.check_output(CMD,shell=1)

def main(INFILE, ofname = None):
    ####
    df = pd.read_csv( INFILE ,header=None,sep='\t',index_col=None)

    #df = pyutil.readData(INFILE,header=None)

    df.columns = ['gene_id','read_count']
    ## drop extra columns
    df = df.loc[~df.gene_id.str.startswith('_')]
    df['CPM'] = df.read_count.astype('float')/df.read_count.sum() * 1E6

    
    s = (df.to_csv(ofname, sep='\t',index=False))
    if ofname is not None:
        return ofname
    else:
        return s
    
    
    
if __name__ ==  '__main__':
    assert len(sys.argv)>=2,'not enough arguments'
    INFILE = sys.argv[1]
    s = main(INFILE)
    sys.stdout.write(s)
