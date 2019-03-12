#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys,os
import pymisca.ext as pyext
def main(fname,ofname=None):
    if ofname is None:
        ofname = os.path.basename(fname)+'.count'
#     if os.path.exists(fname):
        
    dfc = pyext.readData(fname,columns=['gene_id','unstranded','forward','reverse'])
    BEST = dfc.head(4).sum(axis=0).sort_values().index[0]
    # dfcc=  dfc[BEST].tail(-4)
    # dfcc = d
    dfcc = dfc[BEST].to_frame('rawCount')
    dfcc['CPM'] = dfcc['rawCount'] / dfcc['rawCount'].sum() * 10**6
    dfcc.index.name = 'gene_id'
    dfcc.to_csv(ofname,sep='\t')
#     print ('[CMD]:%s/FILE.json'%os.path.dirname(ofname)
    pyext.util__fileDict.main(
        ofname='%s/FILE.json'%os.path.dirname(ofname),
        argD= dict(STAR_BESTCOL=BEST)
    )
    return ofname

if __name__ == '__main__':
    assert len(sys.argv) >=2
    fname = sys.argv[1]
    if len(sys.argv) == 3:
        ofname = sys.argv[2]
    else:
        ofname = None
    main(fname,ofname)

