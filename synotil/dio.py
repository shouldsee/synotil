import CountMatrix as scount
import pymisca.util as pyutil
def readChipPeak(fname,**kwargs):
    '''CountMatrix.geneList is buggy at the moment hence needed to be converted back into 
    CountMatrix.countMatrix
'''
    df = scount.countMatrix.from_DataFrame(fname=fname,ext='tsv',index_col='geneAcc',addFname=0,**kwargs)
    df= df.toGeneList()
    return df
#     df.columns = [pyutil.basename(fname)]
#     chipPeak = df
#     return df