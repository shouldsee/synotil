import CountMatrix as scount
import pymisca.util as pyutil

bedHeader = '''
0:chrom
1:start
2:end
3:acc
4:score
5:strand
6:FC
7:neglogPval
8:neglogQval
9:summit
'''.strip().splitlines()
bedHeader = [x.split(':')[1] for x in bedHeader] 


def readChipPeak(fname,**kwargs):
    '''CountMatrix.geneList is buggy at the moment hence needed to be converted back into 
    CountMatrix.countMatrix
'''
    df = scount.countMatrix.from_DataFrame(fname=fname,ext='tsv',index_col='geneAcc',addFname=0,**kwargs)
    df= df.toGeneList()
    return df


def extract_peak(fname,ext='tsv',header=None,guess_index=0,**kwargs):
    df = pyutil.readData(fname,ext=ext,header=header,guess_index=guess_index,**kwargs)
    df.columns = bedHeader + list(df.columns)[len(bedHeader):]
    return df
    

def npk_expandSummit(df=None,radius=200,fname = None,clip = 1):
    '''
    Expand the summit regions of a .narrowPeak dataFrame
'''
    if df is None:
        df = extract_peak(fname)
    df['abs_summit'] = df.start + df.summit
    st = df.strand
    df.start = (df.abs_summit - radius)
    df.end  = df.abs_summit + radius 
    df.drop('abs_summit',1,inplace=True)
    if clip:
        df.start = df.start.clip_lower(0)
    if fname is not None:
        base = pyutil.basename(fname)
        ofname = '%s_radius=%d.tsv' % (base,radius)
        df.to_csv(ofname,sep='\t',index=None,header=None)
        return ofname
    else:
        return df

#     df.columns = [pyutil.basename(fname)]
#     chipPeak = df
#     return df