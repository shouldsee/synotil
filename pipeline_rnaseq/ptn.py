import re
Rcomp=re.compile
ridL = r'(^|/)(\d{1,4}[RC]'
# ridL = r'(\d{1,4}[RC]'
end = r'[_\/\-]'
runID = Rcomp('.*%s).*'%(ridL,))
runCond = Rcomp('%s%s.*)'%(ridL,end))

baseSpace = Rcomp('(?P<lead>.*)_L(?P<chunk>\d+)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)')
sampleID = Rcomp('_(S\d{1,3})[_/\-\.]')

BdAcc = Rcomp('(Bradi[\da-zA-Z]+)')


import pymisca.util as pyutil
def getRUNid(fname):
#     pt = r'[\^/](\d{3}R)/'
    pt = runID
    res = re.findall(pt,fname)    
    return res[-1][-1] if res else None

def get_sampleID(fname):
    res = sampleID.findall(fname)    
    return res[0] if res else None
def getRUNid_full(fname):
    try:
        pt = r'[\^/](\d{3}R/[^/]+)/'
        res = re.findall(pt,fname)    
        return res[0] if res else None
    except:
#         print fname
        assert 0,fname
def getZT(fname):
    pt = '[/\^-](ZT[-]?\d{1,2})[-_\.]'
    res = re.findall(pt,fname)
    return res[0] if res else None
def getLight(fname):
    pt = '[/\-_]([LS]D)[\-_]'
    res = re.findall(pt,fname)
    return res[0] if res else None
def getGenotype(fname):
    pt = '[/\-_](elf3|ELF3-OX|Bd[^-]+)[\-_]'
    res = re.findall(pt,fname)
    res = res[0] if res else None
    if res is None:
        return res
    else:
        if not res.startswith('Bd'):
            res = 'Bd%s'%res
        res = res.replace('_','')
        res = res.replace('-','')
#         res = res.replace('-','_')
        return res
def getBrachy(fname):

    res= {k:f(fname) for k,f in {
        'RunID':getRUNid,
        'ZTime':getZT,
        'light':getLight,
        'gtype':getGenotype,
        'SampleID':get_sampleID,
    }.items()}
    return res
def add_ZTime_int(meta):
    vals = [ int(re.sub('[^-\d]','',x)) for x in meta['ZTime'] ]
    vals = [24+x if x <0 else x for x in vals]
    meta['ZTime_int'] = vals
#     meta['ZTime'] = pyutil.paste0([ ['ZT'] * len(meta),meta['ZTime_int']])
    return meta

    
    
if __name__=='__main__':
    import tmp
    import subprocess,re,sys
    import pandas as pd
    import pymisca.util as pyutil

#     fnames = 
    Rid2Age={
        '150R':2,
        '148R':2,
        '149R':2,
        '144R':3,
        '143R':4,
        '169R':4,
    }
    
    INDIR = '/media/pw_synology3/BrachyPhoton/Mapped_data'
#     INDIR = [ '%s/%s'%(INDIR,x) for x in Rid2Age.keys()]
    res = subprocess.check_output("find %s -name *.stringtie.count"%INDIR,shell = 1).splitlines()
    assert res[-1]!='[FAIL]'
    fnames = res
    
    meta = {}
    vals = map(getRUNid,fnames)
    meta['RunID']=vals
    print vals

    vals = map(get_sampleID,fnames)
    meta['sampleID']=vals
    print vals
    
    vals = map(getRUNid_full,fnames)
    meta['RunID_full']=vals
    print vals

    vals = map(getZT,fnames)
    meta['ZTime']=vals

    vals = map(getLight,fnames)
    meta['light']=vals

    vals = map(getGenotype,fnames)
    vals = ['Bd21' if x=='BdWT' else x for x in vals]
    meta['gtype'] = vals
    
    


    meta = pd.DataFrame(meta)
    meta.loc[meta['RunID']=='143R','light'] = 'LD' 
    # meta[meta['RunI']]

    meta['Age_int'] = map(Rid2Age.get,meta['RunID'])
    meta['Age'] =[ 'Wk%d'%x for x in meta['Age_int']]

    #### Discussed with Mingjun on June 6th about 143R
#     meta['ZTime'][[4,5]] = 'ZT0'
    meta.loc[ (meta['RunID']=='143R' ) & ( meta['sampleID'].isin(['S5','S6']) ), 'ZTime'] = 'ZT0'
#     meta['ZTime'][[4,5]] = 'ZT0'

    vals = [ int(re.sub('[^-\d]','',x)) for x in meta['ZTime'] ]
    vals = [24+x if x <0 else x for x in vals]
    meta['ZTime_int'] = vals
    meta['ZTime'] = pyutil.paste0([ ['ZT'] * len(meta),meta['ZTime_int']])

    meta['fname'] = list(fnames)


    meta.to_csv('meta.csv')
    meta
#     sys.exit(0)