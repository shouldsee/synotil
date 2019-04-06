
# import pymisca.util as pyutil
import re
Rcomp=re.compile

ridL = r'(^|/)(\d{1,4}[RCQ]{1,2}|SRR\d{7,8}'
# ridL = r'(\d{1,4}[RC]'
end = r'[_\/\-]'
runID = Rcomp('.*%s).*'%(ridL,))
runCond = Rcomp('%s%s.*)'%(ridL,end))
ridPATH = Rcomp('%s%s.*)'%(ridL,end))
runID_RNA = Rcomp('.*%s).*'%(ridL.replace('C',''),))

srr = Rcomp('(?P<lead>SRR.*)_(?P<read>[012]).(?P<ext>.+)')
baseSpace = Rcomp('(?P<lead>.*)_L(?P<chunk>\d+)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)')
baseSpaceSimple = Rcomp('(?P<lead>.*)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)')
sampleID = Rcomp('[_\/](S\d{1,3}|SRR\d{7,8})[_/\-\.]')

BdAcc = Rcomp('(Bradi[\da-zA-Z]+)')

key__fastqRaw = [
    'runID',
    'sampleID',
    'read',
    'chunk']

def get_runID(fname):
#     pt = r'[\^/](\d{3}R)/'
    pt = runID
    res = re.findall(pt,fname)    
    return res[-1][-1] if res else None
###legacy
getRUNid = get_runID

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
        'sampleID':get_sampleID,
    }.items()}
    return res
def add_ZTime_int(meta):
    vals = [ int(re.sub('[^-\d]','',x)) for x in meta['ZTime'] ]
    vals = [24+x if x <0 else x for x in vals]
    meta['ZTime_int'] = vals
#     meta['ZTime'] = pyutil.paste0([ ['ZT'] * len(meta),meta['ZTime_int']])
    return meta


    
    