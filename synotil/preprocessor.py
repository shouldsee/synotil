#!/usr/bin/env python2
# Usage: (python) preprocessor.py /path/to/FASTQ_DIR
# Example: preprocessor.py /media/pw_synology3/PW_HiSeq_data/RNA-seq/Raw_data/testONLY/133R/BdPIFs-32747730/133E_23_DN-40235206
# Purpose: Download fastq files from the supplied path and 
#    combine them into R1.fastq and R2.fastq
# 
# Created:7  OCT 2016, Hui@SLCU map-RNA-seq.py
# Update: 11 Feb 2017. Hui@SLCU
# Update: 29 May 2018. Feng@SLCU preprocessor.py


import tempfile,subprocess
import os, sys, datetime, glob, re
import multiprocessing as mp
import pandas as pd
import ptn


# shellexec = os.system

def check_all_samples(d):
    '''
    [Deprecated] Use LeafFiles(DIR) instead, kept for a reference
    Recursively walking a directory tree
    TBM to return useful output
    '''
    for k in sorted(d.keys()):
        rpath = d[k]
        result = check_directory(rpath)
        if result != '':
            print('%s:%s [%s]' % (k, rpath, result))
            rpath = rpath.rstrip('/')
            if glob.glob(rpath) is []: # rpath does not exist or it is an incomplete path
                path_lst = glob.glob(rpath + '*/')
                assert len(path_lst) == 1,'Found multiple directories under %s'%rpath
                rpath = path_lst[0]
                d[k] = rpath.rstrip('/') # update rpath



def nTuple(lst,n,silent=1):
    """ntuple([0,3,4,10,2,3], 2) => [(0,3), (4,10), (2,3)]
    
    Group a list into consecutive n-tuples. Incomplete tuples are
    discarded e.g.
    
    >>> group(range(10), 3)
    [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    """
    if not silent:
        L = len(lst)
        if L % n != 0:
            print '[WARN] nTuple(): list length %d not of multiples of %d, discarding extra elements'%(L,n)
    return zip(*[lst[i::n] for i in range(n)])

def LinesNotEmpty(sub):
    sub = [ x for x in sub.splitlines() if x]
    return sub

def LeafFiles(DIR):
    ''' Drill down to leaf files of a directory tree if the path is unique.
    '''
    assert os.path.exists(DIR),'%s not exist'%DIR
    DIR = DIR.rstrip('/')
    if not os.path.isdir(DIR):
        return [DIR]
    else:
        cmd = 'ls -LR %s'%DIR
        res = subprocess.check_output(cmd,shell=1)
        res = re.split(r'([^\n]*):',res)[1:]
        it = nTuple(res,2,silent=0)
        DIR, ss = it[0];
        for dd,ss in it[1:]:
            NEWDIR, ALI = dd.rsplit('/',1)
            assert NEWDIR == DIR, 'Next directory %s not contained in %s'%(dd,DIR)
            DIR = dd 
        res = [ '%s/%s'%(DIR,x) for x in LinesNotEmpty(ss)]
        return res                

retype = type(re.compile('hello, world'))
def revSub(ptn, dict):
    '''Reverse filling a regex matcher.
    Adapted from: https://stackoverflow.com/a/13268043/8083313
'''
    if isinstance(ptn, retype):
        ptn = ptn.pattern
    ptn = ptn.replace(r'\.','.')
    replacer_regex = re.compile(r'''
        \(\?P         # Match the opening
        \<(.+?)\>
        (.*?)
        \)     # Match the rest
        '''
        , re.VERBOSE)
    res = replacer_regex.sub( lambda m : dict[m.group(1)], ptn)
    return res

def write_log(fname, s):
    f = open(fname, 'a')
    f.write(s + '\n')
    f.close()

def gnuPara(cmd,debug=0,ncore = 6):
    '''
    [Deprecated] Bad and does not wait for tasks to finish
    '''
    tmp = tempfile.NamedTemporaryFile(delete=True) if not debug else open('temp.sh','w')
    with tmp as tmp:
        print cmd
        tmp.write(cmd)
        E = shellexec('parallel --gnu -j%d <%s &>>parallel.log'%(
            ncore,
            tmp.name
            )
        )
    return E

def mp_para(f,lst,ncore = 6):
    if ncore ==1:
        res = map(f,lst)
    else:
        p = mp.Pool(ncore)
        res = p.map_async(f,lst,)
        res = res.get(10000000)
        p.close()
        p.join()
    return res

datenow = lambda: datetime.datetime.now().strftime("%Y_%m_%d_%H:%M:%S")

#### Regex for downloaded .fastq(.gz) files
# PTN = re.compile('(?P<lead>.*)_S(?P<sample>\d{1,3})_L(?P<chunk>\d+)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)')
PTN = re.compile('(?P<lead>.*)_L(?P<chunk>\d+)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)')


def shellexec(cmd,debug=0):
    print(cmd) 
    if not debug:
        return subprocess.call(cmd,shell=1)
#         return os.system(cmd)

def process_rna_sample(samplePATH, debug=0):
    '''
    Pull together raw reads from an input folder
    Args:
        samplePATH: Folder of .fastq(.gz) fot. be processed
    Comment: Refactored based on Hui's map-RNA-seq.py process_rna_sample().    
    '''
    #     return os.system('/bin/bash -c `%s`'%cmd)
    #     cmd = '/bin/bash -c `%s`'%cmd

    #     return subprocess.call(cmd,env=os.environ,cwd=os.getcwd(),
    #                           shell=True)
    

    samplePATH = samplePATH.rstrip('/')
    shellexec('echo $SHELL')
    
    RNA_SEQ_MAP_FILE = 'some-script.sh'
    DESTINATION_DIR  ='"/path/to/output/"' 
    WORKING_DIR='.'
    
    ### Extract  RunID from samplePATH
    samplePATH = samplePATH.rstrip('/')
#     ptn = '[\^/](\d{1,4}[RC][_/].*)'
#     ridPath = re.findall(ptn,samplePATH)
    ridPath = re.findall(ptn.runCond,samplePATH)
    assert len(ridPath)==1,'[ERROR] Cannot extract RunID from path name:"%s"'%samplePATH
    ridPath = ridPath[0][-1]
    print '[ridPath]',ridPath

    # Create a temporary directory 
    os.system('mkdir -p %s'%WORKING_DIR)
    temp_dir = os.path.join(WORKING_DIR,
                            '%s-%s'%(
                                ridPath.replace('/','-'),
#                                 os.path.basename(samplePATH),
                                datenow(),
                            )
    )
    os.system('mkdir -p %s'%temp_dir)

    

    
    #### Download raw read .fastq from samplePATH
#     print samplePATH
    FILES = glob.glob('%s/*' % samplePATH)
    FILES = sum(map(LeafFiles,FILES),[])
#     ccmd = '%s/* -t %s'%(samplePATH,temp_dir) 
    ccmd = '%s -t %s'%(' '.join(FILES), temp_dir) 
    cmd1 = 'cp -lr %s'%ccmd; 
    cmd2 = 'cp -r %s'%ccmd
    shellexec(cmd1) ==0 or shellexec(cmd2) 
    ODIR = os.getcwd()
    print '[ODIR]',ODIR
    try:
        os.chdir(temp_dir) #     shellexec('cd %s'%temp_dir)

        #### Parse .fastq filenames and assert quality checks
#         print '[MSG] found leaves','\n'.join(FILES)
        if debug:
            FS = [x.rsplit('/')[-1] for x in  FILES]
            print FS[:5]
#             FS = [x[pL+1:] for x in FILES]
#             FS = FILES
    #         assert 0
        else:
            FS = glob.glob('*')
        BUF = '\n'.join(FS)
        PARSED = [{'fname':m.group(0),'data':m.groupdict()} for m in re.finditer(ptn.baseSpace,BUF)]
        for d in PARSED:
            d['data']['fname'] = d['fname']
        data = [x['data'] for x in PARSED]
        meta = pd.DataFrame(data)

        if debug:
            return meta
        else:
            pass
        meta = check_L004(meta)
        meta = meta.sort_values(['lead','read','chunk'])
        unzipAndConcat(meta)
            
        print '[DONE!]:%s'%samplePATH

        os.system('echo %s >OLDDIR'%ridPath)
#         exit(0)
    except Exception as e:        
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        raise e
    finally:
        os.chdir(ODIR)
    print '[[WTFFF1]]'
    #### Stop here
    return temp_dir

def check_L004(meta):
    g = meta.groupby(['lead','read'],as_index=0)
    ct = g.count()

    mout = meta.merge(ct[['lead','read','chunk']] ,on=['lead','read'],suffixes=['','_count'])
    idx = mout['chunk_count'] ==4
    if not idx.all():
        print '[WARN] following reads are discarded due to chunk_count != 4'
        print mout[~idx][['fname','chunk']] 
        mout = mout[idx]
    return mout
def unzipAndConcat(meta,debug= 0):
    idx= [x.endswith('gz') for x in meta['ext']]
    if any(idx):
        #### unzip .gz where applicable
        mcurr = meta.iloc[idx]
        cmds = [cmd_ungzip(x) for x in mcurr['fname']]
        if debug:
            print '\n'.join(cmds[:1])
        else:
            mp_para(shellexec,cmds, ncore=NCORE)            
        #### Remove .gz in DataFrame accordingly
        meta.loc[idx,'ext'] = [ x.rstrip('.gz')  for x in mcurr['ext'] ]

    ### Map metas to fnames after decompression 
    mapper = lambda x: revSub(ptn.baseSpace,x)
    meta['fname'] = meta.apply(mapper,axis=1)
    g = meta.groupby(['lead','read'])
    cmds = [cmd_combineFastq(x[1]['fname']) for x in g]
    if debug:
        print '\n'.join(cmds[:1])
    else:
        mp_para( shellexec,cmds, ncore=NCORE)
#     os.system('sleep 5;')
    return 

def cmd_combineFastq(fnames,run=0):
    fnames = sorted(list(fnames))
    d = ptn.baseSpace.match(fnames[0]).groupdict()
    cmd = 'cat {IN} >{lead}_R{read}_raw.{ext} ; sleep 0; rm {IN} '.format(IN=' '.join(fnames),
                                                 **d)
    return cmd
def cmd_ungzip(F,):
    cmd = 'gzip -d <{IN} >{OUT} ; sleep 0 ; rm {IN} '.format(IN=F,OUT=F.rstrip('.gz'))
    return cmd

assert len(sys.argv) >= 2,'''
    Usage: (python) map-RNA-seq.py /path/to/folder/
        The folder should contains raw reads in .fastq(.gz) format
'''

if __name__=='__main__':
    NCORE = int(os.environ.get('NCORE',6))
    print '[NCORE]=',NCORE
    # NCORE = 1
    samplePATH = sys.argv[1]
    temp_dir = process_rna_sample( samplePATH, )
#     for i in range(10):
#         os.system('sleep 0.5')
#         print i*0.5
    # raise Exception('[WTF]%s'%temp_dir)
    print >>sys.stdout,temp_dir
    sys.exit(0)
