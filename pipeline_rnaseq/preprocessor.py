#!/usr/bin/env python2

# coding: utf-8

# In[272]:


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
    p = mp.Pool(ncore)
    res = p.map_async(f,lst,)
    res = res.get(10000000)
    p.close()
    return res

datenow = lambda: datetime.datetime.now().strftime("%Y_%m_%d_%H:%M:%S")

#### Regex for downloaded .fastq(.gz) files
# PTN = re.compile('(?P<lead>.*)_S(?P<sample>\d{1,3})_L(?P<chunk>\d+)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)')
PTN = re.compile('(?P<lead>.*)_L(?P<chunk>\d+)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)')


def shellexec(cmd,debug=0):
    print(cmd) 
    if not debug:
        return os.system(cmd)

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
    
    # Create a temporary directory 
    os.system('mkdir -p %s'%WORKING_DIR)
    temp_dir = os.path.join(WORKING_DIR,'%s-%s'%(datenow(),
                                                os.path.basename(samplePATH)))
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
    shellexec('cd %s'%temp_dir)
    ODIR = os.getcwd()
    print '[ODIR]',ODIR
    try:
        os.chdir(temp_dir)

        #### Parse .fastq filenames and assert quality checks
        if debug:
            FS = [x.rsplit('/')[-1] for x in  FILES]
            print FS[:5]
#             FS = [x[pL+1:] for x in FILES]
#             FS = FILES
    #         assert 0
        else:
            FS = glob.glob('*')
        BUF = '\n'.join(FS)
        PARSED = [{'fname':m.group(0),'data':m.groupdict()} for m in re.finditer(PTN,BUF)]
        for d in PARSED:
            d['data']['fname'] = d['fname']
        data = [x['data'] for x in PARSED]
        meta = pd.DataFrame(data)
        meta = check_L004(meta)
        
        if debug:
            return meta
        else:
            unzipAndConcat(meta)
#             print '\n\n[BUF]',BUF
#             print PARSED
#             return data
#         meta =  pd.DataFrame(index = data)

#         R1 = [d for d in PARSED if d['data']['read']=='1']
#         R2 = [d for d in PARSED if d['data']['read']=='2']
#         Rboth = R1+R2
#         common_names = set(d['data']['lead'] for d in Rboth)
#         assert len(common_names) == 1,'Common leading strings are not unique: %s' % common_names
#         common_name = common_names.pop()
#         alias = common_name

#         CHU1=sorted(d['data']['chunk'] for d in R1)
#         CHU2=sorted(d['data']['chunk'] for d in R2)
#         assert  CHU1 == CHU2,'Counts of R1/R2 chunks disagree, R1:%s R2:%s '%(CHU1,CHU2) 
#         assert len(CHU1) >= 4,'Counts of R1/R2 chunks: Actual %d Expected: >= 4 ' %len(CHU1) 

#         R1name = sorted((d['fname'] for d in R1),)
#         R2name = sorted((d['fname'] for d in R2),)
#         print R1name,'\n',R2name ### debug printout

#         #### Unzip where required
#         FS = [F for F in [d['fname'] for d in Rboth] if F.endswith('gz')]
#         cmd = '\n'.join(['gzip -d <%s >%s; sleep 0'% (F,F.rstrip('.gz')) 
#                          if not os.path.exists(F.rstrip('.gz')) 
#                          else '## gzip -d skipped since fastq exists for %s' % F for F in FS])

#     #     gnuPara(cmd,debug=0,ncore= 1)
#         mp_para(shellexec,cmd.splitlines(),ncore=NCORE)
#         R1name = [n.rstrip('.gz') for n in R1name]
#         R2name = [n.rstrip('.gz') for n in R2name]    

#         cmd = 'cat %s >%s_R1_raw.fastq' % (' '.join(R1name)  ,alias)
#         cmd +='\ncat %s >%s_R2_raw.fastq' % (' '.join(R2name),alias)
#         shellexec(cmd)

        print '[DONE!]:%s'%samplePATH
        print temp_dir
        samplePATH = samplePATH.rstrip('/')
#         idPath = '/'.join(samplePATH.split('/')[-3:])
        ptn = '[\^/](\d{1,4}[RC][_/].*)'
        idPath = re.findall(ptn,s)[0]
        os.system('echo %s >OLDDIR'%idPath)
        
#         exit(0)
    except Exception as e:        
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        raise e
    finally:
        os.chdir(ODIR)
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
            print cmds[:1]
        else:
            mp_para(shellexec,cmds, ncore=NCORE)
            
        #### Remove .gz in DataFrame accordingly
        meta.loc[idx,'ext'] = [ x.rstrip('.gz')  for x in mcurr['ext'] ]

    ### Map metas to fnames after decompression 
    mapper = lambda x: revSub(PTN,x)
    meta['fname'] = meta.apply(mapper,axis=1)
    # meta['fname'] = meta['
    g = meta.groupby('lead')
    cmds = [cmd_combineFastq(x[1]['fname']) for x in g]
    if debug:
        print cmds[:1]
    else:
        mp_para(shellexec,cmds, ncore=NCORE)
    return 

def cmd_combineFastq(fnames,run=0):
    fnames = list(fnames)
    d = PTN.match(fnames[0]).groupdict()
    cmd = 'cat {IN} >{lead}_R{read}_raw.{ext}; rm {IN};'.format(IN=' '.join(fnames),
                                                 **d)
    return cmd
def cmd_ungzip(F,):
    cmd = 'gzip -d <{IN} >{OUT}; rm {IN}; '.format(IN=F,OUT=F.rstrip('.gz'))
    return cmd


# In[ ]:


assert len(sys.argv) >= 2,'''
    Usage: (python) map-RNA-seq.py /path/to/folder/
        The folder should contains raw reads in .fastq(.gz) format
'''
NCORE = os.environ.get('NCORE',6)
samplePATH = sys.argv[1]
temp_dir = process_rna_sample( samplePATH, )

exit(0)


# In[16]:


s = '/media/pw_synology3/PW_HiSeq_data/RNA-seq/Raw_data/testONLY/133R/BdPIFs-32747730/133E_23_DN-40235206'
import re
ptn = '[\^/](\d{1,4}[RC][_/].*)'
idPath = re.findall(ptn,s)[0]
print idPath


# In[15]:


SELF='preprocessor'
if __name__=='__main__':
    get_ipython().system(u'jupyter nbconvert --to python {SELF}.ipynb')
    get_ipython().system(u' mv {SELF}.py tmp.py; echo \\#!/usr/bin/env python2 >preprocessor.py; ')
    get_ipython().system(u' cat tmp.py>>{SELF}.py; rm tmp.py')


# In[205]:


path = '/home/feng/syno3/PW_HiSeq_data/ChIP-seq/Raw_data/182C/Bd_ELF3-44645602/FASTQ_Generation_2018-06-06_03_43_21Z-101158414/'
res = process_rna_sample(path,debug=1)

NCORE = 4
res = process_rna_sample(path,debug=0)

