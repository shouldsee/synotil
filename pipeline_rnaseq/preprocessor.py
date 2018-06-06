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


def make_directory(dir):
    if not os.path.isdir(dir):
        try:
            os.makedirs(dir)
        except:
            raise Exception('Cannot make directory %s. STOP.' % (dir))
        

# shellexec = os.system
def shellexec(cmd):
    print(cmd) 
    return os.system(cmd)
#     return os.system('/bin/bash -c `%s`'%cmd)
#     cmd = '/bin/bash -c `%s`'%cmd

#     return subprocess.call(cmd,env=os.environ,cwd=os.getcwd(),
#                           shell=True)

def check_all_samples(d):
    '''
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

def process_rna_sample(samplePATH, debug=0):
    '''
    Pull together raw reads from an input folder
    Args:
        samplePATH: Folder of .fastq(.gz) fot. be processed
    Comment: Refactored based on Hui's map-RNA-seq.py process_rna_sample().    
    '''
    

    samplePATH = samplePATH.rstrip('/')
    shellexec('echo $SHELL')
    
    RNA_SEQ_MAP_FILE = 'some-script.sh'
    DESTINATION_DIR  ='"/path/to/output/"' 
    WORKING_DIR='.'
    
    # Create a temporary directory 
    make_directory(WORKING_DIR)
    temp_dir = os.path.join(WORKING_DIR,'%s-%s'%(datenow(),
                                                os.path.basename(samplePATH)))
    make_directory(temp_dir)

    #### Download raw read .fastq from samplePATH
    FILES =  glob.glob('%s/*.fastq*' % samplePATH)
    cmd1 = 'cp -l %s %s'%(' '.join(FILES),temp_dir)    
    cmd2 = 'cp %s %s'%(' '.join(FILES),temp_dir)    
    shellexec(cmd1) ==0 or shellexec(cmd2) 
    shellexec('cd %s'%temp_dir)
    os.chdir(temp_dir)
            
    #### Parse .fastq filenames and assert quality checks
    FS = glob.glob('*')
    PTN = r'(?P<lead>.*)_L(?P<chunk>\d+)_R(?P<read>[012])_(?P<trail>\d+)[^\.]+\.(?P<ext>.+)'
    BUF = '\n'.join(FS)
    PARSED = [{'fname':m.group(0),'data':m.groupdict()} for m in re.finditer(PTN,BUF)]
    R1 = [d for d in PARSED if d['data']['read']=='1']
    R2 = [d for d in PARSED if d['data']['read']=='2']
    Rboth = R1+R2
    common_names = set(d['data']['lead'] for d in Rboth)
    assert len(common_names) == 1,'Common leading strings are not unique: %s' % common_names
    common_name = common_names.pop()
    alias = common_name
    
    CHU1=sorted(d['data']['chunk'] for d in R1)
    CHU2=sorted(d['data']['chunk'] for d in R2)
    assert  CHU1 == CHU2,'Counts of R1/R2 chunks disagree, R1:%s R2:%s '%(CHU1,CHU2) 
    assert len(CHU1) >= 4,'Counts of R1/R2 chunks: Actual %d Expected: >= 4 ' %len(CHU1) 
    
    R1name = sorted((d['fname'] for d in R1),)
    R2name = sorted((d['fname'] for d in R2),)
    print R1name,'\n',R2name ### debug printout
    
    #### Unzip where required
    FS = [F for F in [d['fname'] for d in Rboth] if F.endswith('gz')]
    cmd = '\n'.join(['gzip -d <%s >%s; sleep 0'% (F,F.rstrip('.gz')) 
                     if not os.path.exists(F.rstrip('.gz')) 
                     else '## gzip -d skipped since fastq exists for %s' % F for F in FS])
    
#     gnuPara(cmd,debug=0,ncore= 1)
    mp_para(shellexec,cmd.splitlines(),ncore=NCORE)
    R1name = [n.rstrip('.gz') for n in R1name]
    R2name = [n.rstrip('.gz') for n in R2name]    
    
    cmd = 'cat %s >%s_R1_raw.fastq' % (' '.join(R1name)  ,alias)
    cmd +='\ncat %s >%s_R2_raw.fastq' % (' '.join(R2name),alias)
    shellexec(cmd)
    
    #### Stop here
    return temp_dir


assert len(sys.argv) >= 2,'''
    Usage: (python) map-RNA-seq.py /path/to/folder/
        The folder should contains raw reads in .fastq(.gz) format
'''
NCORE = os.environ.get('NCORE',6)
samplePATH = sys.argv[1]
temp_dir = process_rna_sample( samplePATH, )
print '[DONE!]:%s'%samplePATH
print temp_dir

samplePATH = samplePATH.rstrip('/')
last3path = '/'.join(samplePATH.split('/')[-3:])
os.system('echo %s >OLDDIR'%last3path)