
#### Data I/O
import CountMatrix as scountjo
import pymisca.util as pyutil
import pyBigWig as pybw
pd = pyutil.pd; np = pyutil.np; 
import pymisca.util as pyutil
# import synotil.dio as sdio


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
10:img
'''.strip().splitlines()
bedHeader = [x.split(':')[1] for x in bedHeader] 


def readChipPeak(fname,**kwargs):
    '''CountMatrix.geneList is buggy at the moment hence needed to be converted back into 
    CountMatrix.countMatrix
'''
    df = scount.countMatrix.from_DataFrame(fname=fname,ext='tsv',index_col='geneAcc',addFname=0,**kwargs)
    df= df.toGeneList()
    return df

import StringIO
def guessBedHeader(fname,silent=True,ext='tsv',guess_index=0, 
                   prefix = '', **kwargs):
    cmd = 'head -n5 %s'%fname
    buf = StringIO.StringIO(pyutil.shellexec(cmd,silent=silent))
    df = pyutil.readData(buf,ext=ext,header=None,guess_index=guess_index,**kwargs)
    if len(df.columns) >  len(bedHeader):
        header = bedHeader + list(df.columns)[len(bedHeader):]
    else:
        header = bedHeader[:len(df.columns)]
    if prefix: 
        header = ['%s_%s'%(prefix,x) for x in header]
    return map(str,header)

def header_closest(peak1,peak2):
    '''Get the header for bedtools intersect -wo /closest -d
'''
    header = sum( [guessBedHeader(x,prefix=k) for k,x in 
                   [('',peak1),
                    ('feat',peak2)]
                  ],
                 [])
    header +=['distance']
    return header

def extract_peak(fname,ext='tsv',header=None,guess_index=0,**kwargs):
    df = pyutil.readData(fname,ext=ext,header=header,guess_index=guess_index,**kwargs)
    df.columns = bedHeader[:len(df.columns)] + list(df.columns)[len(bedHeader):]
    return df
    


def npk_expandSummit(fname = None,df=None,radius=200,clip = 1,center_summit= 0):
    '''
    Expand the summit regions of a .narrowPeak dataFrame
'''
    if df is None:
        df = extract_peak(fname)
        
    if 'abs_summit' not in df.columns:
#         if 'summit' not in df.columns and :
        if center_summit:
            df['summit'] = (df.start + df.end )// 2
            
        assert 'summit' in df.columns
        if (df.summit >= df.start).all():
            pass
        else:
            df['summit'] = df.start + df.summit
        df.rename(columns={'summit':'abs_summit'}, inplace=True)
        df.abs_summit = df.abs_summit.astype('int')
#       df = df.rename(columns={'summit':'abs_summit'}, )
        
#     st = df.strand
    df.start = (df.abs_summit - radius)
    df.end  = df.abs_summit + radius 
#     df.drop('abs_summit',1,inplace=True)
    if clip:
        df.start = df.start.clip_lower(0)
        
    if fname is not None:
        base = pyutil.basename(fname)
        ofname = '%s_radius=%d.tsv' % (base,radius)
        df.to_csv(ofname,sep='\t',index=None,header=None)
        return ofname
    else:
        return df

def msg_bwFile(fname):
#     import pyBigWig as pybw
    f = pybw.open(fname)
    print f.chroms()
    f.close()    
    return fname
#     df.columns = [pyutil.basename(fname)]
#     chipPeak = df
#     return df

def parseBedmap(fname = None, df = None, ):
    ''' Parse the output of bedMap
'''
    if df is None:
        df = pyutil.readData(fname,header = None,ext='tsv',guess_index=False)

    df = df.dropna()
    
    df.columns = bedHeader + ['hit']

    res = pyutil.explode(df,'hit','acc',';')
    res = res.merge(df.drop('hit',1),on='acc')
    return res

def extract_closest(fname = None, df = None, ):
    ''' Parse the output of 'bedtools closest'
'''
    if df is None:
        df = pyutil.readData(fname,header = None,ext='tsv',guess_index=False)
#     df = df.dropna()    

    header = bedHeader + pyutil.paste0([['feature_'], bedHeader]).tolist()
    df = df.iloc[:,:18]
    df.columns = header[:17] + ['distance']
    df['hit'] = df['feature_acc']
    return df

parseBedClosest = extract_closest 

def bed_randomise(infile,GSIZE = None,silent=1):
    '''Create a randomly distributed bed file
'''
    ofile = pyutil.basename(infile) + '_type=random.bed'
    assert GSIZE is not None
    LC = pyutil.lineCount(infile)
    cmd = "bedtools random -g {GSIZE} -l 2 -n {LC}  | tee {ofile}".format(**locals())
    pyutil.shellexec(cmd,silent=silent)
    return ofile

def bed__summit(peakFile,GSIZE=None,silent=1,opt='-s -l -0.5 -r -0.5 -pct',
               inplace=True):
    if GSIZE is None:
        GSIZE=pyutil.os.environ.get('GSIZE',None)
    assert GSIZE is not None
    ofname = '%s.summit' % peakFile
    if not inplace:
        ofname = pyutil.os.path.basename(ofname)
    cmd = 'cat {peakFile} \
    | bedtools slop -g {GSIZE} {opt} -i - \
    > {ofname}'.format(**locals())
    pyutil.shellexec(cmd,silent=silent)
    return ofname
def bed__leftSummit(peakFile,opt='-s -pct -l 0. -r -1.0',**kwargs):
    return bed__summit(peakFile,opt=opt,**kwargs)
def bed__rightSummit(peakFile,opt='-s -pct -r 0. -l -1.0',**kwargs):
    return bed__summit(peakFile,opt=opt,**kwargs)
# bed__leftSummit = pyutil.functools.partial(bed__summit,
#                                            opt='-s -pct -l 0. -r -1.0')

def bed__addCol__interval(bed,
                          name='interval',
                          expr = '@pyutil.paste0([chrom, ":",start,["-"],end]).tolist()',
                         ):
    pyutil.df__addCol(bed, name, expr )
    return bed

import StringIO
import numpy as np

def bed__guessWidth(bedFile,silent=1,head=100): 
    res = pyutil.shellexec('head -n{head} {bedFile}'.format(**locals()),
                           silent=silent)
    dfc = extract_peak(StringIO.StringIO(res))
    span = (dfc.end - dfc.start).values.ravel()
    M = np.median(span)
    if span.std() / len(span)**0.5 > 0.1 * M:
        pyutil.sys.stderr.write ('[WARN]:estimation may be unstable\n')
    return int(M) 

def summitDist(peak1,peak2,
               CUTOFF = 400,
            silent = 1,
              GSIZE = None,
               as_fname=0,
              **kwargs):
    '''Find nearby summits within a distance cutoff
'''
    if GSIZE is None:
        GSIZE=pyutil.os.environ.get('GSIZE',None)
    assert GSIZE is not None
    RANGE  = CUTOFF//2 - 1
    infiles = [peak1,peak2]
#     def file_ncol(fname):
#         cmd = 'wc -l %s'%(fname)
#         res = pyutil.shellexec(cmd,silent=silent)
#         ncol = res[0].strip().split('\t')
#     incols = 
    incols = map(pyutil.file_ncol,infiles)

    ### padding/inflate the summit to have radius
    lst = []
    for infile in infiles:
        
        ofile = "{infile}.{RANGE}".format(**locals()).split('/')[-1]
        lst += [ofile]

        cmd = "bedtools slop -g {GSIZE} -b {RANGE} -i {infile} \
          | tee {ofile}".format(**locals())
        _ = pyutil.shellexec(cmd,silent=silent)

    slop1,slop2 = lst
    FOUT = 'infiles:'+ ":".join(map(pyutil.basename,infiles)) \
        + "__cutoff:{}.tsv".format(CUTOFF)

    # ### bed format 1=chrom, 2=start, 3=end
    # cols = ','.join(map(str,[2,3,] + [x + incols[0] for x in [2,3]]))
    # cmd = "bedtools closest -a {slop1} -b {slop2} \
    #   | bedtools overlap -cols {cols} \
    #   | tee {FOUT}".format(**locals())

    cmd = "bedtools intersect -wo -a {slop1} -b {slop2} \
      | tee {FOUT}".format(**locals())


    buf = pyutil.shellexec(cmd,silent=silent)

    ### [TBC]Memory-intensive, Replace with awk mutation in the future
    columns = header_closest(peak1,peak2)
    
    df = pyutil.readData( StringIO.StringIO(buf), header = None, ext='tsv',guess_index=False, 
                         columns = columns)
    df.distance = CUTOFF - df.distance
    df.to_csv(FOUT, sep='\t',index = False)
    if as_fname:
        return FOUT
    else:
        return df
    
def bed__merge(bedFile,silent=1,opt='-c 4 -o first'):
    bname = pyutil.os.path.basename(bedFile)
    path = pyutil.os.path.dirname(bedFile)
    ofname = pyutil.os.path.join(path,'merged__%s'%bname)
    cmd = 'bedtools sort -i {bedFile} |  bedtools merge -i - {opt} > {ofname}'.format(**locals())
    pyutil.shellexec(cmd,silent=silent)
    return ofname

def bed__makewindows(bedFile,windowSize=100,stepSize=None,silent =1 ):
    if stepSize is None:
# windowSize=100
        stepSize = windowSize//2
# bedFile = 'per-score-GT-0dot6_188C_RESEQ-combined.bed'
#         bedbase = pyutil.bname(bedbane)
    ofname = '{bedFile}.w{windowSize}s{stepSize}'.format(**locals())
    cmd = "bedtools makewindows -i srcwinnum -w {windowSize} -s {stepSize} -b {bedFile} > {ofname} ".format(**locals())
    pyutil.shellexec(cmd,silent = silent)
    return ofname


    
def extract_bigwig_worker(lines, bwFile = None,stepSize = 1, stranded= 1, bw = None):
    ''' Helper mapper for querying BigWig
'''
    bw = pybw.open(bwFile)
    chromL = bw.chroms()
    
    lines = [x for x in lines if x]
    nField = lines[0].strip().split('\t').__len__() 
    res = []
    for line in lines:
#     def parse(line, nField = nField):
        if line is None:
            return None
        cols = line.strip().split('\t')
        if nField >= 6:
            chrom, start,end, (id, score, strand) = cols[0],int(cols[1]),int(cols[2]),cols[3:6]
        else:
            strand = '+'
            if nField is 5:
                chrom, start,end, id,_ = cols[0],int(cols[1]),int(cols[2]),cols[3],cols[4]
                
#                 assert 0, 'operation not defined when bedFile has 5 fields:\n%s'%lines[0]
            elif nField is 4:
                chrom, start,end, id = cols[0],int(cols[1]),int(cols[2]),cols[3]
            else:
                chrom, start,end = cols[0],int(cols[1]),int(cols[2])
                id = 'NoID'
                
        if chrom not in bw.chroms():
            o = None
        else:
            start = max(0,start)
            end = min(chromL[chrom],end)
            sec = bw.values(chrom, start, end, numpy=0)
            if strand is not '-' or not stranded:
                vals = sec[::stepSize]
            else:
                vals = sec[::-stepSize]
                
            o = vals
#         return (id,o)
        res+=[(id,o)]
#     res = map( parse, lines) 
    bw.close()
    return res

def extract_bigwig(bwFile,bedFile,
                   stepSize=1,
                  mapChunk = None, 
#                   span = None
                   shift=1,
#                   outIndex = None,
                   stranded=1,
                 ):
    ''' Extracting a signal matrix for each bed region
'''
#     assert NCORE == 1,'Multi-thread is slower here..., so dont! '

#     assert stepSize == 1,'Not implemented'        
    with pybw.open(bwFile) as bw:
        it = open(bedFile)
        worker = pyutil.functools.partial(extract_bigwig_worker,
                                          bwFile =bwFile,
                                         stepSize=stepSize,
                                         stranded=stranded,)
        if  1 == 1:
            res = map(worker,[it])

            
        res = sum(res,[])
#             pass 
        ids, out  = zip(*res)

    #### Replacing "None" and incomplete intervals
    ref = next((item for item in out if item is not None),None)
    assert ref is not None,'Cannot find an reference shape, likely wrong chromosomes.\n\
    bigwigFile:"%s" '%bwFile
#     L = len(ref)
#     L = len(res) if span is None else span //stepSize        
    L = max(map(len,out))
    lst = []
    print '[L]=',L
    for x in out:
        if x is None:
            y = [0.]*L
        else:
            Lx = len(x)
            y = x + [0.] * (L - Lx)            
        lst += [y]
#         out = [[0.]*L if x is None else x for x in out]
    out = np.array(lst)
    out = np.nan_to_num(out)
    
#     MLEN = np.mean([len(x) for x in out]) 
    MLEN='not set'
    assert out.dtype!='O','''Unable to shape the matrix properly: 
    %s, %s '''% (MLEN, [(type(x),x) for x in out if len(x)< MLEN] )
    out = pd.DataFrame(out).set_index([list(ids)])
        
    cols = stepSize * ( np.arange(0, out.shape[-1], )  ) 
    if shift:
        mid = ( L * stepSize )//2
        cols += -mid
    out.columns = cols
    
#     out.columns = (stepSize * np.arange(0, out.shape[-1], ))
            # Do something with the values...
    

#     out = ctMat.countMatrix.from_DataFrame(df=out)
#     out.fname = bwFile
    out.param = pyutil.util_obj()
    out.param['bwFile'] = bwFile
    out.param['bedFile'] = bedFile
    return out

def extract_bigwig_multiple(bedFile= None,peakFile = None,
                            bwFiles = None,fnames=None,
                            radius = None,
                            callback= lambda x:pd.concat(x,axis=1),
                            NCORE=1,
                            stepSize=20,
                            outIndex=pyutil.basename,
                            center_summit = 0,**kwargs):
    '''extract_bigwig() for multiple bwFiles
'''
    #########################
    #### Legacy support #####
    if bedFile is None:
        assert peakFile is not None
        bedFile =  peakFile
    else:
        assert peakFile is None
        
    if bwFiles is None:
        assert fnames is not None
        bwFiles = fnames
    #### Legacy support ####
    #########################
        
    if radius is not None:
        bedFile = npk_expandSummit(fname=bedFile,
                                   radius=radius,center_summit =center_summit)
    print '[L]',pyutil.lineCount(bedFile)

    #### Compute Matrix

    worker = pyutil.functools.partial(
        extract_bigwig,
#         bwFile=fname,
        stepSize=stepSize,
        bedFile=bedFile,
#         outIndex=outIndex,
#         NCORE=1,
        **kwargs) 
    bws = pyutil.mp_map(worker, bwFiles, n_cpu = NCORE)
#     bws = [ for fname in fnames]

    for i,(bwFile,out) in enumerate(zip(bwFiles,bws)):
        replaceCol  = None
        if outIndex is None:
            replaceCol = None
        else:
            if callable(outIndex):
                replaceCol = outIndex(bwFile)
            else:
                replaceCol = outIndex[i]
        if replaceCol is not None:
            #### replace outer index,replace with rename()
            #### [TEMP]
            tmp = out.T
            tmp['ind'] = replaceCol
            tmp.set_index( 'ind', append=1,inplace=True)
            tmp = tmp.reorder_levels(['ind',None])        
            tmp.index.names = ['bwFile','pos']
            out = tmp.T
        out.index.name = 'acc'
        bws[i] = out
        
    dfc =  bws
    
    
#     dfc=pd.concat(bws,)
#     dfc = pd.concat(bws,axis=1) #### abosorbed into callback()
    
    if callback is not None:
        dfc = callback(dfc)
    return dfc

# def bigwig_average_signal(dfc,):
#     L = len(dfc)
#     dfc = pyutil.colGroupMean(dfc)
#     dfc = scount.countMatrix(dfc.T)
#     dfc.param['n_peak'] = L
#     return dfc


sdio = pyutil.sys.modules[__name__]
def job__nearAUG(peakFile = None,featFile =None,
                 peakSummit= None, featSummit = None, 
                 CUTOFF=6000,peakWid=None,GSIZE=None):
    JOB='nearAUG'
    if peakSummit is None:
        assert peakFile is not None
        peakSummit = sdio.bed__summit(peakFile,GSIZE=GSIZE,inplace=0)
    if featSummit is None:
        assert featFile is not None
        featSummit = sdio.bed__leftSummit(featFile,GSIZE=GSIZE,inplace=0)

    if peakWid is None:
        peakWid= sdio.bed__guessWidth(peakFile)//2 
        
    res = sdio.summitDist(
        peakSummit,
        featSummit,
        CUTOFF= 6000 - peakWid,
        GSIZE=GSIZE
    )
    
    out = sdio.extract_peak(peakSummit).merge(
        res[['acc','feat_acc','distance']],
#         res.drop(columns=['chrom','start','end']),
        how='right',
        left_on='acc',right_on='acc').query("distance < %d" % CUTOFF)

    featSummitBase = pyutil.basename(featSummit).replace('_','-')
    peakSummitBase = pyutil.basename(peakSummit).replace('_','-')
    ofname = '\
job_{JOB}__\
peak_{peakSummitBase}__\
cutoff_{CUTOFF}__\
feat_{featSummitBase}.tsv'.format(**locals())
#     pyutil.to_tsv(out,ofname)   
    out.to_csv(ofname,sep='\t')
#     pyutil.to_tsv(out,ofname,header=True,index=1)   
    return ofname


###### [TBC] belongs to LookupTable in the future
def peak2gene(pgf,query,**kwargs):
    res = pgf.merge(query,left_on='acc',right_index=True,**kwargs)
    return res
def gene2peak(pgf,query,**kwargs):
    res = pgf.merge(query,left_on='feat_acc',right_index=True,**kwargs)
    return res

def df__sanitise__Ath(df):
    df.index = df.index.str.extract(".*(AT\dG\d{5}).*",expand=False)
    return df
####

def getChrom(fname):
    f = pybw.open(fname)
    chrom = f.chroms()
    f.close()
    return chrom

def check__bigwig__chrom(fnames,chrom=None):
    
    chroms = map(getChrom, fnames)
    if chrom is None:
        chrom = chroms[0]
    good = map(lambda x: x==chrom, chroms)
    return pd.Series(good,index=fnames,name='check__bigwig__chrom')    


def assign__filename(lst,runID='000R',ext='excel.count',DIR='.',init=1):
    '''
    Internal routine for creating data structure
'''
    out = []
    DIR = DIR.rstrip('/')
    fmt = '{DIR}/{runID}/S{sampleID_int}/{ele}_S{sampleID_int}.{ext}'
    for i,ele in enumerate(lst):
        sampleID_int = i + init
        ele = ele.replace('_','-')
        res = fmt.format(**locals())
        out +=[res]
    return out

def df__deposit(dfc,runID='000R',ext='excel.count',DIR='.',init=1,silent=0,
               sep='\t',header=1):
#     dfc = rnaseq.copy()
    dfc = dfc.copy()
    fnames = assign__filename(dfc.columns,
                              runID=runID,ext=ext,DIR=DIR,init=init)

    dfc.columns = fnames
    dfc.columns.name = 'fname'
    gp = dfc.reset_index().melt(value_name='TPM',id_vars=['gene_id']).groupby('fname')
    for fname,df in gp:
        odf = df.drop(columns='fname').set_index('gene_id')
        pyutil.shellexec('mkdir -p `dirname {fname}`'.format(**locals()),silent=silent,)
        odf.to_csv(fname,sep=sep,header=header)
    return dfc.columns


idKeys = ['runID','sampleID','read']

def rawFile__validateChunk(dfc,):
    '''Validate the df_raw
'''
    idKeys = ['runID',
              'sampleID',
              'read']
    dfc = dfc.sort_values(idKeys + ['chunk'])
    gp  = dfc.groupby(idKeys)
    for (key,df) in gp:
        assert len(df) == 4,key
    if not dfc['ext'].iloc[0].startswith('.'):
        dfc['ext'] = dfc['ext'].map(lambda x:'.%s'%x)
    dfc['fnameCombined'] = pyutil.df__paste0(dfc,
                                              idKeys + ['ext'],
#                                               headerFmt='_',
                                              sep='_',
                                              ).tolist()
    dfc['fnameCombinedSize'] = 0
    return dfc
def rawFile__combineChunk(dfc,silent=0):
    ''' combine chunkedFiles in df_raw
    according to "fname" and "fnameCombined"
'''
    dfc = dfc.copy()
#     dfc = dfc[idKeys + ['fnameCombined','fnameCombinedSize']].drop_duplicates()
    
    fnameFlat = ' \\\n'.join(dfc.fname)
    ofnames = dfc.fnameCombined.unique()
    assert len(ofnames)==1,\
    'contains mulitple fnameCombined!:%s'%ofnames
    ofname = ofnames[0]
    cmd = 'cat {fnameFlat} > {ofname}'.format(**locals())
    res= pyutil.shellexec(cmd,silent=silent)
    dfc['fnameCombinedSize'] = pyutil.os.path.getsize(ofname)
    dfc = dfc[idKeys + ['fnameCombined','fnameCombinedSize']
             ].drop_duplicates()
    return dfc

def worker__rawFile__combineChunk((key,dfc)):
    '''not fast, dont use'''
    return rawFile__combineChunk(dfc,silent=0)

def bam__getHeader(fname,grepKey='SQ',silent=1,head=100):
    cmd = u'samtools view -H %s' % fname 
    if grepKey is not None:
        cmd = u'{cmd} | grep {grepKey}'.format(**locals())
    if head is not None:
        cmd = u'{cmd} | head -n{head}'.format(**locals())
    res = pyutil.shellexec(cmd,silent=silent)
    return res

def listByChip(dfc):
    '''Transform a bwTracks object
'''
    dfcc = pd.DataFrame(iter(dfc.groupby(axis=1,level='bwFile')),
                        columns=['bwFile', 'tracks'],
                       )
    for df in dfcc['tracks']:
        key = df.columns.levels[0][0]
        df.columns = df.columns.droplevel(0)
        df.index =  pd.MultiIndex.from_arrays([
            [key]*len(df),
            df.index],
            names=['bwFile','acc'],
        )
#         df['bwFile'] ='test'
#         df.set_index('bwFile', append=True,inplace=True)
#     dfcc['tracks']  =dfcc['tracks'].map(lambda x:
#                                         x.columns=x.columns.drop_level())
    return dfcc

def listByGene(dfc):
    ''' Transform a bwTracks object
'''
    def pivotRecord(x):
        res = x.melt().pivot_table(index='bwFile',
                                         columns='pos',)
        return res
    dfcc = dfc.groupby(axis=0,level=0).apply(pivotRecord)
#     dfcc.index = dfcc.index.rename(level=0,names='acc')
    return dfcc

def listByGene(dfc,concat=1):
#     idvars = dfc.columns.names
    if not concat:
        raise Exception('Not implemented')
        
    dfcc = dfc.T.reset_index().melt(id_vars = ['bwFile','pos']).pivot_table(
        index=['acc','pos'],
        columns='bwFile')
    return dfcc

# def listByPos(dfc):
#     ''' Transform a bwTracks object
# '''
#     def pivotRecord(x):
#         res = x.melt().pivot_table(index='bwFile',
#                                    columns='pos',
#                                          )
#         return res
#     dfcc = dfc.groupby(axis=0,level=0).apply(pivotRecord).T
#     return dfcc

def file__concat(bedFiles,silent=1,
                ofname = 'concated_file',
                ext = None
#                 opt='-c 4 -o first'
               ):
    '''Concatentat a list of files
'''
    if ext is None:
        sp = bedFiles[0].rsplit('.',1)
        if len(sp) == 2:
            ext = sp[-1]
        else:
            ext = None
    if ext is not None:
        ofname = '.'.join([ofname,ext])
    flatName = ' '.join(bedFiles)
    cmd = 'cat {flatName} >{ofname}'.format(**locals())
#     bname = pyutil.os.path.basename(bedFile)
#     path = pyutil.os.path.dirname(bedFile)
#     ofname = pyutil.os.path.join(path,'merged__%s'%bname)
#     cmd = 'bedtools sort -i {bedFile} |  bedtools merge -i - {opt} > {ofname}'.format(**locals())
    pyutil.shellexec(cmd,silent=silent)
    return ofname

def bed__totalLength(bedFile,silent=1):
    '''Source: https://www.biostars.org/p/68283/#68292
'''
    cmd = "cat %s | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'" % bedFile
    res = pyutil.shellexec(cmd,silent=silent)
    res = int(res.strip())
    return res

def bed__embed(outerBed,innerBed,debug=0,ofname = None,mergeAcc=1):
    if not isinstance(outerBed,pd.DataFrame):
        outerBed=  sdio.extract_peak(outerBed)
    if not isinstance(innerBed,pd.DataFrame):
        innerBed=  sdio.extract_peak(innerBed)
#             lc[key] = sdio.extract_peak(val)
    assert 'acc' in outerBed,'Reference bed file must be named'
    if 'acc' in innerBed:
        if mergeAcc:
            innerBed['acc'] = pyutil.df__paste0(innerBed,
                                                ['chrom','acc'],sep='_').tolist()
    for df in [outerBed,innerBed]:
        if 'strand' not in df.columns:
            df['strand'] = '+'
        df.strand.fillna('+')
        df['strandVal'] = df.strand.isin(['+','*'])
    res = innerBed.merge(outerBed,how='left',
                         left_on='chrom',
                         right_on='acc',
#                          prefixes=['t','1'],
                        suffixes=['Inner','Outer'],)
    res['strandValFinal'] = ~( (res['strandValOuter']) ^ (res['strandValInner']) )
    isNeg = res.strandValFinal==0
    (res.loc[ isNeg,'startInner'], res.loc[isNeg,'endInner']) \
        = (-res.loc[isNeg,'endInner'],-res.loc[isNeg,'startInner'])
    res['valid'] = res.eval('startInner<=endInner')
    assert res.valid.all()
#     print 'acc' in res.columns
    if debug ==2:
        return res

    if debug ==1:
        return res[['valid',
                    'chromInner',
                    'startInner',
                    'endInner',
                    'strandValFinal',
                    'strandValInner',
                    'strandValOuter']]
    res['shift'] = 0
    isOuterNeg = res.strandValOuter == 0
    res.loc[isOuterNeg, 'shift'] = res.loc[isOuterNeg,'endOuter']
    res.loc[~isOuterNeg, 'shift'] = res.loc[~isOuterNeg,'startOuter']
#     res.shift = r
    res['start'] = res['shift'] + res['startInner']
    res['end'] = res['shift'] + res['endInner']
    res['chrom'] = res['chromOuter']
    res['acc'] = res['accInner']
    resDF = res[['chrom',
                    'start',
                    'end',
                     'acc'
#                      'accInner'
                   ]]
    if ofname is not None:
        try:
            pyutil.to_tsv(resDF,ofname,)
            return ofname
        except Exception as e:
            print e
    return resDF

def tsv__getColumns(fname,ext='tsv'):
#     pyutil.readData()
    res = file__header(fname,silent=silent)
    df = pyutil.readData(res,ext=ext)
    return df.columns.tolist()

def file__header(fname,head = 10,silent=1):
    res = pyutil.shellexec('head -n{head} {fname}'.format(**locals()),
                          silent=silent)
    res = pyutil.StringIO.StringIO(res)
    return res

def count__getGeneHeader(fname, 
                         ext='tsv',
                         pipeline=None, silent=1, **kwargs):
    ext = 'tsv' ### hard set
    res = file__header(fname,silent=silent)
    df = pyutil.readData(res,ext=ext,guess_index=0)
    return df.gene_id.tolist()

def count__getGeneHeader__dict(data,head=10):
#     data['ext'] = data['pipeline']
    res =  count__getGeneHeader(head=head,**data)
    return res
    

def wig2bigwig(fname,chromSizes = 'chrom.sizes',silent=1):
    ofbase = pyutil.getBname(fname) 
    ofname = '%s.bw' % ofbase
    cmd = '''wigToBigWig {fname} {chromSizes} {ofname}
    '''.format(**locals())
    res = pyutil.shellexec(cmd,silent=silent)
    return ofname

def clu2bed(segDF, ofname=None):
    '''Must have columns: ('acc','pos','clu')
    '''
    segDF = segDF.reset_index()
#     stdout,isFile = get__stdout(ofname)
    stepSize = np.diff(segDF['pos'].values[:2], axis=0)[0]
    vals = segDF[['clu','acc']].values    
    isDiff = (vals[1:] != vals[:-1]).any(axis=1)
    segDF['isDiff'] = np.concatenate([
                [True],
                isDiff],axis=0)
    it = (pyutil.util_obj(**vars(x)) for x in segDF.itertuples())
    peak = pyutil.collections.OrderedDict(
        (
        ('chrom',None),
        ('start',None),
        ('end',None),
        ('acc',None),
        )
    )
    peaks = []
    def savePeakStart():
        peak['chrom'] = rec.acc
        peak['start'] = rec.pos
        return
    def savePeakEnd():
#         kk = loc
        peak['end'] = oldPos + stepSize
        peak['acc'] ='summitPos%d'%( ( peak['start'] + peak['end']) //2)
        assert peak['end'] > peak['start'], peak
#         pyutil.ppJson(locals())
        peaks.append(peak.copy())
        
#         line = u'\t'.join(map(unicode,peak.values()))
#         stdout.write(u'%s\n'%line)

#         print peak
        return
    def changed():
        if idx != 0:
            if oldClu == 1:
                savePeakEnd()
            if rec.clu == 1:
                if (oldClu == 0) | (oldAcc != rec.acc):
                    savePeakStart()
        else:
            if rec.clu == 1:
                savePeakStart()
        return
    
    #### Starting the loop
    oldClu = 0 
    for idx,rec in enumerate(it):
        if (idx==0):
            changed()
        elif (rec.clu!=oldClu) or (rec.acc!=oldAcc):
            changed()
        oldClu = rec.clu    
        oldPos = rec.pos
        oldAcc = rec.acc
    changed()
    
    resDF = pd.DataFrame(peaks)
    
    if ofname is not None:
        try:
            pyutil.to_tsv(resDF,ofname,)
            return ofname
        except Exception as e:
            print e
    return resDF


# sizeDF.columns = ['chrom','length']
def bed__checkValid(bed, GSIZE, force=0):
    fname = None
    if not isinstance(bed,pd.DataFrame):
        fname = bed
        bed = sdio.extract_peak(bed)
    sizeDF = pyutil.readData(GSIZE,ext='tsv',header=None,guess_index=0)
    sizeDF.columns = ['chrom','length']
    bedDF =  sizeDF.merge(bed)
    bedDF['valid'] = bedDF.eval('start > 0 and end <= length')
    if not force:
        assert bedDF.valid.all()
    else:
        resDF = bedDF.query('valid').drop(columns=['valid','length'])
        if fname is not None:
            ofname = '%s__valid.bed' % fname.rsplit('.',1)[0]
            pyutil.to_tsv(resDF,ofname)
            return ofname
        else:
            return resDF