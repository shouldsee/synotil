# import sy
# import scipy.cluster.hierarchy as sphier
# import scipy.spatial.distance as spdist
import pymisca.vis_util as pyvis
import pymisca.util as pyutil
plt = pyvis.plt; pd = pyutil.pd; np = pyutil.np
import pymisca.models as pymod
import pymisca.ext as pyext

import synotil.modelRoutine as smod
import synotil.CountMatrix as scount
import synotil.PanelPlot as spanel
import synotil.qcplots as sqc
import synotil.dio as sdio
job__nearAUG = sdio.job__nearAUG

import synotil.norm as snorm
sjob= pyutil.sys.modules[__name__]


def figs__peakBW(peakFile,
                 bwFiles,
                 outerRadius = 500,
                 innerRadius=50, 
                 stepSize=10,
                 center_summit = 0,
                 outIndexFunc=None,
                 outIndex=pyutil.basename,
                 ylim = [0,None],
                 NCORE=4,
                 squareSize=(0.2,0.01),
                 detailByChip = None,
                 detailByGene = 0, ### runtime linear in number of genes
                 name = None,
                 **kwargs):
    ### legacy
    if outIndexFunc is not None:
        assert outIndex is None
        outIndex= outIndexFunc
    else:
        assert outIndex is not None
        pass
    #### get the data ready
    bwTracks = sdio.extract_bigwig_multiple(bedFile=peakFile,fnames=bwFiles,
                                            radius=outerRadius,NCORE=NCORE,
                                             stepSize=stepSize,
                                             outIndex=outIndex,
                                            center_summit = center_summit,
                                            **kwargs)
    if detailByChip is None:
        if len(bwTracks) <= 100:
            detailByChip = 1
        else:
            detailByChip = 0
            
    if name is None:
        name = pyutil.getBname(peakFile)
        
    poss = bwTracks.columns.levels[1]
    innerPos = poss[abs(poss) <= innerRadius]

    bwAvg = pyutil.colGroupMean(bwTracks.reindex(columns=innerPos,level=1))
    bwAvg = scount.countMatrix(bwAvg).apply(pyutil.log2p1)

    ##### plotting
#     fig,axs = plt.subplots(1,2,figsize=[12,6])
#     figs = {}
    figs = pyutil.collections.OrderedDict()
    #########
#     plt.figure(figsize=[6,4])
    fig,axs = plt.subplots(1,3,figsize=[18,6])

    ax =axs[0]
    bwAvg.boxplot(rot='vertical',ax=ax)
    ax.set_ylabel('log2(peak intensity)')
    ax.set_ylim(ylim)

    ax =axs[1]
    sqc.qc_pileUp(bwTracks,ax=ax,axLeg=axs[2])
    figs['pileUp-%s' % name ]  = plt.gcf()
    ax.set_ylim(ylim)
    
    #########
    bwAvg.heatmap(figsize=[20,10])
    figs['avgHeatmap-%s' % name] =plt.gcf()
    
    pos = bwTracks.columns.get_level_values('pos')
    cname = 'binding'
    #########
    if detailByChip:
        for key,bwTrack in bwTracks.groupby(axis=1,level=0):
            dfc = bwTrack
            pyvis.heatmap(dfc,
                          transpose=0,
                          squareSize=(0.025,0.2),
                          ytick = dfc.index,
                          xlab = 'distance to %s' % key,
                         vlim=ylim,
                         cname=cname)
            
            ax = plt.gca()
            xticks = pos[ax.get_xticks().astype(int)[:-1]]
            ax.set_xticklabels(xticks,)

#             pyvis.heatmap(bwTrack,
#                           transpose=1,
#                           squareSize=squareSize,
#                           xtick = bwTracks.index,ylab = key)
            figs['detailByChip-%s/%s'%(name,key)] = plt.gcf()
    
    ########
    if detailByGene:

        for key,bwTrack in bwTracks.groupby(axis=0,level=0):
            dfc = bwTrack.melt().pivot_table(index='bwFile',
                                                 columns='pos',)
            pyvis.heatmap(dfc,
                          transpose=0,
                          squareSize=(0.025,0.2),
                          ytick = dfc.index,
#                           xtick = dfc.columns,                          
                          xlab = 'distance to %s' % key,
                          vlim=ylim,
                         cname = cname)
            ax = plt.gca()
            xticks = pos[ax.get_xticks().astype(int)[:-1]]
            ax.set_xticklabels(xticks,)
            
            figs['detailByGene-%s/%s'%(name,key)] = plt.gcf()
            


    
    
    return figs,(bwTracks,bwAvg)

panel_kw_dft =  dict(
    figsize=[14, 7],
    show_axa = 1,
    show_axc = 0,
    showGrid = 0,
    width_ratios = [1,14,0.]
#     title= '',
#     height_ratios = [1,3,3,3,1],
    )
def job__render__panelPlot(
    tracks = None,
    clu=None,
    order=None,
    index=None,
    aliasFmt = '{alias}',
    alias = None,
    baseFile=0,
    figsize = None,
    panel_kw = panel_kw_dft,
    how = 'left',
    debug=0,
    extra = {},
    **kwargs
):
    if figsize is not None:
        panel_kw['figsize'] = figsize
    autoAli = alias is None
    if autoAli:
        alias = ''
    if isinstance(clu,basestring):
        alias += pyext.getBname(clu)
        clu = pyutil.readData(clu,baseFile=baseFile)
    if isinstance(order,basestring):
        alias += pyext.getBname(order)
        order = pyutil.readData(order,baseFile=baseFile)
    if isinstance(tracks,basestring):
        alias += pyext.getBname(tracks)
        tracks = pyutil.readData(tracks,baseFile=baseFile)
        tracks = list(tracks)
    if isinstance(panel_kw,basestring):
        alias += pyext.getBname(panel_kw)
        panel_kw = pyutil.read__buffer( panel_kw,
                                       ext='json',typ='rec',
                                       guess_index=0).to_dict()
    if order is not None:
        clu = order.get(['clu'])
    else:
        assert clu is not None
        order = pd.DataFrame(clu)
        
    if isinstance(index,basestring):
        alias += pyutil.sanitise_query(index)
        locals().update(extra)
        index = eval(index)

    cluTrack = spanel.fixCluster(clu.get(['clu']))
    alias = aliasFmt.format(**locals())
#     cluFile_clean = 'clean_%s.csv' % alias
#     cluc.to_csv(cluFile_clean)    
    tracks = pyext.list__realise(tracks, locals())
    ##### Output heatmap
    pp = spanel.panelPlot(tracks,**panel_kw)
    pp.compile(how=how,
               index=index,
               **kwargs
              )
    pp.compile(order=order)
#     if debug:
#         return pp
    if debug:
        return pp
    fig = pp.render();    
    return (alias,fig)

def job__rawFile__combine(dfc,silent=0):
    '''take a dataframe that identifies raw files and
    convert it to combined fastq.gz ready for downstream.
'''
    dfc = sdio.rawFile__validateChunk(dfc)
    dfcCombined = dfc.groupby('fnameCombined',
                              as_index=False
                             ).apply(sdio.rawFile__combineChunk,
                                            silent=0)
    return dfcCombined


def worker__md5sum(row,column='fname'):
    checkSum_col = '%s_md5sum' % column
    dct= row
    val = dct[column]
    res = pyutil.getMD5sum(val)
    dct[checkSum_col] = res
    return dct

def job__md5sum(dfc,column='fname',n_cpu=1, **kwargs):
    worker = pyutil.functools.partial(worker__md5sum,
                                     column=column)
    it = (x.__dict__ for x in dfc.itertuples())
    res = pyutil.mp_map(worker, it ,n_cpu = n_cpu, **kwargs)
#     res = pd.DataFrame(res)
    return res

def fig__fluffProfile(interval,tracks,ofname = None,
                      annotation= None,
                      scaleOpt = None, 
                      fragmentSize=0,
                      labels = None,
                        silent = 0):
#     annotation = BED12
    # cmd = 'fluff profile'
    trackFlat = u' '.join(tracks)
    if scaleOpt is None:
        scaleOpt = ' -s 1:%d ' % (len(tracks))
    if ofname is None:
        ofname = interval + '.svg'
    cmd =  ''
    cmd += ' fluff profile '
    cmd +=  scaleOpt
    cmd += ' -f {fragmentSize} '
    if annotation is not None:
        cmd += ' -a {annotation} '
    if labels is not None:
        labelFlat = u' '.join(labels)
        cmd += ' -l {labelFlat}'
    cmd += ' -o {ofname} -i {interval} -d {trackFlat} '
    cmd += ' -b white '
    cmd += ' 2>&1 '
    cmd = cmd.format(**locals())
    res=  pyutil.shellexec(cmd,silent=silent)
    return ofname

fig__fluffyProfile = fig__fluffProfile ###legacy


def job__cluster__hpm(tdf,name='test0', 
                      K =40 ,meanNorm=1,
                      threshold=0.,batchSize=500,n_iter=3000,
                      silent=0,
                      NCORE=4,
                      randomState=0,
                      alpha = None,
                      weighted = True,
                     ):
    import pymisca.tensorflow_extra_.hyper_plane_mixture as hpm
    hpm.tf.set_random_seed(randomState)
    np.random.seed(randomState)
    mdl =  hpm.main(K=K,NCORE=NCORE,name=name,meanNorm=meanNorm,threshold=threshold,
                   weighted=weighted,alpha = alpha)
    if batchSize==0 or batchSize is None:
        batchMaker = None
#         batchMaker = hpm.pytfu.batchMaker__random(batchSize=batchSize)
    else:
         batchMaker = hpm.pytfu.batchMaker__random(batchSize=batchSize)
    res = mdl.fit(tdf,
            batchMaker = batchMaker,
            n_iter=n_iter,autoStop=0,
    )
    if not silent:
#         import matplotlib.pyplot as plt
        plt.plot(res)
    cdict = pymod.cache__model4data(mdl,tdf) 
#     assert 0
    mdl.post.__dict__.update(cdict)
    np.save('params.npy', mdl.params)
    res = mdl.params
    res['mdl'] = mdl
    return pyutil.util_obj(**res)

def job__cluster__vmf(tdf,
                      K = 30, init_method='kmeans',
                      weighted= True,
                      n_iter=3000,randomState=None,
                    nStart = 15,
                      min_iters = 50,
                    verbose=1,
                    callback= None,
                    silent=0,
                      sample_weights= 'sd',
):
    import pymisca.model_collection.mixture_vmf as mod    
    np.random.seed(randomState)
    mdl = mod.MixtureVMF(K=K,init_method=init_method,
                         weighted = weighted,
                        )

    histLoss = mdl.fit(tdf,
                       verbose=verbose,
                       callback=callback,
                       nStart=nStart,
                       n_iter = n_iter,
                       min_iters = min_iters,
                       sample_weights =sample_weights,
                      )
    histLoss = -histLoss
    
    if not silent:
#         import matplotlib.pyplot as plt
        plt.plot(histLoss)    
    cdict = pymod.cache__model4data(mdl,tdf) 
    cdict.update(mdl.params)

    np.save('params.npy', cdict)

    cdict['mdl'] = mdl
    return pyutil.util_obj(**cdict)

def job__combinePeak(bwCurr,
                     featSummit='/home/feng/ref/Arabidopsis_thaliana_TAIR10/annotation/genes.gtf.cds.summit',
                     GSIZE='/home/feng/ref/Arabidopsis_thaliana_TAIR10/genome.sizes',
                    CUTOFF=4000,
                     head = 1000,
                    alias = 'testPeaks',
                    center_summit = 1,):

    bwCurr = bwCurr.dropna(subset=['npkFile'])
    bwCurr['npkFileLine'] = bwCurr.eval("npkFile.map(@pyutil.lineCount)")
    print (bwCurr[['bname','npkFileLine']].sort_values('bname'))

    dfs = zip(bwCurr.index,map(sdio.extract_peak,bwCurr.npkFile))
    dfs = dict(dfs)

    dfs = {k:v.sort_values('score',ascending=False) for k,v in dfs.items()}
    bed = pd.concat([df.head(head) for df in dfs.values()],axis=0,sort=False)
    bed = bed.dropna(axis=1)

    ofname = pyutil.to_tsv( bed,'%s__combined.bed' % alias)
    if (len(dfs) > 1) and center_summit :
        ofname = sdio.bed__merge(ofname,silent=0)
    # pyutil.print
    print(ofname,pyutil.lineCount(ofname))

    bedFile= ofname
    
#     peakSummit = sdio.bed__summit(bedFile)
    peakSummit = sdio.npk_expandSummit(fname=bedFile,radius=1,
                                       center_summit=center_summit)
    peak2geneFile = sdio.job__nearAUG(featSummit=featSummit,
                                       peakSummit=peakSummit,
                                       GSIZE=GSIZE,peakWid=1,
                                       CUTOFF=CUTOFF)

    pyutil.fileDict__save('files.json',
                          d=locals(),
                          keys= ['bedFile','peakSummit','peak2geneFile'])
    return 'files.json'


def job__chipTargPaired(
    bwCurr  = None,
    bwMeta = None, control= None, treatment = None,
    xlab = None, ylab = None,
    name = None,
#     bwMeta,
    NCORE = 2,
    params__peakBW = None,
    CUTOFF_FC = 3.0,
    CUTOFF_CHIPDIFF = 0.7,
    innerRadius = 100,
):
    figs = pyutil.collections.OrderedDict()

    
    if control is not None and treatment is not None:
        xlab,ylab = control,treatment
    if xlab is None or ylab is None:
        xlab,ylab = bwCurr.index
    elif bwCurr is None:
        bwCurr = bwMeta.reindex([xlab,ylab])
        
    if params__peakBW is None:
        
        params__peakBW = dict(
            outerRadius=500,
            innerRadius=innerRadius,
            NCORE = NCORE,
            outIndex=bwCurr.header,
        #     detailByCHIP = 0,
        )
    params__peakBW['innerRadius'] = innerRadius        
        
    if name is None:
        name = '{xlab}-{ylab}'.format(**locals())
#     bwCurr = bwMeta
#     bwCurr = bwCurr.loc[[xlab,ylab]]

#     bwCurr.npkFile

    dfs = map(sdio.extract_peak, bwCurr.npkFile,)
    
    fig,ax = plt.subplots(1,1,figsize=[7,7])
#     ax = plt.gca()
    for df in dfs:
        df['per_FC']  = pyutil.dist2ppf(df.FC)
        df.plot.scatter('per_FC','FC',ax=ax)

    fnames = [ pyutil.queryCopy(infile=fname,
                                query='FC>%.3f'%CUTOFF_FC ,
                                reader=sdio.extract_peak,inplace=False) 
        for fname in bwCurr.npkFile]
    # dfs[1]

    peakFlat = ' '.join(fnames)
    ofname = '%s-combined.bed' % ('-'.join(bwCurr.index))
    pyutil.shellexec('cat {peakFlat}>{ofname}'.format(**locals()))
    ofname = sdio.npk_expandSummit(fname=ofname, radius=1)

    pyutil.lineCount(ofname)
    peakFileOrig = peakFile = ofname


    res = sjob.figs__peakBW(
        peakFile = peakFile,
        bwFiles = bwCurr.RPKMFile,
        name = name,
        **params__peakBW
    )
    figs.update(res[0])

    bwTrack,bwAvg = res[1]
    bwAvg.columns = bwCurr.index

    xs,ys = bwAvg[[xlab,ylab]].values.T
#     clu = None
    query = ' val_{ylab} - val_{xlab} > {CUTOFF_CHIPDIFF} '.format(**locals())
    qsans = pyutil.sanitise_query(query)
    
#     peakIndex = pyutil.df__pad(bwAvg).query(query).index
    clu = pd.DataFrame(pyutil.df__pad(bwAvg))
    peakIndex = clu.query(query).index
    clu['clu']= clu.eval('index in @peakIndex')
    
    pyvis.qc_2var(xs,ys,clu=clu.clu,xlab=xlab,ylab=ylab)
    figs['scatterPlot__%s' % name ]= plt.gcf()
    cluFile = ofname = qsans+'.csv'
    clu.to_csv(ofname)
    print (ofname,pyutil.lineCount(ofname))
    peakBase = pyutil.getBname(peakFile)
    ofname =  '{peakBase}-{qsans}.bed'.format(**locals())
    peakFile = pyutil.to_tsv(
        sdio.extract_peak(peakFile).set_index('acc',drop=0).reindex(peakIndex),
                           ofname)
    pyutil.shellexec('mkdir -p output/')
    pyutil.file__link(ofname,'output/%s.bed' % name,force=True)
    
#     peakFile = pyutil.queryCopy(peakFile,
#                                 query='acc in @peakIndex',
#                                 reader=sdio.extract_peak, 
#                                 peakIndex=peakIndex,
#                                )
#     peakFile =  '{peakFile}-{qsans}.bed'
#     pyutil.fileDict__main(ofname='FILE.json',
#                          **pyutil.dictFilter(locals(),
#                                              keys=['cluFile','peakFile',
#                                             'peakFileOrig']
#                                             ))
    
    pyutil.fileDict__save(d=locals(), keys=['cluFile','peakFile',
                                            'peakFileOrig'],fname='FILE.json')
    return figs,clu