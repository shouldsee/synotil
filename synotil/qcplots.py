# import mapplotlib.pyplot as plt
import scipy.cluster.hierarchy as sphier
import scipy.spatial.distance as spdist
import pymisca.vis_util as pyvis
import pymisca.util as pyutil
plt = pyvis.plt
import synotil.modelRoutine as smod
import synotil.CountMatrix as scount

import numpy as np

def qc_PCA(df,xi=0, yi = 1, xlab=None,ylab = None,**kwargs):
    res = smod.fit_PCA(df)
    xs,ys = res.trans_data.T[[xi,yi]]
    if xlab is None:
        xlab = 'PC%d'%xi
    if ylab is None:
        ylab = 'PC%d'%yi
    axs = pyvis.qc_2var(xs,ys,xlab=xlab,ylab=ylab,**kwargs)
    return res,axs

def qc_dist(D,vlim=None,silent=1,axs = None,
            cutoff=0.75,reorder = 1,method='average',
            distance_sort='descending',level=5):
    '''Plot a distance matrix. Perform hierarchical clustering if reorder == True
'''
    if D.ndim == 1:
        D = spdist.squareform(D)
    if reorder:
        D_1d = spdist.squareform(D,checks=0)
        K = sphier.linkage(D_1d,method=method)
        
        if isinstance( cutoff,float):
            par = {'height':cutoff}
        elif isinstance( cutoff,int):
            par = {'n_clusters':cutoff}
            
        dendo = sphier.dendrogram(K,no_plot=1,distance_sort=distance_sort,)
        od = dendo['leaves'][::-1]
        clu = sphier.cut_tree(K,**par).ravel()
        
        #### Sort cluster as dendogram shows
        mapper = {k:i  for i,(k,_) in enumerate(pyutil.itertools.groupby(clu[od]))}        
        clu = map(mapper.get, clu)
#         od = np.argsort(clu)
        D = D[od][:,od]    
    else:
        clu = np.arange(len(D))
    
    if not silent:
        if axs is None:
            fig,axs= plt.subplots(1,3,figsize=[14,4]);
            
        axs=axs.ravel(); i = -1
        
        i+=1;plt.sca(axs[i]);ax=axs[i]
        pyvis.histoLine(D.ravel(),xlim=vlim)
        plt.grid()
        plt.xlabel('distance')
        plt.ylabel('count')


        i+=1;plt.sca(axs[i]);ax=axs[i]
        if reorder:
            sphier.dendrogram(K,ax=ax,orientation='left',distance_sort=distance_sort,
                             truncate_mode='level',p=level,
                             no_labels=1);    
            if isinstance( cutoff,float):         
                plt.vlines(cutoff,*ax.get_ylim(),linestyles='--')

        i+=1;plt.sca(axs[i]);ax=axs[i]
        pyvis.heatmap(D,cname='distance',ax=ax,vlim=vlim)
        pyvis.hide_axis(ax)

    return clu,axs

import PanelPlot as spanel
import pandas as pd
def make_compareClus(tks,stats = []):
    if not isinstance(tks,list):
        tks = [tks]
    def compareClus(clus,tks=tks,stats0=stats, figsize=[14,5],
                    shortName=0,L = None,
                    how= 'left',**kwargs):

#         stats0 = 
        cluTracks = map(spanel.fixCluster,clus)
        
        tracks = cluTracks + tks

        pp =spanel.panelPlot(
            tracks
        )

#         stats = pd.concat(clus + stats,axis=1)
        L = len(clus) if L is None else L
#         L = len(stats.columns)
        figs = []
        for coli in range(L):
            ind = pyutil.escalIndex(L,coli)
            ind = [coli]
            statsc = pd.concat([clus[coli]] + stats0,axis=1)
#             statsc = stats.iloc[:,ind]
        #     statsc = stats.iloc[]
        #     statsc = stats.columns
            pp = spanel.panelPlot(tracks,)
            fig = pp.render(order=statsc,shortName=shortName
                                                  ,figsize=figsize,how=how)
#             fig = plt.gcf()
            figs +=[fig]
        #     pp = spanel.panelPlot(tracks,).render(order=stats.iloc[:,::-1])
        return tracks,figs
    return compareClus


# import scipy.cluster.hierarchy as sphier

import sklearn.metrics as skmet
def worker_cut_tree(nClu,Z=None):
    print nClu,
    if nClu <= 1:
        res = None
    else:
        res = sphier.cut_tree( Z,n_clusters=nClu)
    return res
def worker_silhouette(clu,Ds=None):
    if clu is None:
        res = None
    else:
        nClu = clu.max() + 1
        print nClu,
    #     clu = sphier.cut_tree( Z,n_clusters=nClu)
        S =skmet.silhouette_score(Ds, clu,metric='precomputed')
        res = (nClu,S)
    return res

def qc_silhouette(D,
#                   nClu = 40                 
method = 'complete',
nClu = 40,
NCORE=10,
silent=1,
axs= None
):
    Ds = D

    if D.ndim == 2:
        D = spdist.squareform(D,checks=0)
    if Ds.ndim != 2:
        Ds = spdist.squareform(Ds,checks=0)

    # D = D_mse
    Z = sphier.linkage(D,method= method,)
    lst = []
    clus = []
    nClus = range(0,nClu)
    worker = pyutil.functools.partial(
        worker_cut_tree,
        Z=Z)    
    clus = pyutil.mp_map(worker, 
              nClus,n_cpu=NCORE)
    
    worker = pyutil.functools.partial(
        worker_silhouette,
        Ds = Ds)
#     lst = []
    lst = pyutil.mp_map(worker,
                       clus,
#                        n_cpu=min(NCORE,4)
                       n_cpu=1
                       )

    lst = [x for x in lst if x is not None]
    shind = lst
    shind= np.array(shind)
    if not silent:
        if axs is None:
            fig,axs= plt.subplots(1,1,figsize=[6,4]);
            axs = [axs]
        ax = axs[0]
        X,Y = shind.T
        ax.plot(X,Y)
        # pyvis.abline()
        ax.set_xlim(left=0)
        wid = 3; rid = wid//2
        movavg = map(np.mean,pyutil.window(Y,n=wid,step=1,),)        
        ax.plot(X[rid:-rid], movavg,'x--',
               label = 'Moving average')
        ax.grid(1)
        ax.legend()

    return shind

def qc_cumAvg(X,axis=0,
             silent=1,
             axs=None,
             nMax = int(1E6),
              ybin = None,
#               axiS = [2,None,None]
             ):
    '''Calcuate average and standard error for bootstrapped statistics.
'''
# axis = 0
#     X = egList[:5]
    X = np.array(X)
    X = np.moveaxis(X,axis,0)
    # L = np.shape(X)[axis]
    X = np.reshape(X,(len(X),-1))
    if len(X.T) > nMax:
        ind = np.random.randint(0,len(X.T),nMax)
        X = X[:,ind]

    L = len(X)
    Lx = (1+np.arange(L))[:,None]

    Ex = np.cumsum(X,axis=axis)/Lx
    Ex2 = np.cumsum(np.square(X),axis=axis)/Lx

    M = Ex
    VAR = (Ex2 - np.square(Ex) )
    SD = np.sqrt(VAR)
    SE = SD/np.sqrt(Lx)
    CV = SE/M; CV = abs(CV)
    Lx = np.broadcast_to(Lx,SE.shape) 
    if not silent:
        if axs is None:
            fig,axs= plt.subplots(1,3,figsize=[14,4]);
            axs=[None,None,axs[0]]
        X,Y = Lx,CV
        ybin = np.linspace(*pyutil.span(Y,99.),num=80) if ybin is None else ybin
        xbin = np.arange(0.5,L+1,1)
        pyvis.qc_2var(X,Y,
#                       axs=axs,
                      xbin=xbin,ybin = ybin,
                      ylab=r'$\left| {StdErr} / {Mean} \right |$',
                     xlab='sample size',
                      axs=axs,
                     )
    return (M,SE,CV,Lx),axs


def qc_Avg(C, silent=1,axis=1,
#            nMax = 150, ### depracated size check
           **kwargs):
    
#     if axs is None:
#         if not silent:
#             fig,axs= plt.subplots(1,3,figsize=[14,3])
    C = np.array(C)
#     assert C.shape[axis]<nMax
    MEAN = C.mean(axis=axis,)
    STD = C.std(axis=axis,)
    # plt.hist(X) def parseBedmap
    # plt.hist(X[1])
#     X = MEAN[None,:]
    X = MEAN
    MIN,MAX = X.min(),np.percentile(X,99)
    BINS = np.linspace(MIN,MAX,100)
    CV = STD/MEAN
    if not silent:
        xs,ys = MEAN,STD
        axs = pyvis.qc_2var(xs,ys,xlab='$E(X)$',ylab='$Std(X)$',
                            **kwargs)
    else:
        axs = []
    return (MEAN,STD,CV),axs
def qcAvg(*args,**kwargs):
    '''Legacy support''' 
    return qc_Avg(*args,**kwargs)

def qc_meanVar( C, clu, axs=None,xlim=None,ylim=None,silent=0):
    ''' C of shape (n_gene, n_condition)
    Points colored by cluster
    '''
    clu = np.array(clu).ravel()
    nClu = np.max(clu)+1
    if axs is None:
        fig,axs= plt.subplots(1,3,figsize=[14,3])
    for ci in range(nClu):
        idx = np.where(clu==ci)[0]
        CC = C[idx,:]
        STAT,axs = qcAvg(CC, axs=axs, silent=silent)
    ax = axs[1]
#     ax.set_alpha(0.5)
    MEAN = C.mean(axis=1,keepdims=1).squeeze()
    STD = C.std(axis=1,keepdims = 1).squeeze()
    
    xlim = xlim if xlim is not None else np.span(MEAN,99.9)
    ylim = ylim if ylim is not None else np.span(STD,99.9)
    
    xlim = np.span(MEAN,99.)
    ylim = np.span(STD,99.)
#     for ax in axs:
    ax = axs[0]
    ax.set_xlim(xlim);#ax.set_ylim(ylim)
    ax = axs[1]
    ax.set_xlim(xlim);ax.set_ylim(ylim)
    ax = axs[2]
    ax.set_xlim(xlim);ax.set_ylim(ylim)
    return ((MEAN,STD,STD/MEAN),axs)

# def qc_pileUp(bwt,ax=None,silent =0,
#              sigMax = 20,
#              ):
#     '''bwt: A dataFrame containing the peak regions from a bigWig track
#     sigMax: throw away peaks with average above this maximum
#     sigMax: clip the signal at this maximum 
# '''
#     if not isinstance(bwt,scount.countMatrix):
#         bwt = scount.countMatrix(bwt)
        
#     if bwt.index.duplicated().any():
#         bwt =bwt.reset_index(drop=1,inplace=False)
#     bwt.qc_Avg()
#     index = bwt.summary.query('M<%d'%sigMax).index
#     bwtc = bwt.reindex(index)
#     xs = bwt.columns

#     (M,SD,CV), _ = qc_Avg(bwtc.T,nMax=1E100)
#     SE = SD/len(M)**0.5

#     if not silent:
#         if ax is None:
#             fig,axs= plt.subplots(1,3,figsize=[14,3])    
#             i = -1
#             i += 1;ax=axs[i]; plt.sca(ax)
#         else:
#             axs = None
#         plt.plot(M,'b')
#         plt.plot(M+2*SE, 'b--', alpha=0.5)
#         plt.plot(M-2*SE, 'b--', alpha=0.5)
# #         plt.title(pyutil.basename(fname))
#         plt.grid(1)
        
#     return (M, SE), ax

def qc_pileUp(bwTracks,ax=None,errKey = 'SE',labels=None):
    '''bwTracks is dataFrame with multi-indexed columns (trackname,position), and single-indexed
    rows (peakName,)
'''
    if ax is None:
        fig,axs = plt.subplots(1,1,figsize=[10,6])
        ax =axs
#         fig,ax = plt.subplots(1,1,figsize=[7,7])
    dfc = scount.countMatrix(bwTracks.T)
    
    dfc.qc_Avg()
    summ = dfc.summary.reset_index()
    pileUp = summ.pivot_table(columns='bwFile',index='pos',values=['M','SD'],
                             )
    cols = summ['bwFile'].drop_duplicates()
    pileUp = pileUp.reindex(columns=cols,level=1)
    if labels is not None:
        pileUp.columns = pileUp.columns.set_levels(labels,level='bwFile')

    pyvis.linePlot4DF(pileUp['M'].T,which='plot',ax=ax)
    M,SD = pileUp.M,pileUp.SD
    SE = SD / len(dfc)**0.5
    if errKey is not None:
        err= locals()[errKey]
        ax = pyvis.linePlot4DF(df= M.T - err.T,
                               y2= M.T + err.T,which='fill_between',alpha=0.5,ax=ax)
    ax.legend()
    ax.set_xlabel('relative to peak (bp)')

    return ax


# ?import pymisca
import synotil.CountMatrix as scount
def qc_lineplots(dfcc, ylab = 'average RPGC',
                 xlab = 'bp to the start of peak',
                 L =None):
    if not isinstance(dfcc, scount.countMatrix):
        dfcc = scount.countMatrix(dfcc)
    L = dfcc.param.__dict__.get('n_peak',None) if L is None else L
    
    dfcc.heatmap(cname=ylab,vlim=[0.,6.])

    fig,axs = plt.subplots(1,2,figsize=[12,4])
    ax = axs[0]

    dfcc.plot.line(legend=0,ax=ax)
    ax.set_ylim(0,None)
    ax.grid(1)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)

    ax = axs[1]
    dfcc.plot.line(ax=ax)

    plt.suptitle('N=%d'% L )
    return axs


import synotil.dio as sdio
def qc_narrowPeak(qfile,
                  cutoff = 0.98,
                  ax = None,
                  silent= 1,
                  keyFile = None,
                  ofname = None,
                  cutoff_key = 'per_FC',
#                   cutoff = {'per_FC':0.98}
                 ):
    '''
    Visualise the fold-change distribution and do cutoff
'''
    f=open(qfile)
    fline = f.readline()
    f.close()
    
    if fline.split('\t')[0] == 'chrom':
        qres = pyutil.readData( qfile, guess_index=0)
        pass
    else:
        qres = sdio.extract_peak(qfile)
        
    qres['per_FC'] = pyutil.dist2ppf(qres.FC)
    qres['per_score'] = pyutil.dist2ppf(qres.score)
    
    dfc = qres.query('%s > %.3f'%(cutoff_key,cutoff))
    
    ofname= '%s_chipTarg.tsv'%pyutil.basename(qfile)
    
    dfc.reset_index(drop=1).to_csv(ofname,sep='\t',index=None,header=None)
    print (qres.shape,dfc.shape)
#     print (dfc.shape)
    
    if keyFile is not None:
        keyDF = pyutil.readData(keyFile)
        
        dfcc = dfc.set_index('feature_acc',drop=0)
        dfcc = dfcc.loc[~dfcc.index.duplicated(),]
        
        keyTarg= pd.concat([dfcc[['FC']],keyDF],axis=1,join='inner')
        
        pyutil.ipd.display(keyTarg)

    
    if not silent:
        if ax is None:
            fig,axs = plt.subplots(1,2,figsize=[12,4])
            ax = axs[0]
        plt.sca(ax);
    #     cutoff = 0.98
        raw_key = cutoff_key.split('_')[-1]
        ax.plot( qres[cutoff_key], qres[raw_key] ,'x')

        ax.set_xlim(0.5,1.1)
        ax.grid(1)
        ax.vlines(cutoff,*ax.get_ylim())
        ax.set_xlabel('percentile')
        ax.set_ylabel(raw_key)
        title = 'All=%d'%len(qres) + ', keep=%d'%len(dfc)
        ax.set_title(title)
        
    
    # dfcc =qc_ sutil.tidyBd(dfc.set_index('feature_acc',drop=0))
    
    return ofname,ax

import synotil.dio as sdio;
def qc_summitDist(peak1,peak2,GSIZE,query=None,query1=None,query2=None,
                  xlab = None,ylab=None,
                 CUTOFF=600,
                  axs=None,
#                  ax = None,
                 ):
    '''plot the distribution of inter-summit distance between files
'''
    xbin = np.linspace(2,CUTOFF,50)
    if axs is None:
        fig,axs = plt.subplots(1,3,figsize=[16,4])
#         ax = axs[0]
    i = -1;
    

    infiles = [peak1,peak2]
    
    if query is not None:
        query1 = query2 = query
    if query1 is not None and query2 is not None:
        querys = [query1,query2]
        infiles = [
            pyutil.queryCopy(
                reader = sdio.extract_peak,
                infile = infile,
                query = query,
                inplace=True,
            )
            for query,infile in zip(querys,infiles)
            ]
        
    if xlab is None:
        xlab = pyutil.basename(peak1)
    if ylab is None:
        ylab = pyutil.basename(peak2)

    randFiles = [sdio.bed_randomise(fname,GSIZE=GSIZE)
                for fname in infiles]
    peak1,peak2 = infiles

    i += 1; ax = axs[i];plt.sca(ax)
    plotter = ax.hist
    # plotter = pyvis.histoLine 
    common ={'bins':xbin,
            'density':0,
            'alpha':1.0,
             'histtype':'step',
            }

    df = df__inter = sdio.summitDist(peak1,peak2,GSIZE=GSIZE,CUTOFF=CUTOFF)
    lab = '%s__%s' % (xlab,ylab)
    lab += ',N=%d'%len(df)
    plotter(df.distance,label=lab,**common)
    
    
    
    df = sdio.summitDist(peak1,peak1,GSIZE=GSIZE,CUTOFF=CUTOFF)
    lab = '%s__%s' % (xlab,xlab)
    lab += ',N=%d'%len(df)
    plotter(df.distance,label=lab,**common)
    
    
    
    df = sdio.summitDist(peak2,peak2,GSIZE=GSIZE,CUTOFF=CUTOFF)
    lab = '%s__%s' % (ylab,ylab)
    lab += ',N=%d'%len(df)
    plotter(df.distance,label=lab,**common)

    df = sdio.summitDist(peak1,randFiles[-1],GSIZE=GSIZE,CUTOFF=CUTOFF)
    lab = '%s__randomise-%s' % (xlab,ylab)
    lab += ',N=%d'%len(df)
    plotter(df.distance,label=lab,**common)

    df = sdio.summitDist(peak2,randFiles[0],GSIZE=GSIZE,CUTOFF=CUTOFF)
    lab = 'randomise-%s__%s' % (xlab,ylab)
    lab += ',N=%d'%len(df)
    plotter(df.distance,label=lab,**common)
    
    title = 'query1={query1},query2={query2}'.format(**locals())
    ax.set_title(title)
    ax.set_xlabel('inter-summit distance (bp)')
    ax.legend()
    ax.grid(1)
    
    df = df_inter = df__inter.query('distance>1')
    L1and2 = len(df.acc.unique())
    L1not2 = pyutil.lineCount(peak1) - L1and2

    L2and1 = len(df.feat_acc.unique())
    L2not1 = pyutil.lineCount(peak2)  - L2and1

    Lall  = (L1and2 + L2and1)//2
    i += 1; ax = axs[i];plt.sca(ax)
    indVenn = pyvis.qc_index(subsets=(L1not2,L2not1,Lall),silent=0,xlab=xlab,ylab=ylab,ax=ax)[0]
    return df_inter,indVenn, axs


