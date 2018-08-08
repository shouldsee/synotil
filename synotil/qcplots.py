# import mapplotlib.pyplot as plt
import scipy.cluster.hierarchy as sphier
import scipy.spatial.distance as spdist
import pymisca.vis_util as pyvis
import pymisca.util as pyutil
plt = pyvis.plt

import numpy as np

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
def make_compareClus(tks,stats = None):
    if not isinstance(tks,list):
        tks = [tks]
    def compareClus(clus,tks=tks,stats0=stats, figsize=[14,5],
                    shortName=0,
                    how= 'left',**kwargs):

        cluTracks = map(spanel.fixCluster,clus)
        
        tracks = cluTracks + tks

        pp =spanel.panelPlot(
            tracks
        )

#         stats = pd.concat(clus + stats,axis=1)
        L = len(clus)
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


def qc_Avg(C, silent=1,axis=1,nMax = 150, **kwargs):
#     if axs is None:
#         if not silent:
#             fig,axs= plt.subplots(1,3,figsize=[14,3])
    C = np.array(C)
    assert C.shape[axis]<nMax
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