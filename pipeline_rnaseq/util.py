
# coding: utf-8

# In[ ]:


####### Utilities for Visualising RNASeq. Author: Feng Geng (fg368@cam.ac.uk)
#### Written for BrachyPhoton at SLCU


# In[46]:


if __name__=='__main__':
    get_ipython().system(u'jupyter nbconvert --to python tmp.ipynb')
# !python compile_meta.ipynb && echo '[succ]'


# In[88]:


INDEX = '/media/pw_synology3/BrachyPhoton/raw/index'
OUTDIR = '/media/pw_synology3/BrachyPhoton/Mapped_data'
# ! head {INDEX}

import os,re,sys
import pandas as pd
import numpy as np
import pymisca.vis_util as pyvis
plt=pyvis.plt
import pymisca.util as pyutil
np = pyutil.np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


import sys
modCurr = sys.modules[__name__]



def readLines(fname):
    with open(fname,'r') as f:
        lines = [l.rstrip('\n') for l in f.readlines()]
    return lines


try:
    INDIRs = [os.path.join(OUTDIR,l) for l in readLines(INDEX)]

    ### remove the blacklisted samples
    INDIRs = [ i for i in INDIRs if not bool(re.search("169R_12", i))]
except Exception as e:
    print '[FAIL] to process index file:%s, due to %s'%(INDEX,e)


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map
    Source: https://gist.github.com/jakevdp/91077b0cae40f8f8244a
    """

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    if base_cmap is None:
        base = plt.get_cmap()
    else:
        base = plt.cm.get_cmap(base_cmap)
    
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

    
def histTPM(df,COL='TPM'):
    fig,axs= plt.subplots(1,2,figsize=[14,3])
    vals = df[COL]
    plt.sca(axs[0])
    plt.hist(np.log1p(vals),30,log=0)
    plt.vlines(1,0,6000)
    plt.grid()
#     plt.show()

    plt.sca(axs[1])
    plt.hist(vals,30,log=1)
    plt.grid()
#     plt.show()
    
def timePCA(M,ZTime):
#     COL = meta[colorName]
    COL_RGB,(COL_LAB,COL_LST) = ser2col(pd.Series(ZTime))
    
    fig,axs= plt.subplots(1,2,figsize=[14,4])

    plt.sca(axs[0])    
    labs = np.arange(len(M))
    x = [int(x.lstrip('ZT')) for x in COL]
    y = M[:,0]
    plt.xlabel('time')
    plt.ylabel('PC1')
    l = plt.scatter(x,y,c=COL_RGB,)
    for i,(xx,yy,lab) in enumerate(zip(x,y,labs)):
        plt.annotate(lab,xy=(xx,yy), )
#     plt.legend()
    plt.grid()
    
    plt.sca(axs[1])    
    y = M[:,1]
    plt.xlabel('time')
    plt.ylabel('PC1')
    l = plt.scatter(x,y,c=COL_RGB,)
    for i,(xx,yy,lab) in enumerate(zip(x,y,labs)):
        plt.annotate(lab,xy=(xx,yy), )
#     plt.legend()
    plt.grid()
    return fig    
def histoLine(xs,BINS=None,log= 0,**kwargs):
    ys,edg = np.histogram(xs,BINS)
    ct = (edg[1:] + edg[:-1])/2
    if log:
        ys = np.log1p(ys)
    else:
        pass
    l =plt.plot(ct,ys,**kwargs)
    return l
def qc_Avg(C,axs = None,silent=1):
    if axs is None:
        if not silent:
            fig,axs= plt.subplots(1,3,figsize=[14,3])
    assert C.shape[1]<100
    MEAN = C.mean(axis=1,keepdims=1).squeeze()
    STD = C.std(axis=1,keepdims = 1).squeeze()
    # plt.hist(X) 
    # plt.hist(X[1])
    X = MEAN[None,:]

    MIN,MAX = X.min(),np.percentile(X,99)
    BINS = np.linspace(MIN,MAX,100)
    CV = STD/MEAN
    if not silent:
        BX = np.linspace(*np.span(MEAN,99.9),num=30)
        BY = np.linspace(*np.span(STD,99.9),num=50)
        xlim = np.span(BX)
        ylim = np.span(BY)
        plt.sca(axs[0])
        for i in range(1):
            histoLine(X[i],BINS,alpha=0.4)

        plt.grid(1)
        plt.xlabel('$E(X)$')

        plt.sca(axs[1])

        i = 10000
        idx = np.random.randint(len(MEAN),size=i)
        x,y = MEAN[idx],STD[idx]
        # plt.hist2d(MEAN.squeeze(),STD.squeeze(),(50,50))
        # plt.show()
        plt.plot(x,y,'.')
        plt.grid(1)
        abline()
        plt.xlabel('$E(X)$')
        plt.ylabel('$Std(X)$')
        plt.xlim(xlim)
        plt.ylim(ylim)

        plt.sca(axs[2])

        ct,BX,BY = np.histogram2d(MEAN,STD,(BX,BY))
        plt.pcolormesh(BX,BY,np.log1p(ct).T,)
    #     plt.gca().matshow(np.log1p(ct),aspect='auto')
        plt.xlabel('$E(X)$')
        plt.ylabel('$Std(X)$')
    return (MEAN,STD,CV),axs
def qcAvg(*args,**kwargs):
    '''Legacy support''' 
    return qc_Avg(*args,**kwargs)

def qc_meanVar(C,clu,axs=None):
    ''' C of shape (n_gene, n_condition)
    Points colored by cluster
    '''
    clu = np.array(clu)
    nClu = np.max(clu)+1
    if axs is None:
        fig,axs= plt.subplots(1,3,figsize=[14,3])
    for ci in range(nClu):
        idx = np.where(clu==ci)[0]
        CC = C[idx,:]
        STAT,axs = qcAvg(CC,axs=axs)
    ax = axs[1]
#     ax.set_alpha(0.5)
    MEAN = C.mean(axis=1,keepdims=1).squeeze()
    STD = C.std(axis=1,keepdims = 1).squeeze()
    xlim = np.span(MEAN,99.)
    ylim = np.span(STD,99.)
    ax.set_xlim(0,5);ax.set_ylim(0,2)
    return ((MEAN,STD,STD/MEAN),axs)
    
def matHist(X,idx=None,XLIM=[0,200],nbin=100):    
    plt.figure(figsize=[12,4])
    if idx is not None:
        X = X[idx]
    MIN,MAX = X.min(),np.percentile(X,99)
    BINS = np.linspace(MIN,MAX,nbin)
    for i in range(len(idx)):
        histoLine(X[i],BINS,alpha=0.4,log=1)
    plt.xlim(XLIM)
    plt.grid()
    
def getCV(xs):
    return np.std(xs)/np.mean(xs)

def getCol(dfs,COLUMN='Coverage'):
    ctraw = [df.get(COLUMN).values for df in dfs]
#     ctraw = map(lambda x:getTPM(x,COL=COLUMN),dfs)
    ctraw = np.array(ctraw)
    return ctraw

def abline(k=1,y0=0):
    '''Add a reference line
    '''
    MIN,MAX=plt.gca().get_xlim()
    f = lambda x: k*x+y0
    plt.plot([MIN,MAX],[f(MIN),f(MAX)],'b--')
    print MIN,MAX
    
def subset(dfs,idx):
    return [df.iloc[idx] for df in dfs]


# In[35]:


#### I/O utility

def combine_csv(fnames):
    geneSet = set()
    dfs = []
    CUTOFF = 6

    FCOL = 'Coverage'
    geneAll = set()
    geneAny = set()
    geneRef = pd.DataFrame()
    # for i,fname in enumerate(fnames[:10]):
    for i,fname in enumerate(fnames):
        if not i%10:
            print 'Reading %s'%fname
        df = pd.read_table(fname)
        allGene = df['Gene ID']
        exprGene = allGene.loc[df[FCOL]>CUTOFF]    

        geneDiff = set(allGene).difference(geneAll)    
        appd = df[allGene.isin(geneDiff)]

        geneRef = geneRef.append(appd)

        geneAll.update(allGene)
        geneAny.update(exprGene)
    #     break
    #     if not geneAll:
    #         geneAll.update(df['Gene ID'])
    #     else:
    #         geneSet.intersection_update(df['Gene ID'])
        dfs.append(df)

    geneRef.loc[:,['Coverage','FPKM','TPM']] = 0
    geneRef.sort_values('Gene ID',inplace=1)
    geneRef = geneRef.reset_index()
    geneSet = geneAny
    geneValid = geneRef[geneRef['Gene ID'].isin(geneSet)] 

    print 'Nmber of Genes before filtering:',len(geneAll)
    print 'Nmber of Genes after filtering:',len(geneSet)
    print 'Surviving rate: ',float(len(geneSet))/len(geneAll)
    
    return dfs,(geneRef,geneValid)
def padWithRef(df,ref):
    df = df.append(ref[~ref['Gene ID'].isin(df['Gene ID'])])
    df.sort_values('Gene ID',inplace=1)
    df = df.reset_index()
    return df
def routine_combineCSV(fnames):
    print '[PROG] Starting to readfile'
    dfs,(geneRef,geneValid) = combine_csv(fnames)
    print '[PROG] Finished to readfile'
    
    print '[PROG] Starting to pad'
    f = pyutil.functools.partial(padWithRef,ref=geneRef)
    lst = pyutil.mp_map(f,dfs,n_cpu=1)

    SHP = np.array([df.shape for df in lst])
    assert np.all(SHP == SHP[0:1]),'Arrays not sharing shape:%s'%SHP
    gids = np.array([df['Gene ID'] for df in lst])
    assert np.all(gids == gids[0:1])
    print '[PROG] Finished padding'
    
    dfs = lst
    dfs = [df.iloc[geneValid.index] for df in dfs]
    return dfs,(geneRef,geneValid.reset_index().drop('index',1))


# In[121]:


def preprocess(C,std=1):
    C = np.log1p(C)
#     C = C[clu==1,:][:,meta_wt_LD.index]
    C = (C-C.mean(axis=1,keepdims=1))
    if std:
        STD= C.std(axis=1,keepdims=1)
        nonCst = (STD!=0).squeeze()
        C[nonCst] = C[nonCst]/STD[nonCst]
    return C


def msgGMM(model = None, train_data = None,name='test',**kwargs):
    mdl = model
    s = '''
Name:{name}
Converged:{cvg}
min_logL: {logL}
(lower-bound of) MEAN logL :{mlogL}'''.format(
        name=name,
         cvg=mdl.converged_,
         logL = mdl.lower_bound_,
         mlogL=mdl.lower_bound_/len(train_data) if not train_data is None else mdl.lower_bound_
                                )
    return s

# def sortLabel(Y,X,#return_pos=0
#              ):
#     ### Sorting by maximum
# #         X = vX
#     X  = X-X.mean(axis=1,keepdims=1)
#     coord = np.arange(X.shape[-1])[None]
#     wt_X = (X == X.max(axis=1,keepdims=1))*coord
# #         wt_X = X * coord
#     cis = list(range(max(Y)+1))
#     pos =  [wt_X[Y == ci,:].mean() for ci in cis]
#     sidx = np.argsort(pos)
# #     if return_pos:
# #         return pos[sidx]
# #     else:
#     pj = dict(zip(sidx,cis))
#     Y = np.vectorize(pj.get)(Y)
#     return Y,np.take(pos,sidx)
def qcGMM(model,train_data,name='Test',valid_data = None,pt=None,axs = None,**kwargs):
    mdl = model
    X   = train_data
    if valid_data is None:
        valid_data = X
    vX = valid_data
    if axs is None:
        fig,axs = plt.subplots(1,2,gridspec_kw={"width_ratios": (.1, .9),
                                               'wspace':0.1,
                                                'top':0.8
    #                                             'hspace':0.5
                                               },
                              figsize=[14,3])
    Y = mdl.predict(X)
    s = mdl.score_samples(X)
    

    Y,pos = sortLabel(Y,X)
    
        
    idx = np.argsort(Y)
    if pt is not None:
        ps = np.percentile(s, pt)
        if pt > 50:
            i2 = s[idx] > ps
        else:
            i2 = s[idx] < ps
        im = vX[idx[i2]].T
    else:
        im = vX[idx].T
        
    ax = axs[0]
    plt.sca(axs[0])        
    ax.hist(s,50);
    if pt is not None:
        ax.vlines(ps,0,2000)
    plt.grid()

#     plt.sca(axs[1])
    ax = axs[1]
    ax.matshow(im,aspect='auto')
    ax.xaxis.tick_bottom()
    ax.set_title(name)


# In[124]:


if __name__=='__main__':
    get_ipython().system(u'jupyter nbconvert --to python tmp.ipynb')
# !python compile_meta.ipynb && echo '[succ]'


# In[ ]:


#### PCA utilities

import pymisca.vis_util as pyvis
def fit_PCA(C,n_components=5,**kwargs):
    import sklearn.decomposition.pca as skpca
    mdl = skpca.PCA(n_components=n_components,**kwargs)
    M = mdl.fit_transform(C)
    
    return {'model':mdl,
            'train_data':C,
            'trans_data':M,}

def quickPCA(trans_data=None,model=None,COL_SER=None,index=None,**kwargs):
    M = trans_data
    mdl = model
    nSample =   len(M)
    assert nSample < 100,'Too many samples in the maxtrix: %d>100'%nSample
    
#     labs = np.arange(nSample)
    labs = index
    vara = labs * 0. if mdl is None else mdl.explained_variance_ratio_

    if not isinstance(COL_SER,pd.Series):
        COL_SER = pd.Series(COL_SER)

    
    COL_RGB,(COL_LAB,COL_LST) = ser2col(COL_SER)
    common = {
             }
    fig,axs= plt.subplots(1,2,figsize=[14,4])

    def pc2d(pi,pj):
        x, y  = M[:,pi],M[:,pj]
        l = plt.scatter(x,y,c=COL_RGB,)
        for _,(xx,yy,lab) in enumerate(zip(x,y,labs)):
            plt.annotate(lab,xy=(xx,yy), )
    #     plt.legend()
        plt.grid()
        plt.xlabel('PC%d(%.1f%%)'%(pi+1,vara[pi]*100))
        plt.ylabel('PC%d(%.1f%%)'%(pj+1,vara[pj]*100))
    plt.sca(axs[0])
    pc2d(0,1)
    plt.sca(axs[1])
    pc2d(2,3)
    recs = [pyvis.mpl.patches.Rectangle((0,0),1,1,fc=c) for c in COL_LST]
    fig.legend(recs,COL_LAB)    
def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map
    Source: https://gist.github.com/jakevdp/91077b0cae40f8f8244a
    """

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    if base_cmap is None:
        base = plt.get_cmap()
    else:
        base = plt.cm.get_cmap(base_cmap)
    
    color_list = base(np.linspace(0, 1, N+1))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N+1)
def ser2col(COL_SER):
    COL_VAL, COL_LAB= COL_SER.factorize()
    NCAT = len(COL_LAB)
#     cmap = plt.get_cmap('jet')
    cmap = discrete_cmap(NCAT,'jet')
#     print (COL_VAL.ravel())
    COL_RGB = cmap(COL_VAL.ravel())
    COL_LST = cmap(range(NCAT))
    return COL_RGB,(COL_LAB,COL_LST)

    
def timePCA(ZTime_int,trans_data=None, model = None,COL_SER=None,index=None,**kwargs):
    mdl = model; M = trans_data
#     COL = meta[colorName]
    if COL_SER is None:
        COL_SER = ZTime_int
    COL_RGB,(COL_LAB,COL_LST) = ser2col(pd.Series(COL_SER))
    
    if index is None:
        labs = np.arange(len(M))
        print '[WARN] index not specified'
    else:
        labs = index
    vara = labs * 0. if mdl is None else mdl.explained_variance_ratio_
        
    nPC = 4
    fig,axs= plt.subplots(1,nPC,figsize=[14,4])
        
    for i in range(nPC):
        pi = i
        plt.sca(axs[i])    
        x = ZTime_int
        y = M[:,i]
        plt.xlabel('ZTime')
#         plt.ylabel('PC%d'%(i+1))
        plt.ylabel('PC%d(%.1f%%)'%(pi+1,vara[pi]*100))
        l = plt.scatter(x,y,c=COL_RGB,)

        for i,(xx,yy,lab) in enumerate(zip(x,y,labs)):
            plt.annotate(lab,xy=(xx,yy), )
    #     plt.legend()
        plt.grid()
    recs = [pyvis.mpl.patches.Rectangle((0,0),1,1,fc=c) for c in COL_LST]
    fig.legend(recs,COL_LAB)    

    return fig    


# In[28]:


#### Cluster Profiling Utilities
fname = 'key.gene'
try:
    geneKey = pd.read_table(fname)
except:
    print "[WARN] Cannot find file:%s"%fname
    geneKey = None

def findMarker(df,concise=1,silent = 0,geneKey=geneKey,how='inner'):
    ''' "geneKey" needs to have at least 'Gene Name' and 'Bio Name' 
    '''
    import IPython.display as ipd
    if not isinstance(df,pd.DataFrame):
        df = pd.DataFrame({'Gene Name': df,'Gene ID':df})
    df = df.reset_index().merge(geneKey,how=how).set_index('index')
    df = df.rename(columns={'Gene ID':'Hit ID',
                       'Gene Name':'Query ID'})
    if concise:
#         df = df[['Hit ID',u'Query ID',u'Bio Name', u'Major role', u'Type of gene']]
        df = df[['Hit ID',u'Query ID',u'Bio Name']]
#         df = df[['Hit ID',u'Query ID',u'Bio Name', u'Major role', u'Type of gene']]
    if not silent:
        print '[MARKER] Found %d/%d' %(sum(~df['Hit ID'].isnull()),len(geneKey))
        ipd.display(df)
    return df

def mapTup(lst,n):
    res = [x[:n] for x in lst]
    return res
def isNovo(lst):
    res = map(lambda x:x=='STRG',mapTup(lst,4))
    return res
def countNovo(df,):
    res = pyutil.collections.Counter([])
    return res

def meta2name(meta,keys=['gtype','light','Age','ZTime']):
    res = pyutil.paste0([meta[k] for k in keys],'_')
    return res

def qc_GeneExpr(exprMat,idx=None,
               gene=None,gRef=None,id_col='Gene Name',
                show_ytick = None,
               condName=None,**kwargs):
    if idx is None:
        assert not(gene is None or gRef is None),'Must specify "gene" and "meta" when "idx" not provided'
        ### Query dataframe with id
        qRes = pyutil.gQuery(gene,gRef,id_col=id_col)
        idx = qRes.index
    show_ytick = show_ytick or len(idx)<=100
    if gene is not None:
        ytick = gene.values 
    elif gRef is not None:
        ytick = gRef.loc[idx][id_col]
    else:
        ytick = idx
        print '[WARN] ytick not defined'
    
    if condName is None:
        xtick = None
    elif isinstance(condName,pd.DataFrame):
        xtick = meta2name(condName)
    else:
        xtick = condName
    ax = pyvis.heatmap(exprMat[idx],
                       xlab='Sample ID',ylab='Gene',
                       xtick = xtick,
                       ytick = ytick if show_ytick else None,**kwargs
                      )
    return ax


# In[147]:


def stdNorm(X):
    X = meanNorm(X).copy()
    STD = np.std(X,axis=1,keepdims=1); pos = np.squeeze(STD>0);
    X[pos] = X[pos]/STD[pos]
    return X
def meanNorm(X):
    X = (X-X.mean(axis=1,keepdims=1))
    return X
def ctNorm(X):
    X = (X-X.mean(axis=1,keepdims=1))
    return X
def identityNorm(X):
    return X
def fit_BGM_AllNorm(C,normLst=None,algoLst=None,ALI='Test',**kwargs):
    if normLst is None:
        normLst = [stdNorm,meanNorm,ctNorm,identityNorm]
    if algoLst is None:
        algoLst = ['DPGMM','DDGMM','GMM',]
    mdls = {}
    for normF in normLst:
        mdls[normF.__name__] = fit_BGM(C,normF=normF,
                                       ALI=ALI,
                                       algoLst = algoLst,
                                       **kwargs)
#     np.save(ALI,mdls,)        
    return mdls


def fit_BGM(C,
            ALI = 'Test',
            normF = identityNorm,
            stdPer = 0,
            rowName=None,
            colName=None,
            nClu = 30,
            maxIt = 250,
           algoLst = None):
    if algoLst is None:
        algoLst = ['DPGMM','DDGMM','GMM',]    
    try:
        DIR,ALI = ALI.rsplit('/',1)
    except:
        DIR='.'
    os.system('mkdir -p %s'%(DIR))
    if isinstance(C,pd.DataFrame):
        rowName,colName,C = C.index.values, C.columns, C.values        
        pass
    if stdPer > 0 :
        assert stdPer < 100,'Percentile must < 100, Got %d instead'%stdPer
        (MEAN,STD,CV),_ = qc_Avg(C)
        pIdx = STD > np.percentile(STD, stdPer)        
        rowName = np.array(rowName)[pIdx]; C = C[pIdx]
        
    import sklearn.mixture as skmix
    common = {'n_components': nClu,
          'verbose':2,
         'max_iter':maxIt,}
    alpha = .1
    mdlLst = {'DPGMM':skmix.BayesianGaussianMixture(weight_concentration_prior_type='dirichlet_process',
                                        weight_concentration_prior=alpha,
                                       **common),
          'GMM':skmix.GaussianMixture(**common),
          'DDGMM':skmix.BayesianGaussianMixture(weight_concentration_prior_type='dirichlet_distribution',
                                        weight_concentration_prior=alpha,
                                       **common),
         }
    mdls = {}
    X = normF(C)
#     raise Exception('test')
    for mdlName,mdl in mdlLst.items():
        if mdlName not in algoLst:
            continue
        NAME = '%s_stdPer=%d_norm=%s_genre=%s_nClu=%d_maxIt=%d'%(
            ALI,
            stdPer,
            normF.__name__,
            mdlName,
            mdl.n_components,
            maxIt
        )
        print '[MSG] Now Fitting Model:%s'%NAME
        d = {'name': NAME,
             'train_data':X,
             'colName':colName,
             'rowName':rowName,
             'param':{
                 'stdPer':stdPer,
                 'normF':normF.__name__,
                 'nClu':mdl.n_components,
                 'genre':mdlName,
             },
           }
        try:
            logFile = open('%s/%s.log'%(DIR,NAME),'w')
            with pyutil.RedirectStdStreams(logFile):
                mdl.fit(X)
                d.update({'suc':1,'model':mdl})
#             logFile.close()
            print "[SUCC] to fit Model:%s"%(NAME,)
        except Exception as e:
            print "[FAIL] to fit Model:%s due to :'%s'"%(NAME,e)
            d.update({'suc':0})
        mdls[NAME] = d    
        np.save('%s/%s'%(DIR,NAME),d)
    return mdls
def make_qc_Model(vX,tX=None,normF = None):

    ##### Datasets: Training
#     vX = None
    ##### vX: Validation dataset
#     vX = vX[clu==1,:][:,msort[msort['light']=='SD'][msort['Age_int']==2][msort].index]
#     vX = util.preprocess(vX,std=1)
#     print normF,type(normF)
    def qc_Model(#model,train_data,
                suc=1,
                 name='Test',pt=None,
                normF_override=normF,
                tX=tX,
        **d):
        if not suc:
            print '[]Skipping failed Model %s'%name
            return 
        print msgGMM(name=name,**d)    
        fig,axs = plt.subplots(3,2,gridspec_kw={"width_ratios": (.1, .9),
                                                   'wspace':0.1,
                                                    'hspace':0.5,
                                                    'top':0.8
                                                   },
                                  figsize=[14,6])
        normF = normF_override or getattr(modCurr, name.split('_')[0])
        if tX is None:
            tX = kwargs['train_data']
#         if tX is None:
#             tX = kwargs['X']
        axc= axs[0]
        dname = 'Datasets-Training_'
        qcGMM(valid_data=normF(tX),pt=pt,axs=axc,name=dname+name,**d)
        axc= axs[1]
        dname = 'Datasets-Validation_'
        qcGMM(valid_data=normF(vX),pt=pt,axs=axc,name=dname+name,**d)
        
        axc= axs[2]
        dname = 'Datasets-ValidMinusTraining'
        vD = normF(tX)-normF(vX)
#         vD = -(normF(tX)-normF(vX))
        qcGMM(valid_data=vD,pt=pt,axs=axc,name=dname+name,**d)
        
    return qc_Model


# In[164]:


def qc_Sort(fname=None,df=None,cname = 'test',vmin=-2,vmax=2,title = None,**heatargs):
    if df is None:
        df = pyutil.readData(fname)
        if title is None:
            title = '[file]%s'%fname
    heatargs.update(
        {'vmin':vmin,
         'vmax':vmax,
         'cname':cname,
        }    )
    C = df.values
    (M,V,CV),_= qcAvg(C,silent=0)
    plt.suptitle(title)
    inter = -len(C)//1000
    fig,axs= plt.subplots(3,1,figsize=[14,6],gridspec_kw={'hspace':0.3})
    axs=axs.flat
    pyvis.heatmap(C[V.argsort()][::inter],transpose=1,
                 main='sorted by Varaince',ax=axs[0],**heatargs)

    pyvis.heatmap(C[CV.argsort()][::inter],transpose=1,
                 main='sorted by CV',ax=axs[1],**heatargs)

    pyvis.heatmap(C[M.argsort()][::inter],transpose=1,
                 main='sorted by Average',ax=axs[2],**heatargs)
    return fig,axs


# In[157]:


class countMatrix(pyutil.util_obj):
    def __init__(self,**kwargs):
        super(countMatrix, self).__init__(**kwargs)
    @classmethod
    def from_npy(cls,fname):
        dd = np.load(fname).tolist()
        return cls(**dd)
    def colName_short(self,condName=None):
        condName = self.colName if condName is None else condName        
        try:
            condName = pyutil.meta2flat([ 
                [x[1] for x in y] 
                    for y in pyutil.flat2meta(condName)],
                seps=['_'])
        except Exception as e:
            print '\n \[WARN\] unable to simplify condName. Exception:%s'%e
        finally:
            return condName
    def heatmap(self,C=None,vlim=[-2,2],**kwargs):
        C = self.train_data if C is None else C
        condName = self.colName_short()
#         im = pyvis.heatmap(C[cidx][sidx],
        im = pyvis.heatmap(C,
#                            ylab=(None if not i else 'Gene'),
#                            ytick = (None if not i else gCur['Gene Name']),
                           xlab='Condition',
                           xtick=condName,
                           transpose=1,
                          vmin=vlim[0],vmax=vlim[1],
#                           ax=ax
                          ) 
        return im
def loadMatrix(fname):
    ali, ext = fname.rsplit('.',1)
    if ext == 'npy':
        dd = np.load(fname).tolist()
    return dd
import countMatrix as ctMat
# from countMatrix import countMatrix
countMatrix = ctMat.countMatrix
sortLabel = ctMat.sortLabel
# import countMatrix;reload(countMatrix)
# from countMatrix import countMatrix

def qc_ModelDict(dd=None,fname=None,ali=None,geneKey=None,DIR=None,
                 clu = None,cluMax = 100,
                 vlim= None
                ):
    if isinstance(geneKey,dict ):
        geneKey = pd.Dataframe.from_dict(geneKey)
        geneKey[1] = geneKey.index; 
        geneKey.rename(columns={0:'Bio Name',
                                1:'Gene Name'})
    if dd is None:
        dd = countMatrix.from_npy(fname)
        ali = fname.rsplit('.',1)[0]
        
    if dd.suc ==0:
        print '[WARN] this model is empty due to a failure %s'%dd['name']
        return
    if vlim is None:
        vlim = np.span(dd.train_data,p=99.9)
#         geneKey.rename({})
    sper = 0
#     ali = NBNAME+'_h%d_'%75+
    if ali is not None:
        ali = ali.rsplit('/',1)[-1]
    else:
        ali = dd.__dict__.get('name','test')
        if isinstance(ali,list):
            ali = ':'.join(ali)
    DIR = os.path.abspath(DIR or '.')
    
#     os.system('mkdir -p %s/src'%DIR)
#     os.system('mkdir -p %s'%DIR)
    print '[ALI]',DIR,'/',ali
    
    mdl,tX = dd.model,dd.train_data; tXsd = stdNorm(tX)
    gRef,condName = dd.rowName,dd.colName_short()
    
    #### Process rowName
    gRef = pd.DataFrame({'Gene Name':gRef,'Gene ID':gRef})
    if geneKey is not None:
        gRef = findMarker(gRef, geneKey=geneKey,silent=1,how='left',concise=1)
        gRef['isMarker']=~gRef['Bio Name'].isnull()
        gRef = gRef.rename(columns={'Query ID':'Gene Name'}).drop('Hit ID',1)
        print '[GREF]',len(gRef)
    
    if isinstance(mdl,list):
        print dd.nCol
        tX = tX[:,:dd.nCol[0]]; nidx = np.isnan(tX[:,0])
        mdl = mdl[0]
        if any(nidx):        
            tX,nX = tX[~nidx],tX[nidx]; nn = sum(nidx)
            Y = mdl.predict(tX); s = mdl.score_samples(tX); 
            Y,pos = sortLabel(Y,tX)
            Y = np.hstack([Y,[max(Y)+1]*nn]); s = np.hstack([s,[-1]*nn]); sbin = s> np.percentile(s,sper); 
            print tX.shape,sbin.shape
        else:
            Y = mdl.predict(tX); s = mdl.score_samples(tX); sbin = s> np.percentile(s,sper); 
            Y,pos = sortLabel(Y,tX)
        tX = dd.train_data
    else:
        Y = mdl.predict(tX); s = mdl.score_samples(tX); sbin = s> np.percentile(s,sper); 
        Y,pos = sortLabel(Y,tX)
    
    
    # pcommon= {}
    try:
        os.system('mkdir -p %s/%s'%(DIR,ali))
#         _ , ali = ali.rsplit('/',1)
        CWD= os.getcwd()
        _ = os.chdir('%s/%s'%(DIR,ali))
        os.system('mkdir -p src/')
        OFILE = open('main.md','w')
        ExcelFile= pd.ExcelWriter('main.xlsx', engine='xlsxwriter')
        with pyutil.RedirectStdStreams(OFILE):
#             print dd.param
            parDF = dd.param if isinstance(dd.param, list) else [dd.param]
            parDF = pd.DataFrame(parDF)
#             print '[pAss]'
            print '\n',pyutil.pd2md(parDF)
#             for k,v in .items():
#                 print '%s:%s\n'%(k,v)
            print 'Directory: %s \n \n  Model Name: %s'%(DIR,ali)
            print '\n [.xlsx](main.xlsx)',
            print '[.tar.gz](main.tar.gz)',
            for clu in range(-1,max(Y)+1):
                if clu==cluMax:
                    break 
                fig,axs = plt.subplots(2,1,figsize=[max(7,min(14,len(tX)/3.)),
                                                   max(5,min(18,len(tX.T)/1.5))],
                                       sharex='all',
                                     gridspec_kw={'bottom':0.28,'top':0.8,
                                                 'left':0.2}
                                     )        
                if clu == -1:
#                     Y,pos = sortLabel(Y,tX)
                    cidx = Y>-1                    
                    gCur = gRef[cidx]
                    sidx = np.argsort(Y);gCur = gCur.iloc[sidx]
                    cluName = 'background.gene'
                else:
                    cidx = Y==clu 
                    cidx = cidx & sbin
                    gCur = gRef[cidx]
                    sidx = np.argsort(s[cidx])[::-1];gCur = gCur.iloc[sidx]
                    cluName = 'clu%03d.gene'%(clu)
                print '\n Cluster:%d'%(clu),'\n','Gene Count:%d'%len(gCur)
                if len(gCur)==0:
                    continue
                    
#                 gCur = pd.DataFrame({'Gene Name':gCur})
#                 if clu>=0:
                if geneKey is not None:
                    gMark = gCur[gCur['isMarker']==1]
#                     gMark = findMarker(gCur['Gene Name'],geneKey=geneKey,silent=1,how='right')
                    gMark.to_excel(ExcelFile,cluName,index=True,startcol=2)       
                    print '\n',pyutil.pd2md(gMark)
#                     gCur = gMark
#                     gCur = gMark.rename(columns={'Query ID':'Gene Name'})                
                gCur[['Gene Name']].to_csv('src/%s'%cluName,index=0,)
                gCur[['Gene Name']].to_excel(ExcelFile,cluName,index=True,)                
                
                figList = []
                matDict ={'raw':tX,'stdNorm':tXsd}
                sheet_curr = ExcelFile.sheets[cluName]
                for i,k in enumerate(['raw','stdNorm']):
                    C = matDict[k]; ax =axs[i]                    
                    im = pyvis.heatmap(C[cidx][sidx],
                                       ylab=(None if not i else 'Gene'),
                                       ytick = (None if not i else gCur['Gene Name']),
                                       xlab='Condition',xtick=condName,transpose=1,
                                      vmin=vlim[0],vmax=vlim[1],
                                      ax=ax
                                      )  
                    dd.addBox(ax=ax)
#                     figList +=[FFname]
                    
                    plt.colorbar(im)
                    plt.title(k)
                plt.suptitle('Cluster %d'%clu,y=1)
#                 try:
#                 fig.tight_layout()
#                 except:
#                     print '\n \[WARN\] tight_layout() failed, legend may not display properly'
#                     pass                
                FFname = 'src/clu%03d.png'%(clu,)
                FigMd = pyutil.showsavefig(fname=FFname)
                print '\n',(FigMd) ## remove directory name                
                sheet_curr.insert_image(0, 7, FFname)                
                plt.show()
                plt.close()
        
        ExcelFile.save()
        ExcelFile.close()
        os.system('pdext {fname} html'.format(fname=OFILE.name) )
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        raise e
    finally:
        os.chdir(CWD)
#     os.chdir(sys.path[0])
#     os.system('pdext {ali}.md pdf'.format(ali=ali) )
def log2p1(x):
    res = np.log2(x+1)
    return res


# In[142]:


#### DataSets Management
def dfContrast(dfRef,dfObs):
    ''' Contrast two DataFrames
    '''
    C = dfObs.values - dfRef.values
    df = pd.DataFrame(C); df.set_index(dfObs.index,inplace=1)
    df.columns = pyutil.metaContrast(dfRef.columns,dfObs.columns)
    return df


# In[166]:


if __name__=='__main__':
    get_ipython().system(u'jupyter nbconvert --to python tmp.ipynb')
# !python compile_meta.ipynb && echo '[succ]'


# In[154]:


sorted(['001','005','003','004'])

