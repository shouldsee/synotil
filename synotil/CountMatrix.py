import numpy as np
import pandas as pd
# import util; reload(util); from util import *
import pymisca.util as pyutil
import pymisca.vis_util as pyvis
# from util import sortLabel

##### Utility functions
def sortLabel(Y,X,#return_pos=0
             ):
    ### Sorting by maximum
#         X = vX
    X  = X-X.mean(axis=1,keepdims=1)
    coord = np.arange(X.shape[-1])[None]
    wt_X = (X == X.max(axis=1,keepdims=1))*coord
#         wt_X = X * coord
    cis = list(range(max(Y)+1))
    pos =  [wt_X[Y == ci,:].mean() for ci in cis]
    sidx = np.argsort(pos)
#     if return_pos:
#         return pos[sidx]
#     else:
    pj = dict(zip(sidx,cis))
    Y = np.vectorize(pj.get)(Y)
    return Y,np.take(pos,sidx)

def addBox(dd=None,nCol=None,nRow=None,ax=None):
#     dd = self
    if dd is not None:
        nCol, nRow = dd.nCol, dd.nRow
        if dd.condCount() <= 1:
            return ax    
    if not isinstance( ax, pyutil.mpl.axes.Axes):
        ax = ax.axes        
    x = np.cumsum([0,] + list(nCol))
    for a, b in zip(x[:-1],x[1:]):
#         print a,b,nRow  ##Debug
        r = pyutil.mpl.patches.Rectangle(( -0.5, a - 0.5),
                                         nRow, b-a,
                                         fill=False,color='red')
        ax.add_artist(r)
    return ax

def predictWithCut(mdl,data,CUTOFF=0.5):
    '''Cutoff a BGM by weights of components GMMs
#### This should be added to the BGM class
'''

    W = mdl.weights_.copy()
    idx = W.argsort()[::-1]
    mdl.active = idx[:np.where((W[idx].cumsum()>CUTOFF))[0][0]+1]
    P = mdl.predict_proba(data)
    P = P[:,mdl.active]
    clu = P.argmax(axis=1)
    return clu

def fillNA_withRef(dd, ref):
    naX,naY = np.where(np.isnan(dd.train_data))
    if len(naX)==0:
        pass
    else:
        naX_mapped,naX_uniq = pyutil.reindexZero(naX)
        naY_mapped,naY_uniq = pyutil.reindexZero(naY)
        odf = ref.reindex(index=dd.index[naX_uniq],columns=dd.columns[naY_uniq])
        dd.train_data[(naX,naY)] = odf.values[(naX_mapped, naY_mapped)]
    return dd    


#### !!!This is not working
import decorator


import wrapt
@wrapt.decorator
def wrapDFMethod( f, instance, args, kwargs):
# @decorator.decorator
# def wrapDFMethod(f,*args,**kwargs):
    '''Wrap a pd.DataFrame method to return countMatrix instance
'''
#     self = args[0]
    self = instance
    res = f(*args,**kwargs)
    obj = self.setDF(df=res)
    return obj
#     def newMethod(self,*args,**kwargs):
#         df = DFMethod(self,*args,**kwargs)
#         obj= self.setDF(df=df)
#         return obj
#     newMethod.__name__ = DFMethod.__name__
#     newMethod.__doc__ = DFMethod.__doc__    
#     return newMethod

# def wrapDFMethod(DFMethod):
#     '''Wrap a pd.DataFrame method to return countMatrix instance
# '''
#     def newMethod(self,*args,**kwargs):
#         df = DFMethod(self,*args,**kwargs)
#         obj= self.setDF(df=df)
#         return obj
#     newMethod.__name__ = DFMethod.__name__
#     newMethod.__doc__ = DFMethod.__doc__    
# #     newMethod.func_code = DFMethod.func_code
#     ####  ====TBC====
# #     newMethod.im_func.func_code = DFMethod.im_func.func_code
#     return newMethod

# class countMatrix(pyutil.pd.DataFrame):

import pymisca.bio as pybio
import StringIO

class countMatrix(pd.DataFrame):
    ''' A decorated pandas.DataFrame class that
knows how to plot itself
'''
#     name = None
#     look = None
#     @property
#     def _constructor(self):
#         return countMatrix
    def __init__(self,C=None, rowName=None, colName = None,
                 name = None,
                 look = None,
                 cmap = None,
                 vlim = None,
                fname = None,
                 model = None,
                 colMeta = None,
                 rowMeta = None,
                 height = 1.,
                 **kwargs):
        df = pyutil.init_DF( C=C,rowName=rowName,colName = colName)  
#         print self.__class__,countMatrix, isinstance(self,countMatrix,)
        super(countMatrix, self).__init__(df)
        self.name_ = name
        self.look = look
        self.cmap = cmap
        self.vlim = vlim
        self.fname = fname
        self.model = model
        self.height = height 
#         self.colMeta_ = colMeta
        self.set__colMeta(colMeta)
        self.set__rowMeta(rowMeta)
#         self.rowMeta_ = rowMeta
        self.param = {'normF':'identityNorm',
                     }        
        self.set_config( test=None,**kwargs)
        self.test= None

        assert self.name !='test','Track name cannot be %s'% self.name


    def set__colMeta(self,colMeta,castF = unicode):
        if colMeta is not None:
            colMeta = colMeta.reindex(self.columns).astype(castF)
        self.colMeta_ = colMeta
        return 

    def set__rowMeta(self,rowMeta,castF = unicode):
        if rowMeta is not None:
            rowMeta = rowMeta.reindex(self.index).astype(castF)
        self.rowMeta_ = rowMeta
        return  
    
    def set__vlim(self,vlim=None):
        if vlim is None:
    #         vlim=[-3,3]
            vlim = pyutil.span(self)
            if vlim[0] < 0:
                vlim = - np.mean(vlim) + vlim
        self.vlim = vlim
        return     
    @property
    def isIndexed(self):
        res = self.look not in ['gtf','feat','feats']
        return res
    @property
    def colMeta(self):
        return self.colMeta_
#     .reindex(self.columns)
    @property
    def rowMeta(self):
        return self.rowMeta_
#     .reindex(self.index)
    
    def relabel(self,colLabel=None,rowLabel=None):
        '''Decorate the data frame to be ready for plotting
    '''
        if colLabel is None and rowLabel is None:
            if self.colMeta is not None:
                self.columns = self.colMeta_.index
            if self.rowMeta is not None:
                self.index   = self.rowMeta_.index
#             raise Exception('must supply one argument')
        if colLabel is not None:
            self.columns = self.colMeta[colLabel]             
        if rowLabel is not None:
            self.index   = self.rowMeta[rowLabel] 
        return self        
    
    @property
    def name(self):
        return str(self.name_)
#     .split('/')[-1]
    
    def name4param(self,keys=None):
        d = self.param
        if keys is not None:
            d = pyutil.dictFilter(d,keys) 
        name = pyutil.dict2flat(d)
        self.set_name(name)
        return self.name
    def formatName(self,name=None,**kwargs):
        name = self.name if name is None else name
        res = pyutil.formatName(name,**kwargs)
        return res
    
    def __getstate__(self):
        d = dict(self.__dict__)
#         del d['logger']
        return d

    def __setstate__(self, d):
        self.__dict__.update(d) # I *think* this is a safe way to do it
        
    def set_name(self,name):
        self.name_ = name
        return name
    def set__name(self,*a,**kw):
        return self.set_name(*a,**kw)
    
    def set_config(self, **kwargs):
        if 'name' in kwargs:
            self.set_name(kwargs.pop('name'))
        for k,v in kwargs.items():
#             if v is not None:
            setattr(self,k,v)        
    def get_config(self):
        config = {'name':self.name,
                  'fname':self.fname,
                 'look':self.look,
                 'cmap':self.cmap,
                 'vlim':self.vlim,
                 'model':self.model,}
        return config
    def __repr__(self, ):
        s = '<countMatrix: name=%s>'% self.name
        return s
    
    def setDF(self,df,update=1):
#         config = self.get_config()
        configDict= self.__getstate__()
        [configDict.pop(k) for k in list(configDict) if k.startswith('_')]
        
        if isinstance(df,pd.Series):
            df = df.values
        if isinstance(df, pyutil.np.ndarray):
            
            if update and (df.shape == self.shape):
                self.loc[:,:] = df
                res = self
                
            else:
                res = pyutil.init_DF(df, rowName = self.index)
        else:
            res = df
        if not isinstance(res, self.__class__ ):
            assert isinstance( res, pd.DataFrame)
#             print self.__class__, res.__class__
            res = self.__class__.from_DataFrame(df = res,)            
        res.__setstate__(configDict)
        res.set__colMeta(res.colMeta)
        res.set__rowMeta(res.rowMeta)        
#         res.set_config(**config)
        return res
    def subSampleColumns(self,stepSize=5):
        df = self.reindex(columns = self.columns[::stepSize])
        return df
    

    @classmethod
    def from_DataFrame(cls,df=None,fname=None,name=None,index_col = None, **kwargs):
        if df is None:
            assert fname is not None,'[ERR] must specify one of "df" or "fname" '
            df = pyutil.readData(fname,**kwargs)
            name = pyutil.os.path.basename(fname).rsplit('.',1)[0]
        elif isinstance(df,pd.Series):
            df = df.to_frame()
        if index_col is not None:
            assert index_col in df
            df.set_index(index_col,drop=0,inplace=True
                        
                        )
#         print 3,type(df)
        ins = cls(C=df.values,
                 colName=df.columns,
                 rowName=df.index,
                  name=name,
                  fname=fname,
                 )
        return ins
    
    @classmethod
    def from_npy(cls,fname,**kwargs):
        dd = np.load(fname).tolist()
        ins = cls.from_dict(dd,**kwargs)
        ins.fname = fname
        return ins
    @classmethod
    def from_dict(cls,dd,**kwargs):
        C,colName,rowName = [ dd.pop(key) for key in ['train_data','colName','rowName']]
        dd.update(kwargs)
        df = cls(C=C,colName=colName,rowName=rowName,**dd)
        return df
    @classmethod
    def from_deepcache(cls, d= None, fname = None):
        if fname is not None:
            d = np.load(fname)
            d = dict(d)
        else:
            assert d is not None,'must specify one of variables: "d" or "fname"'
        deepdict = d
        C = deepdict.get('matrix',np.array([[]]))
        colName= deepdict.get('labels',None)
        rowName =deepdict.get('which is row name?',None)
        df  = pyutil.init_DF(C=C,colName=colName,rowName=rowName)
        ins = cls.from_DataFrame(df=df)
        return ins
    def colName_short(self,condName=None):
        condName = self.columns if condName is None else condName        
        try:
            condName = pyutil.meta2flat([ 
                [x[1] for x in y] 
                    for y in pyutil.flat2meta(condName)],
                seps=['_'])
        except Exception as e:
            print '\n [WARN] unable to simplify condName. Exception:%s'%e
        finally:
            return condName
        
    def dropAll(self):
        self.drop(self.columns, 1, inplace=True)
        return self
    def toGeneList(self):
        self.dropAll()
        self[self.fname] = 1
        return self
    def heatmap(self,C=None,vlim=None,
                cname = 'test',
                xlab='Condition',
                ylab = '',
                ytick = None,xtick = None,
                reorder=0,ax=None,transpose=0, ### [buggy] pymisca updated transpose
                tickMax=100,
                short = 1,
                **kwargs):
        vlim = self.vlim if vlim is None else vlim
        reorder and self.reorder();
        C = self.values if C is None else C
        C = C.T ### [buggy] pymisca updated transpose

        if xtick is None:
            xtick = self.colName_short() if short else self.columns
            
        if ytick is None:
            ytick = self.index if len(self) < 500 else None
            
        if transpose:
            xtick,ytick = ytick,xtick
            xlab,ylab = ylab,xlab
        
#         im = pyvis.heatmap(C[cidx][sidx],
        im = pyvis.heatmap(C,
#                            ylab=(None if not i else 'Gene'),
#                            ytick = (None if not i else gCur['Gene Name']),
                           xlab=  xlab,
                           ylab = ylab,

                           transpose=transpose,
                           cname = cname,
                           vlim = vlim,
                           xtick=xtick,
                           ytick =ytick,
                           tickMax = tickMax,
                           ax=ax,**kwargs
                          ) 
        
        if 0:
            addBox(dd= self,ax=im.axes)
        else:
#             print '[DEV] addBox() is disabled fodr now'
            pass
        return im
    ###### [TBC] pd_extra
    def to_bed(self,path_or_buf=None,index=0,header=None,sep='\t',
               escapechar='\\',doublequote=False,
               **kwargs):
        return self.to_csv(path_or_buf=path_or_buf,
                          index=index,
                           escapechar=escapechar,
                           doublequote=doublequote,
                          header=header,sep=sep,**kwargs).replace(escapechar,'')
    def rename_index(self,ind):
        return self.rename(index={x: ind for x in self.index})    
    
    def to_feats(self,**kwargs):
        seqrecs = pybio.read_gtf(StringIO.StringIO(self.to_bed(**kwargs)))
        if len(seqrecs)==0:
            return seqrecs
        else:
            assert len(seqrecs)==1,'more than 1 chromosomes (seqRecord) found'
            feats = seqrecs[0].features
#             feats = [pybio.gene2transcript(x) for x in seqrecs[0].features]
            return feats

    def copy(self,deep = True):
# #         df = super(countMatrix, self).copy(deep=deep)
        df = pd.DataFrame.copy(self,deep=deep)
        dd = self.setDF(df=df)        
        return dd
    
    def qc_Avg(self):
        '''Add summary dataFrame as an attribute
    '''
        import util as sutil
        (M,SD,CV), _ = sutil.qc_Avg( self,silent=1)
        M = np.nanmean(self.values,axis=1)
        SD = np.nanstd(self.values,axis=1)
        CV = SD/M
#         SD = np.nanmedian(self.values,axis=1)
        MAX=np.nanmax(self.values, axis=1)
        MIN=np.nanmin(self.values, axis=1)
        effSize = np.sum( ~np.isnan(self.values), axis = 1)
        SE = SD/effSize ** 0.5
        # df = pd.concat([M,SD,CV],axis=1)
        MSQ = M**2 + SD**2
        df = pd.DataFrame({'M': M,
                           'SD':SD,
                           'MSQ':MSQ,
                           'RMSQ':np.sqrt(MSQ),
                           'MAX':MAX,
                           'MIN':MIN,
                           'CV':CV,
                           'SE':SE,
#                            'MSQ':MSQ,
                          'effSize':effSize,
                          })
#         df = pd.DataFrame({'M':M.values,
#                            'SD':SD.values, 
#                            'CV':CV.values})
        df = df.set_index( self.index)
        df['per_SD'] = pyutil.dist2ppf(df.SD)
        df['per_M'] = pyutil.dist2ppf(df.M)
        df['per_MSQ'] = pyutil.dist2ppf(df.MSQ)
        df['per_CV'] = pyutil.dist2ppf(df.CV)
        df['per_MAX'] = pyutil.dist2ppf(df.MAX)
        df['per_MIN'] = pyutil.dist2ppf(df.MIN)
        df['per'] = df['per_SD']

        self.summary = df
        return self        

#     def fillna(self,fill,inplace=1,**kwargs):
#         pd.DataFrame.fillna(self,fill,inplace=inplace,**kwargs)
#         return self
#     def apply(self,*args,**kwargs):
#         df = pd.DataFrame.apply(self,*args,**kwargs)
    #### ====<
    ### ====> Wrappers for pandas.DataFrame native methods
    @wrapDFMethod
    def sortedLabel(self,clu = None,CUTOFF=None):
        C = self.values
        if clu is None:
            if CUTOFF is not None:
                clu = predictWithCut(self.model,C, CUTOFF=CUTOFF)
            else:
                clu = self.model.predict(C)
                
        clu, pos = sortLabel(clu,C)
        return pd.DataFrame(clu).set_index(self.index)

    @wrapDFMethod
    def predict(self, C= None,CUTOFF=None):
        C = self.values if C is None else C
        if CUTOFF is not None:
            clu = predictWithCut(self.model,C, CUTOFF=CUTOFF)
        else:
            clu = self.model.predict(C)
        return pd.DataFrame(clu).set_index(self.index)
    def predict_proba(self,C=None,asDF=1):        
        C = self.values if C is None else C
        if hasattr(self.model,'predict_proba'):        
            clu = self.model.predict_proba(C)
        elif hasattr(self.model,'predict'):
            clu = self.model.predict(C)
            clu = pyutil.oneHot(clu)
        if asDF:
            clu = pd.DataFrame(clu).set_index(self.index)    
        return clu
    def reorder(self):
        clu = self.sortedLabel()
        clu = clu.sort_values(0,)
        res = self.setDF( self.reindex(clu.index, ))
        return res
#     def condCount(self):
#         if isinstance(self.nCol, pyutil.collections.Iterable):
#             n = sum(1 for x in self.nCol)
#         else:
#             n = 1
#         return n    
    def sanitiseIndex(self,ptn):
        '''Apply a Regex capturer to sanitise/standardise the index
'''
        res = [x[0] for x in self.index.map(ptn.findall)]
        res = self.set_index([res],inplace=True)
        return self
    def flatSubset(self, keep,negate = 0, which = 'columns'):
        vals = getattr(self,which)
        isStr = isinstance(vals,str)
        if isStr:
            vals = [vals]
        vals = map(str,vals)
        vals = pyutil.flatSubset( vals, 
                                         keep = keep, 
                                         negate =negate)
        if isStr:
            vals =vals[0]
        setattr(self, which, vals)
        return self
    @wrapDFMethod
    def transpose(self,):
        return super(countMatrix,self).transpose()

countMatrix.fillNA_withRef = fillNA_withRef
countMatrix.addBox = addBox
countMatrix.mergeByIndex = wrapDFMethod(pyutil.mergeByIndex)
countMatrix.filterMatch = wrapDFMethod(pyutil.filterMatch)
countMatrix.to_tsv = pyutil.to_tsv
countMatrix.melt =  pyutil.melt
countMatrix.makeContrastWithMeta =  pyutil.df__makeContrastWithMeta



countMatrix.fillna= wrapDFMethod(pd.DataFrame.fillna)
countMatrix.apply = wrapDFMethod(pd.DataFrame.apply)
countMatrix.applymap = wrapDFMethod(pd.DataFrame.applymap)
countMatrix.merge = wrapDFMethod(pd.DataFrame.merge)
countMatrix.query = wrapDFMethod(pd.DataFrame.query)
countMatrix.sort_values = wrapDFMethod(pd.DataFrame.sort_values)
countMatrix.reindex = wrapDFMethod(pd.DataFrame.reindex)
countMatrix.astype = wrapDFMethod(pd.DataFrame.astype)

lst = ['clip','mean','max','min','std','dropna','get','rename',
#        'transpose'
      ]
for name in lst:
    setattr(countMatrix,name,
            wrapDFMethod( getattr( pd.DataFrame, name)) )

    
def queryBed(bedDF, qid):
    '''query a bed file with a list of accessions 
'''
    if isinstance(bedDF,str):
        bedFile = bedDF 
        bedDF = countMatrix.from_DataFrame(
            fname = bedFile,ext='tsv',header=None).reset_index()
    resDF = bedDF.set_index(3,drop=0).reindex(qid,).reset_index(drop=1)
    return resDF

# countMatrix.mergeByIndex = wrapDFMethod(pyutil.functools.partial(pd.DataFrame.merge,
#                                                    left_index=True,
#                                                    right_index=True,)
#                                        )

countMatrix.append = wrapDFMethod(pd.DataFrame.append)

def vstack(lst,how='outer', as_index=0, as_df =0, **kwargs):
#     if how = 'left':
    if as_index:
#         lst = [x.setDF( pd.DataFrame(x.index,index=x.index)) for x in lst]
#         lst = [x.setDF( x.iloc[:,:1]) for x in lst]
        lst = [x.drop( x.columns,1) for x in lst]
    res = reduce(lambda x,y :x.mergeByIndex(y,how=how,**kwargs),lst)
    
    if as_index:
        res = res.index    
    return res

def concat(*args,**kwargs):
    df = pd.concat(*args,**kwargs)
    obj= countMatrix.from_DataFrame(df=df)
    return obj


class geneList(countMatrix):
    def __init__(self,C=None,colName=None,rowName=None,**kwargs):
        super(geneList,self).__init__(C=len(rowName) * [1],rowName=rowName, **kwargs)
#         self.astype('int')
        col = self.columns[0]
        self[col]= self[col].astype(int)
        add_predictProba(self)
#         self.add_predictProba()
# geneList.        
def add_predictProba(glist):
    mdl = pyutil.util_obj()
    def f(vals):
        res = np.nan_to_num(vals).astype(int)
        res = pyutil.oneHot(res)
        return res
    mdl.predict_proba =  f
    glist.model = mdl
    return glist

        
#### Obsolete

def makeCountMatrix(arg):
    if isinstance(arg,str):
        res = countMatrix.from_npy(arg)
    elif isinstance(arg,dict):
        res = countMatrix(**arg)
    else:
        res = arg
    return res

def mergeCountMatrix(dd1,dd2,vlim=[-2,2],how='outer'):
    ''' Merge two countMatrix
    '''
    dd1 = makeCountMatrix(dd1); dd2 = makeCountMatrix(dd2) 
#     dd1 = countMatrix(dd1); dd2 = countMatrix(dd2) 
    C1,C2 = dd1.train_data, dd2.train_data
    colName = np.hstack([dd1.colName,dd2.colName])
    condName = dd1.colName_short(colName)
#     print all(sorted(dd1.rowName)==dd1.rowName)
    g1 = set(dd1.rowName); g2 = set(dd2.rowName); 
    g1Only = g1-g2; g2Only = g2-g1;  gAny = g1 | g2;    gAll = (g1 & g2)
    gAll = list(gAll)
    print 'Union:%d, Shared:%d'%(len(gAny),len(gAll))
    
    ### Find union genes
    ### Assume rowNames are sorted
    idx1 = np.in1d( dd1.rowName, gAll)
    idx2 = np.in1d( dd2.rowName, gAll)        
    r1 = dd1.rowName[idx1]; r2 = dd2.rowName[idx2]
    clu1 = dd1.sortedLabel()
    clu2 = dd2.sortedLabel()
    gAll = dd1.rowName[idx1]   
    CAll = np.hstack([C1[idx1],
                      C2[idx2],
                  ])[(clu1)[idx1].argsort()]

    ### Find C1-only genes
    idx1 = np.in1d( dd1.rowName, list(g1Only))
    gLef = dd1.rowName[idx1]
    CLef = np.hstack([ 
        C1[idx1],
        np.zeros(( sum(idx1),len(dd2.colName)) )*np.nan,
    ])[(clu1)[idx1].argsort()]    
    
    idx2 = np.in1d( dd2.rowName, list(g2Only))
    gRht = dd2.rowName[idx2]
    CRht = np.hstack([ 
        np.zeros(( sum(idx2),len(dd1.colName)) )*np.nan,
        C2[idx2],
    ])[(clu2)[idx2].argsort()]
    
    C   = np.vstack([CLef,CAll,CRht])
    gRef= np.hstack([gLef,gAll,gRht])
    
    dd = {'train_data':C,
                    'colName':colName,
                    'rowName':gRef}

    for k in ['suc', 'name', 'param', 'model','nCol']:
        v1 = getattr(dd1,k);
        v2 = getattr(dd2,k)
        if not isinstance(v1,list):
            v1 = [v1,]
        if not isinstance(v2, list):
            v2 = [v2,]
        dd[k] = v1 + v2
    
    dd = countMatrix(**dd)
    if how =='outer': 
        pass 
    elif how=='left':
        idx = np.hstack([gLef,gAll])
        dd.setDF(dd.df.loc[idx])
#         dd.df = dd.df.reindex(index = ,copy=0)
    elif how=='right':
        idx = np.hstack([gAll,gRht])
        dd.setDF(dd.df.loc[idx])
    elif how=='inner':
        idx = np.hstack([gAll])
        dd.setDF(dd.df.loc[idx])
    
    return dd
