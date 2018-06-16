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
def fillNA_withRef(dd, ref):
    naX,naY = np.where(np.isnan(dd.train_data))
    if len(naX)==0:
        pass
    else:
        naX_mapped,naX_uniq = pyutil.reindexZero(naX)
        naY_mapped,naY_uniq = pyutil.reindexZero(naY)
        odf = ref.reindex(index=dd.rowName[naX_uniq],columns=dd.colName[naY_uniq])
        dd.train_data[(naX,naY)] = odf.values[(naX_mapped, naY_mapped)]
    return dd    
class countMatrix(pyutil.util_obj):
    name = None
    def __init__(self,arg1=None,**kwargs):
        super(countMatrix, self).__init__(**kwargs)
#         print self.train_data[:5]
#         df = pd.DataFrame(self.train_data,index=None).set_index(self.rowName)#.drop(['index'],1)
#         self.setDF(df)
    def setDF(self,df):
        self.colName = df.columns #= self.colName
        self.rowName = df.index
        self.df = df 
        self.train_data = self.df.values
        if not isinstance(self.name, list):
            self.nCol =  len(self.colName)
        self.nRow = len(self.train_data)
        return self
    @classmethod
    def from_DFrame(cls,df=None,fname=None):
        if df is None:
            assert fname is not None,'[ERR] must specify one of "df" or "fname" '
            df = pyutil.readData(fname)
        ins = cls()
        ins.setDF(df)
        return ins
    @classmethod
    def from_npy(cls,fname):
        dd = np.load(fname).tolist()
        ins = cls.from_dict(dd)
        return ins
    @classmethod
    def from_dict(cls,dd):
        C,colName,rowName = [ dd.pop(key) for key in ['train_data','colName','rowName']]
        df = pd.DataFrame(C,index=None).set_index(rowName);df.columns=colName
        #.drop(['index'],1)
        return cls(**dd).setDF(df)
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
    def heatmap(self,C=None,vlim=[-2,2],cname = 'test',
                reorder=0,
                **kwargs):
        reorder and self.reorder();
        C = self.train_data if C is None else C
        condName = self.colName_short()
#         im = pyvis.heatmap(C[cidx][sidx],
        im = pyvis.heatmap(C,
#                            ylab=(None if not i else 'Gene'),
#                            ytick = (None if not i else gCur['Gene Name']),
                           xlab='Condition',
                           xtick=condName,
                           transpose=1,
                           cname = cname,
                          vmin=vlim[0],vmax=vlim[1],
#                           ax=ax
                          ) 
        
        addBox(dd= self,ax=im.axes)
        return im
    def sortedLabel(self):
        C = self.train_data
        clu = self.model.predict(C)
        clu, pos = sortLabel(clu,C)
        return clu
    def reorder(self):
        clu = self.sortedLabel()       
        self.setDF(self.df.iloc[clu.argsort()])
        return self
    def condCount(self):
        if isinstance(self.nCol, pyutil.collections.Iterable):
            n = sum(1 for x in self.nCol)
        else:
            n = 1
        return n
countMatrix.fillNA_withRef = fillNA_withRef
countMatrix.addBox = addBox

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