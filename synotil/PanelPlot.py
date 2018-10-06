# reload(scount)

import synotil.CountMatrix as scount
from synotil.CountMatrix import *
import pymisca.vis_util as pyvis; reload(pyvis)
import pymisca.util as pyutil; reload(pyutil)
plt = pyvis.plt
import textwrap
import xlsxwriter
# reload(pyutil)
# pyutil.cluMap = pyutil.mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])


def guessLook(obj):
    '''Guess how DataFrame should be plotted from its shape
'''
    look = 'patch'
    if obj.shape[-1] > 1:
        look = 'matrix'
    else:                
        dtype = obj.dtypes.iloc[0] 
        if dtype == 'int':
            look = 'patch'
        elif dtype == 'bool':
            look = 'tick'
        elif dtype == 'object':
            look = 'text'
        else:
            look = 'line'
    obj.look = look
    return look
import copy
def fixCluster(ele):
    ''' Fix a integer-type track into a color track
'''
    val = ele.fillna(0).values
    N = val.max() - val.min() + 1 
    if N==2:
        shuffle = 0
    else:
        shuffle = 1
    cmap = pyvis.discrete_cmap(N,'Spectral',shuffle = shuffle,seed=0)
    cmap.set_bad('black',1.)
    C = ele.values; 
#                     print C.shape
    C = np.ma.array ( C, mask=np.isnan(C))
#                     print C.shape
    cmat = cmap( C ).squeeze()
#     print cmat.shape
#     ele  =pd.DataFrame(cmat).set_index( ele.index)
#     ele.look = 'patch'
    elenew = ele.setDF(pd.DataFrame(cmat).set_index( ele.index) )
    elenew['clu'] = ele.values; ele =elenew
    ele.cmap = cmap
    ele.look = 'patch'
    return ele

# from matplotlib.ticker import NullFormatter

hide_axis = pyvis.hide_axis
hide_frame= pyvis.hide_frame

# ax.xaxis.set_major_formatter(NullFormatter())

def getPatchTicks(df,tol = 10):            
    gp = df.reset_index().groupby(df.columns[-1])
    m = gp.apply(lambda x: 

             pd.Series(
                 (x.index.values.mean(),len(x.index)) + pyutil.span(x.index.values) 
             ,index=['M','LEN','MIN','MAX'],
             ) 
            )
    m['tickthis'] = m.MAX - m.MIN + 1 < m.LEN + 10
    return m.loc[m.tickthis]

class panelPlot(list):
    '''Class to systematically render plots
'''
    @property
    def shape(self):
        return [x.shape for x in self]
    
    def __init__(self,lst, validate= 1, 
                 _compiled = False, 
                 orig = None,
                 index =None,
                 geneName = None,
                 names  =None,
                 figsize = None,
                 *args,**kwargs):
        first = list(lst[0])
        if isinstance(first[0],str):
            if len(first) == 2:
                names, lst = zip(*lst)
        lst = [ scount.countMatrix(x)
               if not isinstance(x,scount.countMatrix)
               else x
               for x in lst
              ]
        super(panelPlot, self).__init__(lst)
        self.set_names(names)
#         for x in self:
#             x.cmap = None
        self.index = None
        self.validate() if validate else None
        self._compiled = _compiled
        self.orig = orig
        self.index = index
        self.geneName = geneName 
        self.DIR = '.'
        self.fig = {}
        self.figsize = figsize
    def __getslice__(self, *args):
        '''To be extended to allow for easier subIndex
'''
        dct = self.__dict__
        self.__init__(
            super(panelPlot, self).__getslice__(*args),
            **dct)
        return self
    def setAttr(self,**kwargs):
        for k,v in kwargs.items():
            setattr(self,k,v)
        return self
        
    def validate(self,):
        for x in self:
            pass
#             assert hasattr(x,'render'),'element "%s" does not support render()' % x.__repr__()
        return True
    
    def _render(self, obj, axs= None, look = 'patch',
                figsize=[14,6],
                tickMax = 100,
                shortName=0, **kwargs):
        '''[CORE]: Contains rendering methods for different types of data
'''
        if axs is None:
            fig, axs = plt.subplots(1,1,figsize=figsize,sharex='all')
#             ax = axs
            ax =axs
        else:
            ax = axs[1]
        assert isinstance(obj, scount.countMatrix)
        axa,axb,axc = axs
        axa.get_shared_y_axes().join(axa, axb)

        hide_axis(axa); 
        hide_frame(axc);hide_axis(axc); 
        
        vlim = obj.vlim
#         xticks = np.linspace(0,len(obj),20+1)
        xlines = np.linspace(0, len(obj), 20+1)
        axb.set_xticks(xlines, minor = True)
        xticks = None
        xticklabs = None
        if look =='matrix':
            res = obj.heatmap(ax=axb,cname = None)   
            plt.colorbar(mappable=res, cax = axc )
            axcy = axc.get_yaxis(); axcy.tick_right(); axcy.set_visible(True)
    
        elif look =='patch':
            df = obj.fillna(0); 
            res = pyvis.heatmap(
                df.values[None,:,:3],
                ax=axb,
                cname=None,
#                 vlim = [None,None],
#                 vlim = [0,1]
            )           

            dftick = getPatchTicks(df,tol =10)
            if not dftick.empty:
                xticks = dftick.M.values
                xticklabs = dftick.index.values
    
        elif look =='line':
            res = axb.plot(obj.values, 'x')
            if vlim is not None:
                axb.set_ylim(vlim)
                
        elif look == 'tick':
            res = axb.vlines(np.nonzero(obj.fillna(False).values),
                              -1, 1,
                             linewidth = 0.5,
                             )
        elif look == 'text':
            objc = obj.reset_index(drop=True).dropna()
#             axb_ymid =  sum(axb.get_ylim())/2.
            axb_ymid = 0.
            axb.set_ylim(-1,1)
            def plotter(row):
                i,val = row.name, row[0]
                res = axb.text( i + 0. , axb_ymid, 
                         val,
                        horizontalalignment='center',
                        verticalalignment='center',
                        rotation='vertical',
                        size='large',
# #                          fontsize=1,
#                         transform=axb.transAxes
                        ) 
            objc.apply(plotter,axis=1)
#         print axb.get_xlim(),axb.get_ylim()
            
        axb.set_xlim(0-0.5,len(obj)-0.5)
        if look not in ['tick']:
            axb.grid(color='black',axis='x',linestyle='--',which='minor')
        if xticks is not None:
            axb.set_xticks(xticks)
            axb.set_xticklabels(xticklabs)
        else:
            hide_axis(axb,which = 'x')

        axb.set_ylabel('')
        hide_axis(axb,which='y');



        ############################################
        #### Add row label within each track #######
        ############################################
        colnames = obj.colName_short() if shortName else obj.columns
#         colnames = obj.columns
        if len(colnames) > tickMax:
            pass
        else:
            for i,col in enumerate(colnames):
                axa.text(1., i,
                     str(col)[:25],
                     horizontalalignment='right',
                      verticalalignment='center',
                      clip_on=True
                        )
        
        ### Add track name
        trackName = pyutil.formatName(obj.name)
        axa.text(-.0, sum(axa.get_ylim())/2.,
                 trackName,
                 horizontalalignment='right',
                  verticalalignment='center',                 
                )            
        axa.yaxis.set_visible(True);axa.yaxis.set_ticks([])
    
        return axs
    
    def namedIndex(self,index=None,geneName = None):
        index = self.index if index is None else index
        geneName = self.geneName if geneName is None else geneName
        if geneName is None:
            ind = list(index)
        else:
            assert isinstance(geneName,pd.DataFrame)
            res = geneName.reindex(index=index)
            ser = res.iloc[:,0]
            nidx = ser.isnull()
            ser[nidx] = ser.index[nidx]
            ind = list(ser)
#             res.loc[nidx, res.columns[0]] = res.index[nidx]
#             ind = res.iloc[:,0]
        return ind
            
        
    def autoFigsize(self):
        xsize = min(len(self.index)*0.25, 14)
        ysize = (len(self) + 1 )*1.5
        return [xsize,ysize]
        
    def render(self,  figsize= None,
              index=None, how = 'outer', order = None,
               shortName=0,
              silent= 1):
        
        if not self._compiled:
            self.compile(how = how,index = index, order =order)
        else:
            self.joinIndex(how=how) if how is not None else None
            self.orderBy(order=order) if order is not None else None
            
        figsize = figsize or self.figsize or self.autoFigsize()

        fig,axsLst = plt.subplots(len(self),3,figsize=figsize,
#                                   sharex='col',
#                                   frameon = False,
                                  gridspec_kw={'width_ratios':[1, 8, .15],
                                               'wspace':0.05
                                              })
        if len(self)==1:
            axsLst = axsLst[None,:]
#         if not isinstance(axsLst,pyutil.np.ndarray) else axsLst
        
#         plt.suptitle(title,y=0.95)     
        for i in range(len(self)):
            axs, ele = axsLst[i], self[i]            
#             ax = axs[1]
            self._render(ele, axs=axs, look = ele.look,silent = silent,
                        shortName=shortName)
    
        title = 'N=%d'%len(ele); 
        axsLst[0][1].set_title(title)
        
        if self.index.__len__()<=100:
            xticks = self.namedIndex()
            ax = axsLst[-1][1]; plt.sca(ax)
#             ax.xaxis.set_major_formatter(mticker.FixedFormatter())
            plt.xticks(range(len(xticks)),xticks, rotation = 'vertical')
        self.fig['handle'] = fig
        return fig 
    
    def cache(self, cacheName = 'test', fig= None):
        fig = self.fig.get('handle',None) if fig is None else fig
        assert fig is not None
        if self.DIR == './':
            print '[WARN] self.DIR is defaulting to current directory'
        figName = '%s/src/%s.png'%(self.DIR,cacheName)
        self.check_DIR( pyutil.os.path.dirname(figName) )
        figMd = pyutil.showsavefig(fig= fig,fname=figName)
        self.fig['name'], self.fig['md'] = figName,figMd
        return (figName, figMd)
    def copy(self,deep = False):
        new = copy.copy(self)
        return new
    def compile(self,
                index = None, how = 'outer', order = None,
               silent=1):
        '''
index: Specify the data to be included. Inferred from joinIndex() if is None
    how: passed to joinIndex()
    order: passed to joinIndex()
'''
        if self.orig is None:
            self.orig = self.copy()
#         self.joinIndex(how=how) if how is not None else None
        if index is None:
            self.index = self.index if self.index is not None else self.joinIndex(how=how, )
        else:
            self.index = index            
        self.orderBy(order=order) if order is not None else None        
            
        looks = [guessLook(ele) if not ele.look else ele.look for ele in self ]
        
        for i, ele in enumerate(self):
            if not silent:
                print 'Index:%d'% i,
                print 'look:', ele.look, len(ele)
            ele = ele.setDF(ele.reindex( self.index ))
            if looks[i] == 'patch':
                if ele.cmap is None:
                    ele = fixCluster(ele)
            self[i] =  ele
        
        self._compiled = True
        return self
    def recompile(self,*args,**kwargs):
        orig = self.reset()
        orig.compile(*args,**kwargs)
        orig._compiled=False
        return orig
#         return {'index':index}
    
    def reset(self, orig =None):
        orig = orig or self.orig
        if orig is None and not self._compiled:
            print '[WARN] Trying to reset() a raw <panelPlot>'
            orig = self
        return orig
    
    def makeFirst(self,i=0):
        '''
Promote ith element to the starting of list
'''
        self.insert(0,self.pop(i))
        return self
    def joinIndex(self, how = None,):
        '''
order: A DataFrame which will be sorted to obatin the ordering      
'''
        if how is None:
            how = 'outer'
        idx = vstack(self,as_index = 1, how=how)
        self.index = idx
        return idx
    def orderBy(self, order = None, how=None, head = None ):
        if self.index is None:
            idx = self.joinIndex(how=how)
        else:
            idx = self.index
            
        if order is not None:
            df = order.reindex(idx)
            pyutil.reset_columns(df)
            df = df.sort_values(by=list(df.columns),axis=0)
            self.index = idx = df.index
            
        if head is not None:
            idx = self.index[:head]
        self.index = idx
        return self.index
    def set_names(self,names):
        if names:
            it = iter(names)
            for i in range(len(self)):
                name = next(it,None)
                if name:
                    self[i].set_name( name)
        return self
    @property
    def bigTable(self,how = 'outer'):
        lst = [x.copy() for x in self]
        types = pyutil.flatSubset([x.name for x in self],keep='type')
        for i,dd in enumerate(lst):
#             print '[look]=',dd.look
            if dd.look =='patch':
#                 print '[A patch]'
                lst[i] = dd.setDF(dd[dd.columns[-1:]])
    
        for i,dd in enumerate(lst):
            if dd.name.find('type')!=-1:
                TYPE = pyutil.flatSubset([dd.name],keep='type')[0]
            else:
                TYPE = 'unknown'
            dd.columns = ['%s_%s'%(TYPE, col) for col in dd.columns]
            
        res = vstack(lst, how=how)
        return res
    
    def splitTable(self,index=None):
        index = self.index if index is None else index
        lst = [x.reindex(self.index) for x in self]
        return lst
    
    def check_DIR(self,DIR = None):
        DIR = self.DIR if DIR is None else DIR
        pyutil.shellexec('mkdir -p %s'%DIR)
        return DIR
    def emptyExcel(self,fname ='main.xlsx'):
        ExcelFile= pd.ExcelWriter('%s/%s'%(self.DIR,fname), engine='xlsxwriter')
        return ExcelFile
    def dump(self, ExcelFile = 'main.xlsx', sheetName = 'test',figName = None,
             saveClose = 1,
            ):
        figName = figName or self.fig.get('name', None) 
        self.check_DIR()
        if isinstance(ExcelFile, str):
            ExcelFile = self.emptyExcel(ExcelFile)
        self.bigTable.to_excel( ExcelFile, sheetName,index=True,)
        sheet_curr = ExcelFile.sheets[sheetName]
        if figName is not None:
            sheet_curr.insert_image(0, 7, figName)
        if saveClose:
            ExcelFile.save()
            ExcelFile.close()
        return ExcelFile
    
class geneList(scount.countMatrix):
    def __init__(self,C=None,colName=None,rowName=None,**kwargs):
        super(geneList,self).__init__(C=len(rowName) * [1],rowName=rowName, **kwargs)
        self.astype('int')
        

def renderByCluster(panel,clu,DIR = None, addClu = 0,ExcelFile = None):
    assert isinstance(panel,panelPlot)
    assert isinstance(clu, scount.countMatrix)
#     if panel._compiled:
#         panel = panel.reset()
#     if clu not in panel:
    clu = clu.reindex(panel.index)
    if addClu:
        assert 0,"addClu=1 is buggy"
        fixClu = fixCluster(clu)
        fixClu.set_name( 'type=cluster_%s'%fixClu.name)
        panel.append(fixClu)
    if DIR is not None:
        panel.check_DIR(DIR)
        panel.DIR = DIR
    if ExcelFile is not None:
        if isinstance(ExcelFile,str):
            if not ExcelFile.endswith('.xlsx'):
                ExcelFile = ExcelFile +'.xlsx'
            ExcelFile = panel.emptyExcel(ExcelFile) 
        else:
            pass
    it = pyutil.itertools.chain( [(-1, panel)], clu.groupby( clu.columns[0], ) )
    for ci, df in it:
#         print ci
        pp = panel.copy()
        print pp
#         print 'cmap',pp[-1].cmap
        pp.index = df.index
        pp.compile(index = df.index)
        pp.render()
        if ExcelFile is not None:
            sheetName =  'clu%03d'%ci
            pp.cache(sheetName)
            ExcelFile = pp.dump(ExcelFile = ExcelFile,sheetName = sheetName, saveClose=0)
    ExcelFile.save(); ExcelFile.close();
panelPlot.renderByCluster = renderByCluster
#         pp.bigTable()
#         plt.show()