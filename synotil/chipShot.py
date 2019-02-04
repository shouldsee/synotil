# %load chipShot.py
#!/usr/bin/env python2
import matplotlib as mpl
mpl.use('agg')
import pymisca.util as pyutil
import synotil.CountMatrix as scount
# import synotil.PanelPlot as spanel
import synotil.util as sutil
import synotil.filterByCDS
import synotil.dio as sdio
import synotil.jobs as sjob
# spanel.plt.ioff()
import time
stderrLine = lambda x:pyutil.sys.stderr.write(x+'\n')


def query__gtf(gtf,accs):
    val = gtf[gtf.columns[-1]]
    qres = val.str.contains('|'.join(accs))
    res =  gtf.reindex(gtf.index[qres])
    return res

def prepare_chipTrack(x,vlim = None,**kwargs):
    ind  = x.columns
    x.columns = x.columns.droplevel(0) 
    x.name = ind.levels[0][0]
    y  = scount.countMatrix(x,look='fill',name = x.name,vlim=vlim,**kwargs) 
    
    return y
# reload(pybio)
# reload(spanel)
# figsize = 
# def worker__drawPeak(peakAcc,DIR=None,drawPeak=None):
#     ofname = '%s/%s.svg'%(DIR,peakAcc)
#     fig = drawPeak(peakAcc)
#     fig.savefig(ofname)
#     return ofname

def worker__drawPeak(peakAcc,
                     df_near, 
                     gtfs,
                     chipTracks,
                     radius,
                     figsize=None,
                     DIR=None,
                     debug=0,
                     ylim = None,
                    ):
    '''For parallel running'''
    import synotil.PanelPlot as spanel
#     def worker():
    if 1:
        df_feats = df_near.query('acc=="%s"'%peakAcc)
        if df_feats.empty:
            start = 0
            end = radius * 2
#             return None
        else:
            first = df_feats.iloc[0]
            start= first.start - radius
            end = first.end + radius

        def prepare__gtf(gtf,):
            if df_feats.empty:
                return None
            else:
                gcurr = query__gtf(gtf, accs = df_feats.feat_acc.values)
                #### shift to appropriate coordinate for gtf
                gcurr[3]  += -start   
                gcurr[4]  = gcurr[4] -start 
            #     gcurr = gcurr.reindex(gcurr.index[(gcurr[3]>0) & (gcurr[4]>0)])
                gcurr = gcurr.reindex(gcurr.index[(gcurr[4]>0)])
                return gcurr

        tracks = chipTracks +[ x for x in map(prepare__gtf,gtfs) if x is not None]
        if debug:
            val = prepare__gtf(gtfs[0])
            print('[GCURR]')
            if val is None:
                print (val)
            else:
                print(val.look)
                print (pyutil.mat2str(val.values[:,:8]))
        
        # gcurr = gcurr.reindex(gcurr.index[(gcurr[3]>0) & (gcurr[4]>0)])
        span = 2*radius
        title = '''\
{peakAcc}\nstart={start}  end={end}, span={span}
ylim={ylim}
'''.format(**locals())
        ofname = '%s/%s.svg'%(DIR,peakAcc)
        while True:
            try:
                pp = spanel.panelPlot(tracks ,
                                      index=[peakAcc],
                                      figsize=figsize,
                                      debug=debug, 
                                      xlim = [0,radius * 2],
                                      show_axa=False, 
            #                           vlim = ylim,
                                         title=title)
                fig = pp.render(silent=1);
                fig.savefig(ofname)
                break
            except RuntimeError as e:
                print ('[WARN] caught exception when plotting:%s'%e)
                time.sleep(0.05)
#         fig = pp.render(silent=1);

        print('[MSG] plotting to: %s' % ofname)

        return ofname
    
def worker__fluff(rec,):
    rec = pyutil.util_obj(**rec)
    DIR = getattr(rec,'DIR','.')
    ext = getattr(rec,'ext','svg')
    labels = getattr(rec,'labels',None)
#     ofname  = rec.acc + '.svg'
    ofname  = '%s/%s.%s' % (DIR,rec.acc,ext) 
    interval = rec.interval
    tracks = rec.tracks
    annotation = rec.annotation
#             ofname   = bed.acc[i]  + '.svg'
#             interval = bed.interval[i]
    ofname   = sjob.fig__fluffProfile(
                        interval,
                        tracks, 
                        ofname = ofname,
                        annotation = annotation,
                        labels = labels)
    return ofname        
# worker = pyutil.functools.partial(sjob.fig__fluffyProfile,
# #                                          interval = 
#                                  )            
def main(
    #### necessary
    bedFile = None,
    bwFiles = None,
    ####  
    DIR=None,
    figsize=[14,14],
    debug = 0,
    ylim = [0,10],
    radius = 2000,
    stepSize = 10,
    NCORE= 4,
    silent = 0,
    gtfFile = None,
    cdsFile = None,
    annotation = None,
    GSIZE= None,
    center_summit = 0,
    trackNames = None,
    backend = 'fluff',
    ext = 'png',
    **kwargs
):
#     vlim = ylim
    figsize= map(int,figsize)
    # for peakAcc in df_near.acc.unique()[:1]:   

    prefix = 'PROG=chipShots_bedFile='
#     bname  = pyutil.basename(bedFile)
    bname = pyutil.os.path.basename(bedFile)
    odname = prefix + bname 
    if DIR=='inplace':
        DIR = pyutil.os.path.dirname(bedFile) + odname
    elif DIR is None:
        DIR = odname
    pyutil.shellexec('mkdir -p %s' % DIR,silent=silent)
    DIR = pyutil.os.path.abspath(DIR)
#     odname = pyutil.ospath

    if cdsFile is None:
        cdsFile = gtfFile + '.cds'
    if backend == 'synotil':
        # nearFile = './DE2017/type=closest_bed=lux22_radius=1_feat=genes.gtf.cds.tsv'
        # import synotil.filterByCDS
        nearFile = synotil.filterByCDS.main(
            peakFile=bedFile,
            cdsFile=cdsFile,
            downStream=radius,
            upStream=radius,
            peakRadius = 1,
            GSIZE=GSIZE,
            center_summit=center_summit,
            )
        df_near = pyutil.readData(nearFile,)

        stderrLine('[MSG]Loading bed intervals from bigwig tracks....')

        chipTracks = sutil.extract_bigwig_multiple(fnames=bwFiles,
                                                  bedFile=bedFile,
                                                  radius=radius,
                                                  stepSize=stepSize,
                                                   callback=None,
                                                   outIndex=trackNames,

    #                                               callback=callback,
                                                   center_summit=center_summit,
                                                   shift = 0, #### use positive coordinate
                                                   stranded=False,
                                                   NCORE=NCORE)
        if ylim is None:
            ylim = pyutil.span(pyutil.np.hstack([x.values.flat for x in chipTracks]),99)
            ylim = list(ylim)
            ylim[0] = 0.
        callback = lambda x:[prepare_chipTrack(ele,
                                               vlim=ylim) 
                             for ele in x] 
        chipTracks = callback(chipTracks)

        if debug:
            stderrLine(chipTracks[0].columns)

        gtf = pyutil.readData(gtfFile,ext='tsv',header=None,guess_index=0)
        gtf = scount.countMatrix(gtf,look='gtf')
        gtfs= [gtf]


    #     uniqPeak = df_near.acc.unique()
    #     bedDF = pyutil.readData(bedFile,header=None,guess_index=0)
    #     bedDF.columns = sutil.bedHeader[:len(bedDF.columns)]
        bedDF = sutil.extract_peak(bedFile)
    #     uniqPeak
    #     uniqPeak = bedDF[bedDF.columns]

        worker=  pyutil.functools.partial(worker__drawPeak, 
                                          DIR=DIR,
                                          chipTracks = chipTracks,
                                          df_near=df_near,
                                          gtfs = gtfs,
                                          radius = radius,
                                          figsize=figsize,
                                          ylim = ylim,
                                          debug=debug,)

        ofnames = pyutil.mp_map(worker, bedDF.acc, n_cpu=NCORE,)
    elif backend=='fluff':
        bedDF = sdio.extract_peak(bedFile)
        
        argDF = bedDF.copy()
        argDF = sdio.bed__addCol__interval(argDF)
        tracks = list(bwFiles)
        argDF['tracks'] = [tracks ] * len(bedDF)
        argDF['annotation'] = annotation
        argDF['DIR'] = DIR
        argDF['ext'] = ext
        if trackNames is not None:
            argDF['labels'] = [list(trackNames)]  * len(bedDF)
            
        ofnames = pyutil.mp_map(
#         ofnames = map(
            worker__fluff, 
            ( vars(x) for x in argDF.itertuples() ), 
            n_cpu=NCORE,
         )
#         ofnames = 
        
    
    bedDF['img'] = ofnames
    indexFile = '%s/%s.index.tsv'%(DIR,bname)
    pyutil.to_tsv(bedDF,indexFile)
    indexFile = '%s/figureIndex.tsv'%(DIR)
    pyutil.to_tsv(bedDF,indexFile)
    
    try:
        import synotil.shot2html as shot2html
        htmlFile = shot2html.shot2html(indexFile,localPath=True)
    except Exception as e:
        stderrLine ('[WARN]:cannot produce html :%s'%e)
        htmlFile = None
#     print ('[OUTPUT]:',)
#     print ('html:',htmlFile)
#     print ('index:',indexFile)
    print (indexFile) 
    print (htmlFile)
    return (indexFile,htmlFile)


testDict = {
    'bedFile': './DE2017/lux22_radius=1.tsv',
    'bwFiles': ['/home/feng/envs/Fig_Brachy/DE2017/bw/122C/EC-1167_S5/EC-1167_S5_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_genomenorm.bw',
     '/home/feng/envs/Fig_Brachy/DE2017/bw/122C/EC-1167_S6/EC-1167_S6_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_genomenorm.bw',
     '/home/feng/envs/Fig_Brachy/DE2017/bw/122C/EC-1167_S7/EC-1167_S7_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_genomenorm.bw',
     '/home/feng/envs/Fig_Brachy/DE2017/bw/122C/EC-1167_S8/EC-1167_S8_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_genomenorm.bw'],
    'gtfFile': '/home/feng/ref/Arabidopsis_thaliana_TAIR10/annotation/genes.gtf',
    'GSIZE': '/home/feng/ref/Arabidopsis_thaliana_TAIR10/genome.sizes',
    'NCORE' :10,
    'radius':2000,
}


import argparse

parser= argparse.ArgumentParser()
parser.add_argument('-v','--silent',default=0,type = int)

parser.add_argument('-b','--bedFile',default=None,)
parser.add_argument('bwFiles',nargs='+')
parser.add_argument('-r','--radius',default=2000, type=int,)

# parser.add_argument('-c','--cdsFile',default=pyutil.os.environ.get('GTF','none')+'.cds')
parser.add_argument('-a','--gtfFile',
                    default=pyutil.os.environ.get('GTF','none')+'.cds')
parser.add_argument('-c','--cdsFile',
                    default=None)

parser.add_argument('-g','--GSIZE',
                    default=pyutil.os.environ.get('GSIZE',None))
parser.add_argument('-o','--DIR',default=None)
parser.add_argument('-j','--NCORE',
                    default=pyutil.os.environ.get(4,None),type=int)
parser.add_argument('-s','--center_summit',
                    default=0, type=int)
parser.add_argument('-d','--debug',
                    default=0, type=int)
parser.add_argument('-f','--figsize',
                    default=[14,14], type=int,nargs=2)
# parser.add_argument('-y','--ylim',
#                     default=[0., 10.], type=float,nargs=2)
parser.add_argument('-y','--ylim',
                    default=None, type=float,nargs=2)

defaults = {act.dest:act.default for act in parser._actions}
for key in ['bedFile','bwFiles']:
    defaults.pop(key)

main = pyutil.functools.partial(main,**defaults)
# argparser
if __name__=='__main__':
    args = parser.parse_args()
    pars = vars(args)
    assert pars['bedFile'] is not None
    print (pyutil.ppJson(pars))
    main(**pars)
    stderrLine('[Done]')