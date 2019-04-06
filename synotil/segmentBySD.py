import hmmlearn.hmm as hlhmm
import synotil.dio as sdio
import pymisca.util as pyutil
import numpy as np
import pandas as pd

def main(bedFile,
#          = None, 
         bwFiles=None,
         bwTrackFile=None,
         relistByGene = 0, #### potentially takes a long time
         stepSize=50,
         
         radius = None,
         center_summit = 0,
         NCORE=1,
    ):
    if bedFile is not None:
        refBed = sdio.extract_peak(bedFile)
    if bwTrackFile is not None:
        dfc =  pyutil.readData(bwTrackFile,)
    else:
        sdio.extract_bigwig_multiple(bedFile=bedFile,
                                     bwFiles=bwFiles,
                                     center_summit=center_summit,
                                    NCORE=NCORE,
                                    stepSize=stepSize,
                                     radius=radius,
                                    )
        relistByGene=1
    if relistByGene:
        dfc =   sdio.listByGene(dfc)
    dfc0 = dfc
    
    tdf =  np.std(dfc.values,axis=1,keepdims=1)
    lr = hlhmm.GaussianHMM(n_components=2, 
                          covariance_type="diag",
                          init_params="cmt", params="cmt")
    lr.fit(tdf)
    seg = lr.predict(tdf)


    segDF  = pd.DataFrame(seg,
                          index = dfc0.index,columns=['clu'])
    ofname = pyutil.getBname(bedFile)+'__HVPeak.bed'
    ofname = sdio.clu2bed(segDF,ofname)
    return ofname

if __name__ == '__main__':
    pass