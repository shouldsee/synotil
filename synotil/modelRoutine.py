
import sklearn.mixture as skmix
import pymisca.util as pyutil
np = pyutil.np; pd = pyutil.pd
import os
import qcmsg
import CountMatrix as scount

def fit_BGM(C,
            ALI = 'Test',
#             DIR = ''
#             normF = identityNorm,
            stdPer = 0,
            rowName=None,
            colName=None,
            nClu = 25,
            maxIt = 1000,
            algo = 'DPGMM',
            DIR= '.',
#            algoLst = ['DPGMM'],
            alpha = .1,
            covariance_type = 'diag',
            fixMean = 0,
            reorder=1,
            model_only =0,
            random_state= None,
#             covariance_type = None,
#             **kwargs
           ):
    '''
Fit an BayesianGaussianMixture() model from sklearn
'''
#     if algoLst is None:
#         algoLst = ['DPGMM','DDGMM','GMM',]    
    try:
        DIR,ALI = ALI.rsplit('/',1)
    except:
        DIR= DIR
    os.system('mkdir -p %s'%(DIR))
    
    
    ###### Manage meta attributes of the model ########
    param = {
            'fixMean':fixMean,
             'stdPer':stdPer,
             'nClu':nClu,
             'genre':algo,
             'covarianceType': covariance_type,
              'maxIt' : maxIt,
              'randomState':random_state,
         }
    param.update(getattr(C,'param',{}))
    
    ####### Convert to numpy arrary ######
    if isinstance(C,pd.DataFrame):
        if  ALI=='Test':
            ALI = getattr(C,'name','Test') 

        rowName,colName,C = C.index.values, C.columns, C.values
        pass
    
    ##### Old routine that filter by STD ###########
    if stdPer > 0 :
        assert stdPer < 100,'Percentile must < 100, Got %d instead'%stdPer
        (MEAN,STD,CV),_ = qc_Avg(C)
        pIdx = STD > np.percentile(STD, stdPer)        
        rowName = np.array(rowName)[pIdx]; C = C[pIdx]
    print '[ALI]=',ALI
    nFeat = C.shape[-1]
        
    #####====== Defnitions of fitters=========#######
    
    ###### Arguments shared among fitters ######
    common = {'n_components': nClu,
          'verbose':2,
         'max_iter':maxIt,
             'covariance_type':covariance_type,
              'random_state':random_state,
             }
    if fixMean:
        mean_precision_prior = 1E-128
        mean_prior = [0.]*nFeat
    else:
        mean_precision_prior  = None
        mean_prior = None
        
    ####### List of fitters ######
    mdlLst = {'DPGMM':skmix.BayesianGaussianMixture(weight_concentration_prior_type='dirichlet_process',
                                        weight_concentration_prior=alpha,
                                        mean_precision_prior = mean_precision_prior,
                                        mean_prior = mean_prior,
                                       **common),
          'GMM':skmix.GaussianMixture(**common),
          'DDGMM':skmix.BayesianGaussianMixture(weight_concentration_prior_type='dirichlet_distribution',
                                        weight_concentration_prior=alpha,
                                        mean_precision_prior = mean_precision_prior,
                                        mean_prior = mean_prior,
                                       **common),
         }
    
    
    ############# Select model by "algo"####
    X = C
    print pyutil.qc_matrix(X)
    mdl = mdlLst.get(algo,None)
    assert mdl is not None, 'Algorithm %s not found '%algo
    
    NAME = '%s_%s'%(ALI,pyutil.dict2flat(param))    
    print '[MSG] Now Fitting Model:%s'%NAME
    
    
    
    ####### Meta data of the training Data #######
    d = {'name': NAME,
         'train_data':X,
         'colName':colName,
         'rowName':rowName,
         'param':param,
       }
    
    
    ##### Fitting model and caching the result to specified DIR/NAME ####
    try:
        logFile = open('%s/%s.log'%(DIR,NAME),'w',0)
        with pyutil.RedirectStdStreams(logFile):
            mdl.fixMean= fixMean
            mdl.fit(X)
#             reorderByMSQ(mdl)
            if reorder:
                mdl.reorderByMSQ()
            d.update({'suc':1,'model':mdl})
#             logFile.close()
        print "[SUCC] to fit Model:%s"%(NAME,)
        print qcmsg.msgGMM(mdl)
    except Exception as e:
        print "[FAIL] to fit Model:%s due to :'%s'"%(NAME,e)
        d.update({'suc':0})
    if model_only:
        d['train_data'] = None
        d['rowName'] = None
        d['colName'] = None
    np.save('%s/%s'%(DIR.rstrip('/'),NAME),d)
    d = scount.countMatrix.from_dict(d)
    return d
