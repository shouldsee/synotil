import numpy as np
import pymisca.util as pyutil

import sklearn.decomposition as skdecomp
pd = pyutil.pd
import CountMatrix as scount
# import util as util


def stdNorm(X):
    X = meanNorm(X).copy()
    STD = np.std(X,axis=1,keepdims=1); pos = np.squeeze(STD>0);
    X[pos] = X[pos]/STD[pos]

    param = getattr(X,'param',{})    
    param['normF'] = 'stdNorm'
#     X.param = param
    
    return X
def meanNorm(X):
    deco = 0
    if isinstance(X,scount.countMatrix):
        deco =1
        df = X 
        
    if isinstance(X,pd.Series):
        X = X.values
        X = (X - X.mean())
    else:
        if isinstance(X, pd.DataFrame):
            X = X.values
        assert isinstance(X,np.ndarray)
        X = (X-X.mean(axis=1,keepdims=1))
        
    param = getattr(X,'param',{})    
    param['normF'] = 'meanNorm'
#     X.param = param
    
    if deco:
        X = df.setDF(X)
        X.param['normF'] = 'meanNorm'        
    return X

def meanNormProj(X):
    deco = 0
    if isinstance(X,scount.countMatrix):
        deco =1
        df = X
        
    X = meanNorm(X)
    if isinstance(X, pd.DataFrame):
        X = X.values    
        
    D = X.shape[1]
    Wn  = pyutil.meanNormBasis(D,orthonormal =1 )
    X = X.dot(Wn.T)
    
    param = getattr(X,'param',{})    
    param['normF'] = 'meanNormProj'
#     X.param = param

    if deco:
        X = df.setDF(X)
        X.param['normF'] = 'meanNormProj'        
    return X
import sklearn.decomposition as skdecomp
def meanNormPCA(X,cutoff=1.0):
    name = 'meanNormPCA'
    deco = 0
    if isinstance(X,scount.countMatrix):
        deco =1
        df = X
        
    X = meanNorm(X)
    X = X.values if isinstance(X, pd.DataFrame) else X
    
    mdl = skdecomp.PCA()
    mdl.fit(X)
    
    X = mdl.transform(X)
    
    idx = np.cumsum(mdl.explained_variance_ratio_) >= cutoff
    X = X[:,:np.argmax(idx) + 1]
        
    param = getattr(X,'param',{})    
    param['normF'] = name
#     X.param = param

    if deco:
        X = df.setDF(X)
        X.param['normF'] = name        
    return X
def ctNorm(X):
    X = (X-X.mean(axis=1,keepdims=1))
    return X
def identityNorm(X):
    param = getattr(X,'param',{})    
    param['normF'] = 'identityNorm'
#     X.param = param    
    return X
