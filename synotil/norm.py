import numpy as np
import pymisca.util as pyutil

import sklearn.decomposition as skdecomp
pd = pyutil.pd
import CountMatrix as scount
# import util as util


def stdNorm(X):
    deco = 0
    if isinstance(X,scount.countMatrix):
        deco =1
        df = X 
        
    X = meanNorm(X)
    if isinstance(X, pd.DataFrame):
        X = X.values
    X = X.copy()
    assert isinstance(X,np.ndarray)
    STD = np.std(X,axis=1,keepdims=1); pos = np.squeeze(STD>0);
    X[pos] = X[pos]/STD[pos]

    if deco:
        X = df.setDF(X)
        X.param['normF'] = 'stdNorm'
    return X

def meanNorm(X,axis=1):
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
        X = (X-X.mean(axis=axis,keepdims=1))
        
    param = getattr(X,'param',{})    
    param['normF'] = 'meanNorm'
#     X.param = param
    
    if deco:
#         cols = X.columns
        X = df.setDF(X)
#         X.columns = cols
        X.param['normF'] = 'meanNorm'        
    return X

def sumNorm(X,axis=1):
    deco = 0
    if isinstance(X,scount.countMatrix):
        deco =1
        df = X 
        
    if isinstance(X,pd.Series):
        X = X.values
        SUM = X.sum()
        if SUM != 0:
            X = (X / SUM)
    else:
        if isinstance(X, pd.DataFrame):
            X = X.values
            
        assert isinstance(X,np.ndarray)
        SUM =X.sum(axis=axis,keepdims=1)
        SUM[SUM==0] = 1.
        X = (X / SUM)
        
    param = getattr(X,'param',{})    
    param['normF'] = 'sumNorm'
#     X.param = param
    
    if deco:
#         cols = X.columns
        X = df.setDF(X)
#         X.columns = cols
        X.param['normF'] = 'sumNorm'        
    return X

def diffNorm(X,axis=1):
    deco = 0
    if isinstance(X,scount.countMatrix):
        deco =1
        df = X 
        
    if isinstance(X,pd.Series):        
        X = X.values
        X = np.diff(X,axis=0)
    else:
        if isinstance(X, pd.DataFrame):
            X = X.values            
        assert isinstance(X,np.ndarray)
        X = np.diff(X,axis=axis)
#         X = (X-X.mean(axis=axis,keepdims=1))
        
    param = getattr(X,'param',{})    
    param['normF'] = 'diffNorm'
#     X.param = param
    
    if deco:
#         cols = X.columns
        X = df.setDF(X)
#         X.columns = cols
        X.param['normF'] = 'diffNorm'        
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
def meanNormPCA(X,cutoff=1.0, withModel=0):
    name = 'meanNormPCA'
    deco = 0
    if isinstance(X,scount.countMatrix):
        deco =1
        df = X
        
    X = meanNorm(X)
    X = X.values if isinstance(X, pd.DataFrame) else X
    
    X0 = X
    mdl = skdecomp.PCA()
    mdl.fit(X)
    
    X = mdl.transform(X)
    
    first = np.cumsum(mdl.explained_variance_ratio_) >= cutoff
    idx   =np.argmax(first) + 1
    mdl.components_ = mdl.components_[:idx]
#     mdl.mean_ = mdl.mean_[:,idx]

    X = X[:,:idx]
        
    param = getattr(X,'param',{})    
    param['normF'] = name
#     X.param = param

    if deco:
        X = df.setDF(X)
        X.param['normF'] = name     
    if withModel:
        return (X,mdl,X0)
    else:
        return X
def ctNorm(X):
    X = (X-X.mean(axis=1,keepdims=1))
    return X
def identityNorm(X):
    param = getattr(X,'param',{})    
    param['normF'] = 'identityNorm'
#     X.param = param    
    return X
