import models as mym;reload(mym)
import pymisca.vis_util as pyvis
import numpy as np
import pymisca.util as pyutil; reload(pyutil)
import pymisca.vis_util as pyvis; reload(pyvis)
# if pyutil.is_ipython():
#     %matplotlib inline
def toyData(N0=500,pi=None,K=3,Cs=None,mus=None,D=2):
    np.random.seed(10)
#     D = 2
    lst = []
    if pi is not None:
        idx = np.argsort(pi)[::-1]
    for k in range(K):
        
        if pi is not None:
            k = idx[k]
            N = int(pi[k]*N0*K)
        else:
            N= N0
        N = max(N,1)
        
        if Cs is not None:
            C = Cs[k]
        else:
            C = pyutil.random_covmat(D=D)
            
        if mus is not None:
            mu = mus[k]
        else:
            mu = (np.random.random(size = D) - 0.5)*5
#         if pi is None:
#         mu = (0,)*D
        print mu,C,N
        X = np.random.multivariate_normal(mu, C,size=N)
        lst += [X]
    X = np.vstack(lst)
    return X

def plotModel(m,X):
    clu = m.predict(X)
    print sorted(m.weights_)[::-1]
    pyvis.qc_2var(X.T[0], X.T[1], clu=clu)
    _ = m.predict(X)
    Y = m.x_post.sample(len(X)).eval()
#     Y = toyData(Cs=   m.covariances_,
#                 mus = m.means_,
#                 pi  = m.weights_,
#                 K=3)
    pyvis.qc_2var(Y.T[0],Y.T[1],)
    
# if __name__!='__main__':
#     exit()
#     import sys
#     sys.exit()
  