
# coding: utf-8

# In[1]:


get_ipython().run_cell_magic(u'writefile', u'tutil.py', u"import models as mym;reload(mym)\nimport pymisca.vis_util as pyvis\nimport numpy as np\nimport pymisca.util as pyutil; reload(pyutil)\nimport pymisca.vis_util as pyvis; reload(pyvis)\n# if pyutil.is_ipython():\n#     %matplotlib inline\ndef toyData(N0=500,pi=None,K=3,Cs=None,mus=None,D=2):\n    np.random.seed(10)\n#     D = 2\n    lst = []\n    if pi is not None:\n        idx = np.argsort(pi)[::-1]\n    for k in range(K):\n        \n        if pi is not None:\n            k = idx[k]\n            N = int(pi[k]*N0*K)\n        else:\n            N= N0\n        N = max(N,1)\n        \n        if Cs is not None:\n            C = Cs[k]\n        else:\n            C = pyutil.random_covmat(D=D)\n            \n        if mus is not None:\n            mu = mus[k]\n        else:\n            mu = (np.random.random(size = D) - 0.5)*5\n#         if pi is None:\n#         mu = (0,)*D\n        print mu,C,N\n        X = np.random.multivariate_normal(mu, C,size=N)\n        lst += [X]\n    X = np.vstack(lst)\n    return X\n\ndef plotModel(m,X):\n    clu = m.predict(X)\n    print sorted(m.weights_)[::-1]\n    pyvis.qc_2var(X.T[0], X.T[1], clu=clu)\n\n\n    Y = toyData(Cs=   m.covariances_,\n                mus = m.means_,\n                pi  = m.weights_,\n                K=3)\n    pyvis.qc_2var(Y.T[0],Y.T[1],)\n    \n# if __name__!='__main__':\n#     exit()\n#     import sys\n#     sys.exit()\n  ")


# In[ ]:


import models as mym;reload(mym)
import pymisca.vis_util as pyvis
import numpy as np
import pymisca.util as pyutil; reload(pyutil)
import pymisca.vis_util as pyvis; reload(pyvis)
if pyutil.is_ipython():
    get_ipython().magic(u'matplotlib inline')
    
import tutil; reload(tutil)
from tutil import *


# In[ ]:


m = mym.GMM_VIMAP(D=2)

X =  np.random.random(size=(500,2))
m.fit(X)
plotModel(m,X);
m_diag = m



# In[ ]:


X = toyData(K=3)
# print X.shape
pyvis.qc_2var(X.T[0],X.T[1],)


# In[ ]:



m = mym.GMM_VIMAP(D=2,K=3,name='testB')
m.fit(X)
plotModel(m,X)


# In[ ]:


mi = 5
m = mym.GMMLRP_VIMAP(name ='t%d'%mi,D=2,K=3).init_model()
m.fit(X=X)
plotModel(m,X)


# In[ ]:


mi = 6
m = mym.GMMLRP_VIMAP(name ='t%d'%mi,D=2,K=5).init_model()
m.fit(X=X)
plotModel(m,X)


# In[ ]:


# np.sum(np.exp(logP),axis=1)
# mym.tf.nn.log_softmax(logP,axis =1)
logP = m.predict_proba(X)
score = pyutil.logsumexp( logP, axis =1,keepdims=1)
# logPcond = logP - score
# logPcond = 
pyvis.heatmap(logP.T)
pyutil.logsumexp(logP,axis=1,log=0)


# In[ ]:


# m.emission.
if __name__=='__main__':
    import sys
    sys.exit(0)
    get_ipython().system(u'jupyter nbconvert test__gmm.ipynb --to python')


# In[ ]:


# %matplotlib inline
if 0 :
    from __future__ import absolute_import
    from __future__ import division
    from __future__ import print_function

    import edward as ed
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import numpy as np
    import six
    import tensorflow as tf

    from edward.models import (
        Categorical, Dirichlet, Empirical, InverseGamma,
        MultivariateNormalDiag, Normal, ParamMixture)

    plt.style.use('ggplot')

    def build_toy_dataset(N):
        pi = np.array([0.4, 0.6])
        mus = [[1, 1], [-1, -1]]
        stds = [[0.1, 0.1], [0.1, 0.1]]
        x = np.zeros((N, 2), dtype=np.float32)
        for n in range(N):
            k = np.argmax(np.random.multinomial(1, pi))
            x[n, :] = np.random.multivariate_normal(mus[k], np.diag(stds[k]))

        return x


    N = 500  # number of data points
    ed.set_seed(42)

    x_train = build_toy_dataset(N)


    plt.scatter(x_train[:, 0], x_train[:, 1])
    plt.axis([-3, 3, -3, 3])
    plt.title("Simulated dataset")
    plt.show()


# We will use the joint version in this analysis.
