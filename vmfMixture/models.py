# import edward as ed
# import edward.models as edm

import tensorflow_probability.python.distributions as tfdist
# import tensorflow_probability.python.edward2 as tfed
import pymisca.tensorflow_extra as pytf
ed  = edm = pytf.ed; 


import tensorflow as tf
import numpy as np
import pymisca.numpy_extra as pynp


class BaseModel(object):
    def __init__(self,name='test'):
        self.name = name
        print type(self),self.__dict__
        pass
    def sanitise(self,X):
        X  = np.asarray(X,np.float32)        
        return X
    def fit(self,X,n_iter = 1000, n_print=100,
            **kwargs):
        X = self.sanitise(X)
        self._fit(X,n_iter = n_iter,n_print = n_print,
                  **kwargs)
        pass
    
    def predict_proba(self,X,norm = 1,log=1):
        X = self.sanitise(X)
        prob = self._predict_proba(X)
        if norm:
            prob = prob - pynp.logsumexp( prob, axis =1,keepdims=1)
        if not log:
            prob = np.exp(log)
        return prob
    def score(self,X,keepdims=0):
        prob = self.predict_proba(X,norm=0,log=1)
        score = pynp.logsumexp( prob, axis =1,keepdims=keepdims)
        return score
    
    def predict(self,X):
        proba = self.predict_proba(X,norm = 0)
        clu = np.argmax(proba,axis = 1)
        return clu
    
    def expand_input(self,X):
        N = len(X)
        X_bdc = tf.tile(
            tf.reshape(
                X, [N, 1, 1, self.D]),
                   [1, 1, self.K, 1])
        return X_bdc
    
    
class ModelA(BaseModel):
    def __init__(self,):
        super(ModelA,self).__init__()
        print self,self.__dict__
        pass
    def fit(self,X):
        pass
m = ModelA()
m.fit(1)# from models import *
from gamma_radial_theta import *
from gamma_radial_affine_lrp import *



class GMM_VIMAP(BaseModel):
    def __init__(self,D=None,K=20,*args,**kwargs):
        super(GMM_VIMAP,self).__init__(*args,**kwargs)
        self.K = K
        self.D = D
        self.N_post = None
        self.emission = None
        self.prior = None
        if D is not None:
            self.init_param(D=D)    
        
        print self,self.__dict__
#         pass
    def init_param(self,D=None,K = None,
#                    N=None
                  ):
        self.D = D = self.D if D is None else D
        assert D is not None
        self.K = K = self.K if K is None else K
        assert K is not None
  
        self.pi = edm.Dirichlet(tf.ones(K))
    
        self.mu = edm.Normal(tf.zeros(D), tf.ones(D), sample_shape=K)
        self.sigmasq = edm.InverseGamma(tf.ones(D), tf.ones(D), sample_shape=K)
#         name = self.name if self.'test'
        name = self.name

        try:
            tf.get_variable(name+'/test',[1])
            reuse = None
        except:
            reuse = True
        print reuse
        
        with tf.variable_scope(name, reuse=reuse):
            self.qmu = ed.models.PointMass(
            tf.get_variable("qmu", [K,D])
        )
            self.qsig = ed.models.PointMass(
            tf.nn.softplus(tf.get_variable("qsig", [K,D]))
        )
            self.q_pi = ed.models.PointMass(
                tf.nn.softmax(
                    tf.get_variable("q_pi", [K])
                )
            )
        self.pDict = {self.mu: self.qmu, self.sigmasq: self.qsig,
                      self.pi:self.q_pi}        
        
        return self
    @property
    def means_(self):
        res = self.qmu.eval()
        return res
    @property
    def covariances_(self):
        res = map(np.diag,self.qsig.eval())
        return res
    @property
    def weights_(self):
        res = self.q_pi.eval()
        return res    
    def _fit(self,X):
        X  = np.asarray(X,np.float32)        
        K = self.K
        N = len(X)
        
        self.cat = edm.Categorical(probs=self.pi, sample_shape=N)
        self.components = [
            edm.MultivariateNormalDiag(
                self.mu[k],
                tf.sqrt(self.sigmasq[k]
                                          ), sample_shape=N)
            for k in range(K)]
        self.emission = edm.Mixture(
            cat=self.cat, 
            components=self.components,
            sample_shape=N)
        
        self.dDict = {self.emission: X}
        self.inference = ed.MAP(self.pDict,self.dDict)
        self.inference.run(n_iter=500, n_print=100, )
        return self
    
    
    def _predict_proba(self,X, N=None, norm = 0):
        ''' self.emission does not work, manually build posterior
'''
        if X is None:
            assert N is not None
        else:
            N = len(X)
        if not self.N_post is N:        
            ones = tf.ones([N, 1, 1, 1])
            self.x_post = em = edm.Normal(
                    loc = ones * self.qmu,
                  scale = ones * self.qsig)
            self.N_post = N
        #     _make_posterior(self,X=X)
#         self._make_posterior(X=X)
        N = len(X)
        #### (num_sample, num_post, num_component, dimensionality)
        X_bdc = tf.tile(
            tf.reshape(
                X, [N, 1, 1, self.D]),
                   [1, 1, self.K, 1])

        X_loglik = ll = self.x_post.log_prob(X_bdc)
        ll = tf.reduce_mean(ll,axis=1)  ### over posterior samples
        ll = tf.reduce_sum(ll,axis=-1)  ### over dimensions
        logP = ll.eval()
        return logP
    
#     def predict(self,X):
#         logP = self.predict_proba(X)
#         clu  = logP.argmax(axis=1)
#         return clu# name = 'test6'
# reload(mym)
# from models import *
class GMMLRP_VIMAP(BaseModel):
#     self.emDist = edm.MultivariateNormalDiagPlusLowRank    
#     emDist = edm.MultivariateNormalDiagPlusLowRank    
    def emDist(self,*args,**kwargs):
        emd = edm.as_random_variable(
            tfdist.MultivariateNormalDiagPlusLowRank(*args,**kwargs)
        )
        return emd
#     emDist = 
    def __init__(self,D=None,K=20,
                 *args,**kwargs):
        super(GMMLRP_VIMAP,self).__init__(*args,**kwargs)
        self.K = K
        self.D = D
        self.initialised = False
        if D is not None:
            self.init_model(D=D)   
        
            
    def init_model(self,D=None,K = None,alpha = 1.0):
        self.D = D = self.D if D is None else D
        assert D is not None
        self.K = K = self.K if K is None else K
        assert K is not None
        name = self.name
        try:
            tf.get_variable(name+'/test',[1])
            reuse = None
        except:
            reuse = True
        print reuse
        with tf.variable_scope(name, reuse=reuse):
            
            ##### Prior
            mu = edm.Normal(tf.zeros(D), tf.ones(D), sample_shape=K)            
            scale_diag =  edm.Uniform(sample_shape=(K,D))
            scale_perturb_factor =  edm.Uniform(sample_shape=(K,D,1))            
#             scale_perturb_factor = edm.Normal(
#                     loc=tf.zeros(1),
#                     scale=tf.ones(1),
#                     sample_shape=(K,D)
#                 )
            
#             self.pi = edm.Dirichlet(tf.ones(K))
            self.pi = pi = edm.Dirichlet( float(alpha)/K * tf.ones(K) )
            
        
            ##### Posterior
            self.q_pi = ed.models.PointMass(
                tf.nn.softmax(
                    tf.get_variable("q_pi", [K])
                )
            )
            self.q_mu = ed.models.PointMass(
                tf.get_variable("q_mu", [K,D])
            )
            self.q_scale_diag  = edm.PointMass(
                tf.nn.softplus(
                    tf.get_variable('q_scale_diag',shape=[K,D])
                              ),
            )
            self.q_scale_perturb_factor = ed.models.PointMass(
                (
                    tf.get_variable("q_scale_perturb_factor", [K,D,1])
                )
            )
            
    
        ##### Dictonary for constructing self.emDist(**self.param)
        self.param = {'loc':(mu,
                             self.q_mu),
                     'scale_diag':(scale_diag,
                                   self.q_scale_diag),
                     'scale_perturb_factor':(scale_perturb_factor,
                                             self.q_scale_perturb_factor
                                             ),
                     'weight':(self.pi,
                               self.q_pi),
                     }
        self.emKey = ['loc','scale_diag','scale_perturb_factor']
        self.priorDict = {v[0]:v[1] for v in self.param.values()}
#         self.priorDict.update({self.pi:self.q_pi})
        
        self.postDict = {k:v[1] for k,v in self.param.items()}
        
        ### Prior components
        cDicts = [
            {key:v[0][k] 
             for key,v in self.param.items() 
             if key in self.emKey} 
            for k in range(K)]
        self.components = [self.emDist(**d) for d in cDicts]
        
        ### Posterior generative
        self.x_post = em = self.emDist(**{k:v for k,v in self.postDict.items()
                                         if k in self.emKey})

        self.initialised = True; return self
            
    def _fit(self,X,n_iter=1000, n_print=100, **kwargs):
#         X  = np.asarray(X,np.float32)        
        K = self.K
        N = len(X)

        self.cat = edm.Categorical(probs=self.pi, 
                                   sample_shape=N
                                  )
        self.emission = edm.Mixture(
            cat=self.cat, 
            components=self.components,
            sample_shape=N
        )
        
        
        
        print ('hiModel')
        self.dDict = {self.emission: X}
        self.inference = ed.MAP(self.priorDict,self.dDict)
        self.inference.run(n_iter = n_iter, n_print=n_print,*kwargs)
        return self
    
    @property
    def means_(self):
        res = self.x_post.mean().eval()
        return res
    @property
    def covariances_(self):
        res = self.x_post.covariance().eval()
        return res
    @property
    def weights_(self):
        res = self.q_pi.eval()
        return res
    def _predict_proba(self,X, N=None, norm = 0):
        ''' self.emission does not work, manually build posterior
'''
        
        N = len(X)
        X_bdc = self.expand_input(X)
        # em.log_prob(X_bdc).eval()
        X_loglik = ll = self.x_post.log_prob(X_bdc)
        ll = tf.reduce_mean(ll,axis=1)  ### over posterior samples
    #     ll = tf.reduce_sum(ll,axis=-1)  ### over dimensions
        logP = ll.eval()   
        return logP
    
class GLRP_VIMAP(BaseModel):
#     self.emDist = edm.MultivariateNormalDiagPlusLowRank    
    def emDist(self,*args,**kwargs):
        emd = edm.as_random_variable(
            tfdist.MultivariateNormalDiagPlusLowRank(*args,**kwargs)
        )
        return emd

    def __init__(self,D=None,K=20,
                 *args,**kwargs):
        super(GMMLRP_VIMAP,self).__init__(*args,**kwargs)
        self.K = K
        self.D = D
        self.initialised = False
        if D is not None:
            self.init_model(D=D)   
        
            
    def init_model(self,D=None,K = None,alpha = 1.0):
        self.D = D = self.D if D is None else D
        assert D is not None
        self.K = K = self.K if K is None else K
        assert K is not None
        name = self.name
        try:
            tf.get_variable(name+'/test',[1])
            reuse = None
        except:
            reuse = True
        print reuse
        with tf.variable_scope(name, reuse=reuse):
            
            ##### Prior
            mu = edm.Normal(tf.zeros(D), tf.ones(D), sample_shape=K)            
            scale_diag =  edm.Uniform(sample_shape=(K,D))
            scale_perturb_factor =  edm.Uniform(sample_shape=(K,D,1))            
#             scale_perturb_factor = edm.Normal(
#                     loc=tf.zeros(1),
#                     scale=tf.ones(1),
#                     sample_shape=(K,D)
#                 )
            
#             self.pi = edm.Dirichlet(tf.ones(K))
#             self.pi = pi = edm.Dirichlet( float(alpha)/K * tf.ones(K) )
            
        
            ##### Posterior
            self.q_pi = ed.models.PointMass(
                tf.nn.softmax(
                    tf.get_variable("q_pi", [K])
                )
            )
            self.q_mu = ed.models.PointMass(
                tf.get_variable("q_mu", [K,D])
            )
            self.q_scale_diag  = edm.PointMass(
                tf.nn.softplus(
                    tf.get_variable('q_scale_diag',shape=[K,D])
                              ),
            )
            self.q_scale_perturb_factor = ed.models.PointMass(
                (
                    tf.get_variable("q_scale_perturb_factor", [K,D,1])
                )
            )
            
    
        ##### Dictonary for constructing self.emDist(**self.param)
        self.param = {'loc':(mu,
                             self.q_mu),
                     'scale_diag':(scale_diag,
                                   self.q_scale_diag),
                     'scale_perturb_factor':(scale_perturb_factor,
                                             self.q_scale_perturb_factor
                                             ),
                     'weight':(self.pi,
                               self.q_pi),
                     }
        
        self.emKey = ['loc','scale_diag','scale_perturb_factor']
        self.priorDict = {v[0]:v[1] for v in self.param.values()}
#         self.priorDict.update({self.pi:self.q_pi})
        
        self.postDict = {k:v[1] for k,v in self.param.items()}
        
        ### Prior components
        cDicts = [
            {key:v[0][k] 
             for key,v in self.param.items() 
             if key in self.emKey} 
            for k in range(K)]
        self.components = [self.emDist(**d) for d in cDicts]
        
        ### Posterior generative
        self.x_post = em = self.emDist(**{k:v for k,v in self.postDict.items()
                                         if k in self.emKey})

        self.initialised = True; return self
            
    def _fit(self,X,n_iter=1000, n_print=100, **kwargs):
#         X  = np.asarray(X,np.float32)        
        K = self.K
        N = len(X)

        self.cat = edm.Categorical(probs=self.pi, 
                                   sample_shape=N
                                  )
        self.emission = edm.Mixture(
            cat=self.cat, 
            components=self.components,
            sample_shape=N
        )
        
        
        
        print ('hiModel')
        self.dDict = {self.emission: X}
        self.inference = ed.MAP(self.priorDict,self.dDict)
        self.inference.run(n_iter = n_iter, n_print=n_print,*kwargs)
        return self
    
    @property
    def means_(self):
        res = self.x_post.mean().eval()
        return res
    @property
    def covariances_(self):
        res = self.x_post.covariance().eval()
        return res
    @property
    def weights_(self):
        res = self.q_pi.eval()
        return res
    def _predict_proba(self,X, N=None, norm = 0):
        ''' self.emission does not work, manually build posterior
'''
        
        N = len(X)
        X_bdc = self.expand_input(X)
        # em.log_prob(X_bdc).eval()
        X_loglik = ll = self.x_post.log_prob(X_bdc)
        ll = tf.reduce_mean(ll,axis=1)  ### over posterior samples
    #     ll = tf.reduce_sum(ll,axis=-1)  ### over dimensions
        logP = ll.eval()   
        return logP

    
    
# import models; reload(models); from models import *
import pymisca.tensorflow_extra as pytf; reload(pytf)
import pymisca.util as pyutil

class SunifLRP_VIMAP(BaseModel):
#     self.emDist = edm.MultivariateNormalDiagPlusLowRank    
#     emDist = edm.MultivariateNormalDiagPlusLowRank    
    def emDist(self,*args,**kwargs):
        emd = edm.as_random_variable(
            tfdist.sphereUniformLRP(*args,**kwargs)
        )
        return emd
    
    def __init__(self,D=None,K=20,
                 *args,**kwargs):
        super(SunifLRP_VIMAP,self).__init__(*args,**kwargs)
        self.K = K
        self.D = D
        self.initialised = False
        if D is not None:
            self.init_model(D=D)   
        
            
    def init_model(self,D=None,K = None,alpha = 1.0):
        self.D = D = self.D if D is None else D
        assert D is not None
        self.K = K = self.K if K is None else K
        assert K is not None
#         print (K)
        
        uspan = [-1E5,1E5]
        name = self.name
        try:
            tf.get_variable(name+'/test',[1])
            reuse = None
        except:
            reuse = True
        print reuse
        prior = pyutil.util_obj()
        post = pyutil.util_obj()
        with tf.variable_scope(name, reuse=reuse):
            
            ##### Prior
            prior.mu = edm.Normal(tf.zeros(D), tf.ones(D), sample_shape=K)            
            prior.scale_diag =  edm.Uniform(*uspan,sample_shape=(K,D))
            prior.scale_perturb_factor =  edm.Uniform(*uspan,sample_shape=(K,D,1))   
            prior.concentration =  edm.Uniform(*uspan,sample_shape=(K,1))
            prior.rate =  edm.Uniform(*uspan,sample_shape=(K,1))
#             prio
#             scale_perturb_factor = edm.Normal(
#                     loc=tf.zeros(1),
#                     scale=tf.ones(1),
#                     sample_shape=(K,D)
#                 )
#             prior.weight = edm.Dirichlet(tf.ones(K))
            prior.weight = pi = edm.Dirichlet( float(alpha)/K * tf.ones(K) )
        
            ##### Posterior
            post.weight = ed.models.PointMass(
                tf.nn.softmax(
                    tf.get_variable("q_pi", [K])
                )
            )
            post.mu = ed.models.PointMass(
                tf.get_variable("q_mu", [K,D])
            )
            
            post.scale_diag  = edm.PointMass(
                tf.nn.softplus(
                    tf.get_variable('q_scale_diag',shape=[K,D])
                              ),
            )
            
            post.scale_perturb_factor = ed.models.PointMass(
                (
                    tf.get_variable("q_scale_perturb_factor", [K,D,1])
                )
            )
            post.concentration  = edm.PointMass(
                tf.nn.softplus(
                    tf.get_variable('concentration',shape=[K,1])
                              ),
            )
            post.rate  = edm.PointMass(
                tf.nn.softplus(
                    tf.get_variable('rate',shape=[K,1])
                              ),
            )
            
            
        self.prior = prior
        self.post = post
        
        ##### Dictonary for constructing self.emDist(**self.param)
        self.em_key =[
            'scale_diag',
            'scale_perturb_factor',
            'concentration',
            'rate',
        ]
        self.mix_key = [
            'weight',
        ]
        self.param_key = (self.em_key + 
                          self.mix_key)


#         self.emKey = ['loc','scale_diag','scale_perturb_factor']
        self.paramDict = {getattr(prior,name):
                          getattr(post,name) for name in self.param_key}
#         self.paramDict = {}
#         self.priorDict = {v[0]:v[1] for v in self.param.values()}
#         self.priorDict.update({self.pi:self.q_pi})
        
#         self.postDict = {k:v[1] for k,v in self.param.items()}
        
        ### Prior components
        cDicts = [
            {key: v[k] 
             for key,v in prior.__dict__.items() 
             if key in self.em_key} 
            for k in range(K)]
        self.components = [self.emDist(**d) for d in cDicts]
        
        ### Posterior generative
#         edm.Mixture
        cDicts = [
            {key: v[k] 
             for key,v in post.__dict__.items() 
             if key in self.em_key} 
            for k in range(K)]
        self.postComponents = [self.emDist(**d) for d in cDicts]
        

        
        
#         edm.ParamMixture
#         self.x_post = em = self.emDist(**{k:v for k,v in self.post.__dict__.items()
#                                          if k in self.em_key})

        self.initialised = True; return self
            
    def _fit(self,X, n_iter=1000, n_print=100, **kwargs):
#         X  = np.asarray(X,np.float32)        
        K = self.K
        N = len(X)

        self.cat = edm.Categorical(probs=self.prior.weight, 
                                   sample_shape=N
                                  )
        self.emission = edm.Mixture(
            cat=self.cat, 
            components=self.components,
            sample_shape=N
        )
        
        
        
        print ('hiModel')
        self.dDict = {self.emission: X}
        self.inference = ed.MAP(self.paramDict,self.dDict)
        self.inference.run(n_iter = n_iter, n_print=n_print,*kwargs)
        return self
    
    @property
    def means_(self):
        res = self.x_post.mean().eval()
        return res
    @property
    def covariances_(self):
        res = self.x_post.covariance().eval()
        return res
    @property
    def weights_(self):
        res = self.q_pi.eval()
        return res
    def _predict_proba(self,X, N=None, norm = 0):
        ''' self.emission does not work, manually build posterior
'''
        
        N = len(X)
        X_bdc = self.expand_input(X)
        # em.log_prob(X_bdc).eval()
        self.cat = edm.Categorical(probs=self.post.weight, 
                                   sample_shape=N
                                  )
        self.x_post = edm.Mixture(
            cat = self.cat, 
            components = self.postComponents,
            sample_shape=N
        )        
        ll = tf.concat([ comp.log_prob(X)[:,None]
                        for comp in self.postComponents],axis=1) + tf.log( self.post.weight)
        
#         X_loglik = ll = self.x_post.log_prob(X_bdc)
##         ll = tf.reduce_mean(ll,axis=1)  ### over posterior samples
    #     ll = tf.reduce_sum(ll,axis=-1)  ### over dimensions
        logP = ll.eval()   
        return logP
    
# i += 1