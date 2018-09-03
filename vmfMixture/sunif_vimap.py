# import models; reload(models);from models import * 
import pymisca.tensorflow_extra as pytf
edm = pytf.edm; tf = pytf.tf;
import edward as ed
from models import BaseModel
import pymisca.util as pyutil

class SunifLRP_VIMAP(BaseModel):
#     self.emDist = edm.MultivariateNormalDiagPlusLowRank    
#     emDist = edm.MultivariateNormalDiagPlusLowRank    
#     emDist = edm.sphereUniformLRP
    def __init__(self,D=None,K=20,
                 *args,**kwargs):
        super(SunifLRP_VIMAP,self).__init__(*args,**kwargs)
        self.K = K
        self.D = D
        self.initialised = False
        if D is not None:
            self.init_model(D=D)   

    em_key =[
                'scale_diag',
                'scale_perturb_factor',
                'concentration',
                'rate',
            ]
    def emDist(self,*args,**kwargs):
        kwargs = {k:v for k,v in kwargs.items() if k in self.em_key}
        return edm.sphereUniformLRP(*args,**kwargs)
        
        
    def make_prior(self):
        D = self.D
        K = self.K
        alpha = self.D/2.
        name = self.name
        self.prior = prior = pyutil.util_obj()
        try:
            tf.get_variable(name+'/prior',[1])
            reuse = None
        except:
            reuse = True
        print ('reuse',reuse)

        with tf.variable_scope(name, reuse=reuse):
            
            uspan = [-1E5,1E5]
            ##### Prior
            prior.mu = edm.Normal(tf.zeros(D), tf.ones(D), sample_shape=K)            
            prior.scale_diag =  edm.Uniform(*uspan,sample_shape=(K,D))
            prior.scale_perturb_factor =  edm.Uniform(*uspan,sample_shape=(K,D,1))   
            prior.concentration =  edm.Uniform(*uspan,sample_shape=(K,1))
            prior.rate =  edm.Uniform(*uspan,sample_shape=(K,1))
            prior.weight = pi = edm.Dirichlet( float(alpha)/K * tf.ones(K) )
        return prior
    
    def make_post(self):
        D = self.D
        K = self.K
        alpha = self.D/2.
        name = self.name
        self.post = post = pyutil.util_obj()
        try:
            tf.get_variable(name+'/post',[1])
            reuse = None
        except:
            reuse = True
        print ('reuse',reuse)

        with tf.variable_scope(name, reuse=reuse):
            
            uspan = [-1E5,1E5]
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
        return post
    
    def init_model(self,D=None,K = None,alpha = 1.0):
        self.D = D = self.D if D is None else D
        assert D is not None
        self.K = K = self.K if K is None else K
        assert K is not None
        
        prior = self.make_prior()
        post = self.make_post()
        
        ##### Dictonary for constructing self.emDist(**self.param)
#         self.
        self.mix_key = [
            'weight',
        ]
        self.param_key = (self.em_key + 
                          self.mix_key)


        self.paramDict = {getattr(prior,name):
                          getattr(post,name) 
                          for name in self.param_key}
        
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
#             sample_shape=N
        )        
        ll = tf.concat([ comp.log_prob(X)[:,None]
                        for comp in self.postComponents],axis=1) + tf.log( self.post.weight)
        
#         X_loglik = ll = self.x_post.log_prob(X_bdc)
##         ll = tf.reduce_mean(ll,axis=1)  ### over posterior samples
    #     ll = tf.reduce_sum(ll,axis=-1)  ### over dimensions
        logP = ll.eval()   
        return logP
    