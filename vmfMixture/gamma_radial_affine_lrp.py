# import models; reload(models);from models import * 

# import tensorflow_probability.python.edward2 as tfed
import tensorflow_probability.python.distributions as tfdist
import pymisca.tensorflow_extra as pytf
ed = edm = pytf.ed
tf = pytf.tf;

from models import BaseModel

import pymisca.util as pyutil
np = pyutil.np

class GammaRadialAffineLRP_VIMAP(BaseModel):
    
    bkd = tfdist
    
    def __init__(self,D=None,K=20,debug=False,
                 *args,**kwargs):
        super(GammaRadialAffineLRP_VIMAP,self).__init__(*args,**kwargs)
        self.K = K
        self.D = D
        self.initialised = False
        self.sess = None
        self.feed_dict = None
        self.debug = debug
        
        if D is not None:
            self.init_model(D=D)  
                        
    em_key =[
                'scale_diag',
                'scale_perturb_factor',
#                 'loc',
                'concentration',
                'rate',
            ]
    mix_key = [
            'weight',
        ]
    
    def random(self,size):
        '''
        random initialiser
'''
        
        return np.random.normal(0.,1.,size=size).astype(np.float32)
    
    def emDist(self,*args,**kwargs):
        #### Construct emission from purely tf.Variable()
        #### Avoid tfdist.Distritbution()
        
        D = self.D
        
        dft = {'concentration':
#                    tfdist.Deterministic
                   (
                       tf.ones([])* self.D/2
                   ) ,
                'rate':
#                tfdist.Deterministic
               (tf.ones([]))
#                    tf.constant(1.)
              }
        kw = { k: kwargs.get( k, v) for k,v in dft.items()}
        
#         kw = {k:v for k,v in kwargs.items() if k in keys}
        dist_rsq = self.bkd.Gamma(**kw)
        dist_xyz = pytf.AsRadial(
            distribution = dist_rsq,
            D = self.D,
            debug = self.debug,)
        
        mdl = dist_xyz

        dft = {
            'loc':
#             tfdist.Deterministic
            (tf.zeros([D,])),
            'scale_diag': 
#             tfdist.Deterministic
            (tf.ones([D,])),
#             'scale_perturb_factor':edm.Deterministic(tf.ones([D,1])),
        }
        kw = {k:kwargs.get( k, v) for k,v in dft.items()}
        
        dist_aff = pytf.AffineTransformDiagPlusLowRank(
#         dist_aff = pytf.AffineTransformDiag(
            distribution=dist_xyz,
            **kw)

        mdl = dist_aff
        return mdl
        
    def make_prior(self):
        D = self.D
        K = self.K
        alpha = self.D/2.
        diriAlpha = 0.1
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
            prior.loc = edm.Normal(tf.zeros(D), tf.ones(D), sample_shape=K)            
#             prior.loc =  edm.Uniform(*uspan,sample_shape=(K,D))
            prior.scale_diag =  edm.Uniform(*[0.001,10.],sample_shape=(K,D))
            prior.scale_perturb_factor =  edm.Uniform(*uspan,sample_shape=(K,D,1))
#             prior.scale_perturb_factor = edm.Normal(tf.zeros([1]), tf.ones([1]), 
#                                                     sample_shape=(K,D,))            
#           
            prior.concentration =  edm.Uniform(*[0.01,10.],sample_shape=(K,))
#             prior.concentration =  edm.Uniform(*uspan,sample_shape=(K,))
            prior.rate =  edm.Uniform(*[0.01,10.],sample_shape=(K,))
            
            prior.weight = pi = edm.Dirichlet( float(diriAlpha)/K * tf.ones(K) )            
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
            post.weight = edm.PointMass(
                tf.nn.softmax(
                    tf.get_variable("q_pi",[K]),
#                     tf.Variable(name="q_pi",initial_value = self.random([K]) ),
                )
            )
            
            post.loc = edm.PointMass(
                tf.get_variable("q_mu", [K,D])
            )
        
            
            post.scale_diag  = edm.PointMass(
                tf.nn.softplus(
                    tf.get_variable('q_scale_diag',shape=[K,D])
                              ),
            )
            
            post.scale_perturb_factor = edm.PointMass(
                (
                    tf.get_variable("q_scale_perturb_factor", [K,D,1])
                )
            )
            
            post.concentration  = edm.PointMass(
                tf.nn.softplus(
#                     tf.Variable(name="concentration",initial_value = self.random([K]) ),
                    tf.get_variable('concentration',shape=[K,])
                              ),
            )
            
            post.rate  = edm.PointMass(
                tf.nn.softplus(
#                     tf.Variable(name="rate",initial_value = self.random([K]) ),
                    tf.get_variable('rate',shape=[K,])
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

        self.param_key = (self.em_key + 
                          self.mix_key)

        self.paramDict = {getattr(prior,name,None):
                          getattr(post,name,None) 
                          for name in self.param_key}
        
        self.paramDict = {k:v 
                          for k,v in self.paramDict.items() 
                          if k is not None and v is not None}
        
        
        ### Prior components
        cDicts = [
            {key: v[k] 
             for key,v in prior.__dict__.items() 
             if key in self.em_key} 
            for k in range(K)]
        self.prior.components = [self.emDist(**d) for d in cDicts]

        
        ### Posterior generative
#         edm.Mixture
        cDicts = [
            {key: v[k] 
             for key,v in post.__dict__.items() 
             if key in self.em_key} 
            for k in range(K)]
        self.post.components = [self.emDist(**d) for d in cDicts]
        

    
        self.initialised = True; return self
    
    def build_likelihood(self,X,env=None):
        if env is None:
            env = self.prior
        elif isinstance(env,str):
            env = getattr(self,env)
        K = self.K
        N = len(X)
#         env,cat = bkd
        env.cat = self.bkd.Categorical(probs=env.weight, 
#                                    sample_shape=N
                                  )
        env.emission = self.bkd.Mixture(
            cat=env.cat, 
            components=env.components,
#             sample_shape=N,
        )        
        return env.emission
    
    def _fit(self,X, n_iter=1000, n_print=100, env=None,**kwargs):
        
        K = self.K
        N = len(X)
        
        emission = self.build_likelihood(X,env=env)

        print ('hiModel')
        self.dDict = {emission: X}
        self.inference = ed.MAP( self.paramDict, self.dDict)
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
        res = self.post.weight.eval()
        return res
    def _predict_proba(self,X, N=None, norm = 0):
        ''' self.emission does not work, manually build posterior
'''
        assert self.sess is not None,"\"self.sess\" is None, please fit the model first with a tf.Session()"
    
        N = len(X)
        X_bdc = self.expand_input(X)
        

        ll = tf.concat([ comp.log_prob(X)[:,None]
                        for comp in self.post.components],axis=1) 
        ll = ll + tf.log( self.post.weight )
        
##         ll = tf.reduce_mean(ll,axis=1)  ### over posterior samples
    #     ll = tf.reduce_sum(ll,axis=-1)  ### over dimensions
        logP = ll.eval(session=self.sess,
                      feed_dict=self.feed_dict)   
        return logP
    
#     import tensorflow as tf
    def freeVarDict(self, env):
        if isinstance(env,dict):
            idict = env
        else:
            if isinstance(env,str):
                env = getattr(self,env)
            else:
                raise Exception('"env" not recognised:%s'%env)
            idict = env.__dict__
        odict = {k:
#                  x.value() ### for edward.PointMass()
                 x
                 for k,x in 
                       idict.iteritems() 
                       if not isinstance(x,list) and k in self.param_key}            
        return odict
        
            
    def _fit_MAP(self,x_obs):
        mdl = self.build_likelihood(X=x_obs,env='post')
        x_place = tf.placeholder(dtype='float32')
        
        ### prior likelihood
        lp_param = [
            tf.reduce_sum(k.log_prob(v.value()))  ### ed.RandomVariable.value()
            if v.__class__.__name__ == 'Uniform' 
            else 0.
            for k,v in self.paramDict.items()]
#         print (tf.reduce_sum(lp_param))
        ### data likelihood
        lp_data = tf.reduce_mean(
            mdl.log_prob(x_place)
        )
        lp = tf.reduce_sum(
            map(tf.reduce_sum,[lp_param,lp_data])
        )
        loss = -lp
#         loss = 0.
#         loss += tf.reduce_sum(-lp_param) + tf.reduce_sum()

        self.feed_dict = {x_place:x_obs}
#         fitted_vars_dict = {k:x.value() for k,x in 
#                        self.post.__dict__.iteritems() 
#                        if not isinstance(x,list) and k in self.param_key}
        # fitted_vars
        reload(pytf)
        LEARNING_RATE = 0.1
        optimizer = tf.train.AdamOptimizer(learning_rate=LEARNING_RATE)
        sess, last_vars, hist_loss, opt = pytf.op_minimise(
            loss,
            self.freeVarDict('post').values(),
            optimizer,
            feed_dict = self.feed_dict
        )
        self.sess = sess
        return mdl,last_vars
    