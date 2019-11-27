# import models; reload(models);from models import * 

# import tensorflow_probability.python.edward2 as tfed
import tensorflow_probability.python.distributions as tfdist
import pymisca.tensorflow_extra as pytf
ed = edm = pytf.ed
tf = pytf.tf;

from models import BaseModel

import pymisca.util as pyutil
np = pyutil.np

# from pymisca
# VonMisesFisher = 
def GammaRadialTheta( self, 
                      gamma_concentration = None,
                      gamma_rate = None,
                      vm_concentration = None,
                      vm_direction = None,
                      D = None
                    ):
    D = self.D if D is None else D
    if gamma_concentration is None:
        gamma_concentration = tf.ones([]) * self.D/2.
#     mean_direction = tf.convert_to_tensor( (mean_direction),dtype='float32')

    dist_rsqTheta = mdl=  pytf.JointDist(
        [
            tfdist.Gamma(gamma_concentration, 
                         gamma_rate,
                         name='r_squared'),
            pytf.VonMisesFisher(
                concentration = vm_concentration,
                mean_direction = vm_direction
                                 )
        ])

    #### Use a bijector to calculate P(x) from P(r^2)
    dist_xyz = mdl = pytf.AsRadialTheta(
        distribution=mdl,
        D=D)

    ### Allow an affine transformation y = M x + x_0
#     dist_aff = mdl = pytf.AffineTransformDiag(distribution=dist_xyz,
#                                         scale_diag=sigma,)
    # gaussian_dist = tf.contrib.distributions.Normal(loc=mu, scale=sigma)
    # emission  = dist_rsq
    
    return mdl

class GammaRadialTheta_VIMAP(BaseModel):
    
    bkd    = tfdist
    emDist = GammaRadialTheta
    
    def __init__(self,D=None,K=20,
                 debug=False,
                 *args,**kwargs):
        super(
            GammaRadialTheta_VIMAP,
            self).__init__(*args,**kwargs)
        self.K = K
        self.D = D
        self.initialised = False
        self.sess = None
        self.feed_dict = None
        self.debug = debug
        
        if D is not None:
            self.init_model(D=D)  
                        
    em_key =[
        'gamma_concentration',
        'gamma_rate',
        'vm_concentration',
        'vm_direction',
        ]
    mix_key = [
            'weight',
        ]
    
    def random(self,size):
        '''
        random initialiser
'''
        
        return np.random.normal(0.,1.,size=size).astype(np.float32)
    


        
    def make_prior(self):
        D = self.D
        K = self.K
        alpha = self.D/2.
#         diriAlpha = self.K /10.
#         diriAlpha = 1.

        diriAlpha = 0.001
#         diriAlpha = 0.00001
#         diriAlpha = 0.0000000000000000000000000000000001        
#         diriAlpha = 10.

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
#             prior.gamma_concentration = edm.Normal(tf.zeros(D), tf.ones(D), sample_shape=K)            
#             prior.loc =  edm.Uniform(*uspan,sample_shape=(K,D))
            prior.gamma_concentration =  edm.Uniform(*[0.001,1000.],sample_shape=(K,))
            prior.gamma_rate =  edm.Uniform(*[0.001,100000.],sample_shape=(K,))
            prior.vm_concentration =  edm.Uniform(*[0.001,100000.],sample_shape=(K,))
            prior.vm_direction =  edm.Uniform(*[0.001,100000.],sample_shape=(K,D))
#             prior.vm_direction = edm.Normal(tf.zeros(D), tf.ones(D), sample_shape=K)      
 
#             prior.weight =  edm.Uniform(*[0.001,100000.],sample_shape=(K,))
            prior.weight = pi = edm.Dirichlet( float(diriAlpha)/K * tf.ones(K) )            

#             prior.cat = edm.Categorical(weight = post.weight)
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
            i = -1
            i += 1
#             post.weight = edm.PointMass(
#                 tf.nn.softmax(
#                     tf.get_variable(str(i), shape=[K]), 
#                     name = 'weight',
# #                     tf.Variable(name="q_pi",initial_value = self.random([K]) ),
#                 )
#             )
            
            post.weight = edm.PointMass(
                tf.square(
                    tf.nn.l2_normalize(
                        tf.get_variable(str(i), shape=[K]), 
                    ),
                    name = 'weight',
#                     tf.Variable(name="q_pi",initial_value = self.random([K]) ),
                )
            )     


#             post.cat = edm.PointMass(
#                 tf.nn.softmax(
#                     tf.get_variable("cat",[K]),
# #                     tf.Variable(name="q_pi",initial_value = self.random([K]) ),
#                 )
#             )            
        
            i += 1            
            post.gamma_concentration  = edm.PointMass(
                tf.nn.softplus(
#                     tf.Variable(name="concentration",initial_value = self.random([K]) ),
                    tf.get_variable(str(i),shape=[K,]),
                    name = 'gamma_concentration',
                              ),
            )

            i += 1
            post.gamma_rate  = edm.PointMass(
                tf.nn.softplus(
#                     tf.Variable(name="concentration",initial_value = self.random([K]) ),
                    tf.get_variable(str(i),shape=[K,]),
                    name = 'gamma_rate',
                              ),
            )
            
            i += 1
            post.vm_concentration  = edm.PointMass(
                0.0 + tf.nn.softplus(
                    10. -  tf.nn.softplus(
    #                     tf.Variable(name="concentration",initial_value = self.random([K]) ),
                        tf.get_variable(str(i),shape=[K,]),
                        
                                  ),
                name = 'vm_concentration',)
            )            
            
#             post.vm_concentration  = edm.PointMass(
#                 10. *  tf.nn.sigmoid(
# #                     tf.Variable(name="concentration",initial_value = self.random([K]) ),
#                     tf.get_variable('vm_concentration',shape=[K,])
#                               ),
#             )            
            i += 1
            post.vm_direction = edm.PointMass(
                tf.nn.l2_normalize(
                    tf.get_variable(str(i), [K,D]),
                    axis = -1,
                    name = "vm_direction",
                ),
            )
            
#             post.rate  = edm.PointMass(
#                 tf.nn.softplus(
# #                     tf.Variable(name="rate",initial_value = self.random([K]) ),
#                     tf.get_variable('rate',shape=[K,])
#                               ),
#             )
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

        ### build 
        env.cat = self.bkd.Categorical(
            probs = env.weight,
#             sample_shape=N
        )
#         env.cat = edm.PointMass(
#             tf.nn.softmax(
#                 tf.get_variable("cat",[N,K]),
# #                     tf.Variable(name="q_pi",initial_value = self.random([K]) ),
#             )
#         )            


        env.emission = self.bkd.Mixture(
            cat = env.cat, 
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
        
            
    def _fit_MAP(self,x_obs,**kwargs):
        mdl = self.build_likelihood(X=x_obs,env='post')
        x_place = tf.placeholder(dtype='float32')
        
        ### prior likelihood
        self.lp_param = lp_param = [
            tf.reduce_sum(
                k.distribution.log_prob(v)  ### ed.RandomVariable.value()
            ) 
            if k.distribution.__class__.__name__ != 'Uniform' 
            else 0.
            for k,v in self.paramDict.items()]
#         print (tf.reduce_sum(lp_param))
        ### data likelihood
        self.lp_data = lp_data = tf.reduce_sum(
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
#         reload(pytf)
#         LEARNING_RATE = 0.1
#         LEARNING_RATE = 1.
#         optimizer = tf.train.AdamOptimizer(learning_rate=LEARNING_RATE)

        LEARNING_RATE = 1.
#         LEARNING_RATE = 0.1
#         LEARNING_RATE = 10.
        optimizer = tf.train.AdadeltaOptimizer(learning_rate=LEARNING_RATE)

#         LEARNING_RATE = 0.0001
#         optimizer = tf.train.GradientDescentOptimizer(learning_rate=LEARNING_RATE)
        
        sess, last_vars, hist_loss, opt = pytf.op_minimise(
            loss,
            self.freeVarDict('post').values(),
            optimizer,
            feed_dict = self.feed_dict,
            **kwargs
            
        )
        self.sess = sess
        return mdl,(last_vars, hist_loss, opt)
    