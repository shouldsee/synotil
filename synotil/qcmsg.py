
def msgGMM(model = None, train_data = None,name='test',**kwargs):
    mdl = model
    s = '''
Name:{name}
Converged:{cvg}
min_logL: {logL}
(lower-bound of) MEAN logL :{mlogL}'''.format(
        name=name,
         cvg=mdl.converged_,
         logL = mdl.lower_bound_,
         mlogL=mdl.lower_bound_/len(train_data) if not train_data is None else mdl.lower_bound_
                                )
    return s