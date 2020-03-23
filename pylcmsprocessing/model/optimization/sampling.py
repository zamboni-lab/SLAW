import math
import numpy as np
import inspect
import concurrent.futures

from model.rsm import bbd
import model.optimization.optimizer as moo



class samplingOptimizer:
    def __init__(self,sampler,optimizer,bounds,fixed_arguments=None):
        if fixed_arguments is None:
            self.fixed_arguments = {}
        else:
            self.fixed_arguments = fixed_arguments
        self.bounds=bounds
        self.sampler=sampler
        self.optimizer=optimizer

    def optimize(self,func,relative_increase = 0.02,max_its=3):

        ###We add some points
        global_best_point = None
        global_best_value = 0
        num_its = 0
        self.sampler.sample(bounds=self.bounds,func=func,fixed_arguments=self.fixed_arguments)
        current_best_point,current_best_value = self.optimizer.get_maximum(self.sampler.get_points(),self.sampler.get_values())
        ###We always test the best point
        while abs(global_best_value-current_best_value)/global_best_value < (1+relative_increase) and num_its < max_its:
            ###We restrain thge constraints using the newly determined best points
            self.bounds = self.bounds.contract_bound(current_best_point, contraction=0.5, extension=0.1,
                                     extreme=0.02, only_positive=True)
            self.sampler.sample(bounds=self.bounds,func=func,fixed_arguments=self.fixed_arguments)
            current_best_point, current_best_value = self.optimizer.get_maximum(self.sampler.get_points(), self.sampler.get_values())
            if current_best_value>global_best_value:
                global_best_point = current_best_point
            num_its += 1

        ##We pick the bset sampled points
        all_values = self.sampler.get_values()
        all_points = self.sampler.get_points()
        pindex = all_values.index(max(all_values))
        return all_points[pindex,:]




class bounds:
    def __init__(self,lb,ub=None):
        if ub is None:
            self.ub = [x[0] for x in lb]
            self.lb = [x[1] for x in lb]
        else:
            self.ub = ub
            self.lb = lb

    def lower_bound(self):
        return self.lb

    def upper_bound(self):
        return self.ub

    def contract_bound(self,best_point, contraction=0.5, extension=0.1,
                                 extreme=0.02, only_positive=True):
        tnlb = [0.0] * len(self.lb)
        tnub = [0.0] * len(self.lb)
        crange = [(cub - clb) / 2 for cub, clb in zip(self.ub, self.lb)]

        for ip in range(len(self.lb)):
            nlb = best_point[ip] - crange[ip] * contraction
            nub = best_point[ip] + crange[ip] * contraction

            if (abs((best_point[ip] - self.lb[ip]) / crange[ip]) - 1) <= extreme:
                nlb = self.lb[ip] - crange[ip] * contraction * extension * 2
            else:
                if nlb < self.lb[ip]:
                    nlb = self.lb[ip]

            if (abs((best_point[ip] - self.ub[ip]) / crange[ip]) - 1) <= extreme:
                nub = self.ub[ip] + crange[ip] * contraction * extension * 2
            else:
                if nub > self.ub[ip]:
                    nub = self.ub[ip]

            if only_positive:
                if nlb < 0: nlb = 0

            tnlb[ip] = nlb
            tnub[ip] = nub
        self.lb = tnlb
        self.ub = tnub


class boundedSampler:
    def __init__(self):
        ###Self. point should be a list
        self.points = []
        self.values = []

    ##This method jsut need to inherit the first one
    def sample(self,bounds,func,fixed_arguments=None):
        points,values=self.sample_points(bounds,func,fixed_arguments)
        self.append_points(points,values)


    def append_points(self,points,values):
        self.points.append(points)
        self.values.append(values)

    def get_points(self):
        return np.concatenate(self.points)

    def get_values(self):
        return [y for x in self.values for y in x]



def wrap_func_dic(x):
  func = x.pop("fun")
  ###We compute the function value
  return func(**x)

def do_uniform_random_sampling(lb,ub,func,num_points=1000,num_cores = None,fixed_arguments=None):
  '''
  :param constraints: a scipy.optimize.Bounds object
  :param func: The function ot be optimized
  :param max_call: The maximum number of call to the function allowed INCLUDING the initial points
  :param initial_points: The number of initial point sused for bounds estimation
  :return: The optimized paramters
  '''
  ###We count the number of arguments
  if num_cores is None:
    num_cores = 1

  if fixed_arguments is None:
    fixed_arguments = {}
  args_name = inspect.getfullargspec(func)[0]
  to_optimize = [name for name in args_name if name not in fixed_arguments]

  ##We read the dimension from the contraints
  ndim = len(lb)
  if len(to_optimize) != ndim:
    raise Exception("Arguments don't match.")

  val_par = [None]*ndim

  ###We just sample every parameters across the different
  for idx, (ilb,iub) in enumerate(zip(lb,ub)):
    par_seq = np.random.uniform(ilb, iub, num_points)
    val_par[idx] = par_seq

  ###We reshape the data frame in a list of dictionnary
  all_dict = [{**dict(zip(to_optimize, cargs)),**fixed_arguments,"fun":func} for cargs in zip(*val_par)]
  list_for_numpy = [cargs for cargs in zip(*val_par)]
  vres = map(wrap_func_dic, all_dict)
  # with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
  #   vres = list(executor.map(wrap_func_dic, all_dict))

  ###We build a tuple of parameters
  np_table = np.array(list_for_numpy)
  return np_table,list(vres)

class uniformBoundedSampler(boundedSampler):
    def __init__(self):
        super().__init__()

    def sample_points(self,limits,func,fixed_arguments=None):
        points,values = do_uniform_random_sampling(limits.lower_bound(), limits.upper_bound(), func, num_points=num_points,
                                                   num_cores=None, fixed_arguments=fixed_arguments)
        return points,values


def do_bdd(lb, ub, func, num_cores=None, fixed_arguments=None):
        '''
        :param constraints: a scipy.optimize.Bounds object
        :param func: The function ot be optimized
        :param max_call: The maximum number of call to the function allowed INCLUDING the initial points
        :param initial_points: The number of initial point sused for bounds estimation
        :return: The optimized paramters
        '''
        ###We count the number of arguments
        if num_cores is None:
            num_cores = 1

        if fixed_arguments is None:
            fixed_arguments = {}
        args_name = inspect.getfullargspec(func)[0]
        to_optimize = [name for name in args_name if name not in fixed_arguments]

        ##We read the dimension from the contraints
        ndim = len(lb)
        print(lb,ub)
        bbd = bbdesign(ndim, center=1)
        for it in range(ndim):
            bbd[:, it] = bbd[:, it] * (ub[it] - lb[it]) / 2 + (ub[it] + lb[it]) / 2
        lpar = bbd.tolist()

        ###We reshape the data frame in a list of dictionnary
        all_dict = [{**dict(zip(to_optimize, cargs)), **fixed_arguments, "fun": func} for cargs in lpar]
        list_for_numpy = [cargs for cargs in zip(*lpar)]
        # with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        #     vres = list(map(wrap_func_dic, all_dict))
        vres = list(map(wrap_func_dic, all_dict))
        np_table = np.array(list_for_numpy)
        return np_table.T,vres

class bbdSampler(boundedSampler):
    def __init__(self):
        super().__init__()

    def sample_points(self,bounds,func,fixed_arguments):
        return do_bdd(bounds.lower_bound(), bounds.upper_bound() , func, num_cores=None, fixed_arguments=fixed_arguments)



class boundedSamplerLIPO(boundedSampler):
    def sample_with_limits(self, func, max_call=50, initial_points=3, fixed_arguments=None):
        vp = moo.LIPO(self.lb, self.ub, func, max_call=max_call, initial_points=initial_points, fixed_arguments=fixed_arguments)



if __name__=="__main__":
  lb=[-1.5,-2.5,0.05]
  ub=[3.5,2.5,10]
  def test(x,y,z,k):
    return -5*(10*(x+y*z)+2*(x+y)**2)+k+3*z**2+0.5*z*x

  ob = bounds(lb,ub)
  bbsampler = bbdSampler()
  rsm_optimizer = rsmOptimizer()
    sO = samplingOptimizer(bbsampler,rsm_optimizer,ob,fixed_arguments={"k":3})
    best_point = sO.optimize(test)

        def optimize(self,func,relative_increase = 0.02,max_its=3):