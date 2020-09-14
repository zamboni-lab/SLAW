##Sampling theory wth sumlatd annealing
import numpy as np
import pickle
import logging
from scipy.spatial import distance,KDTree
from scipy.stats import hmean
from scipy.optimize import basinhopping,Bounds,shgo

import pylcmsprocessing.common.references as cr

class boundedGrid:
    def __init__(self,bounds=None,ub=None,lb=None,weight=None,limit = True):
        if bounds is not None:
            self.lb = bounds.lb
            self.ub = bounds.ub
        else:
            self.lb = lb
            self.ub = ub
        self.weight = [1]*len(self.ub)
        # if weight is not None:
        #     self.weight=weight
        self.limit = limit

    def sample(self,npoints,mode="kd"):
        ##We first check if the number of point exits
        with open(cr.DATA["OPTIMIZATION"]["BALANCED_POINTS"], "rb") as f:
            stored_data = pickle.load(f)
            try:
                logging.debug("Retrieved grid ndim:"+str(len(self.ub))+" npoints:"+str(npoints))
                saved_data = stored_data[(len(self.ub),npoints)][0]
                for idx in range(len(self.ub)):
                    saved_data[idx, :] = saved_data[idx, :] * (self.ub[idx] - self.lb[idx]) + self.lb[idx]
                return saved_data
            except KeyError:
                pass

        if npoints>60:
            ####Given the fact that the exponential point
            logging.info("More than 50 points picked with an optimal grid sampling only 50 points will be picked")
            npoints = 50
        coords = [0]*len(self.ub)
        for idx in range(len(self.ub)):
            coords[idx] = np.random.uniform(0,1,npoints)*self.weight[idx]
        factor = 1
        if self.limit:
            factor = int(npoints/10)
        if factor==0:
            factor=1
        final_table = np.stack(coords)
        ndim = len(self.lb)

        def summed_distance_kdtree(unrolled_points,ndim = 3):
            points = unrolled_points.reshape((ndim, int(unrolled_points.shape[0] / ndim)))
            kd = KDTree(points.T)
            dist, ind = kd.query(points.T, k=2, p=1)
            return -np.mean(dist[:,1])

        def summed_distance(unrolled_points,ndim = 3):
            points = unrolled_points.reshape((ndim, int(unrolled_points.shape[0] / ndim)))
            cdist = distance.cdist(points.T, points.T, p=2.)
            return -np.mean(np.sort(cdist, axis=1)[:,1])

        vfun = summed_distance_kdtree
        if mode != "kd":
            vfun = summed_distance

        unrolled_table = final_table.reshape((1,np.prod(final_table.shape)))
        nlb = np.zeros(npoints*ndim)
        nub = np.repeat(self.weight,npoints)
        bounds = Bounds(lb=nlb, ub=nub)
        minimizer_kwargs = {"method":'L-BFGS-B', "bounds":bounds,"args":ndim}
        niter= int(200/factor)
        niter = 10
        ###Let try it with the shgo algorithm
        bound_shgo = [(lb,ub) for lb,ub in zip(nlb,nub)]

        # eval = shgo(vfun,bound_shgo,n=100,iters=2,args=(ndim,))
        eval = basinhopping(vfun,x0 = unrolled_table,niter=10,stepsize=0.1,minimizer_kwargs=minimizer_kwargs)
        ###We now reshape the array
        res = eval.x.reshape(ndim,npoints)
        for idx in range(len(self.ub)):
            res[idx,:]=(res[idx,:]/self.weight[idx])*(self.ub[idx]-self.lb[idx])+self.lb[idx]
        return res

##Test a very stupid idea with rtree to sample the datasts.
if __name__=="__main__":
    lb = np.array([0.005,0.01,10])
    ub = 2*lb
    bg = boundedGrid(lb=lb,ub=ub,weight=[1,1,1])

    import time
    t1 = time.time()
    tt = bg.sample(15)
    t2 = time.time()

    import time
    t1 = time.time()
    tt = bg.sample(25,"euclidian")
    t2 = time.time()
    from matplotlib import pyplot
    ax = pyplot.axes(projection='3d')
    ax.scatter3D(tt[0,:], tt[1,:], tt[2,:])