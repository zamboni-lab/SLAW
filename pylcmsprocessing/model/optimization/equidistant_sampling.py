##Sampling theory wth sumlatd annealing
import numpy as np
from scipy.spatial import distance
from scipy.optimize import basinhopping,Bounds

class boundedGrid:
    def __init__(self,bounds=None,ub=None,lb=None,weight=None,limit = True):
        if bounds is not None:
            self.lb = bounds.lb
            self.ub = bounds.ub
        else:
            self.lb = lb
            self.ub = ub
        self.weight = [1]*len(self.ub)
        if weight is not None:
            self.weight=weight
        self.limit = limit

    def sample(self,npoints):
        if npoints>50:
            ####Given the fact that the exponential point
            print("More than 50 points picked with an optimal grid sampling only 50 points will be picked")
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
        def summed_distance(unrolled_points,ndim = 3):
            points = unrolled_points.reshape((ndim, int(unrolled_points.shape[0] / ndim)))
            cdist = distance.cdist(points.T, points.T)
            return -np.sum(np.sort(cdist, axis=1)[:,1])
        unrolled_table = final_table.reshape((1,np.prod(final_table.shape)))
        nlb = np.zeros(npoints*ndim)
        nub = np.repeat(self.weight,npoints)
        bounds = Bounds(lb=nlb, ub=nub)
        minimizer_kwargs = {"method":'L-BFGS-B', "bounds":bounds,"args":ndim}
        eval = basinhopping(summed_distance,x0 = unrolled_table,niter=int(200/factor),stepsize=0.01,minimizer_kwargs=minimizer_kwargs)
        ###We now reshape the array
        res = eval.x.reshape(ndim,npoints)
        for idx in range(len(self.ub)):
            res[idx,:]=(res[idx,:]/self.weight[idx])*(self.ub[idx]-self.lb[idx])+self.lb[idx]
        return res

##Test a very stupid idea with rtree to sample the datasts.
if __name__=="__main__":
    lb = np.array([0.005,0.01,10,1000,500,2150])
    ub = 2*lb
    bg = boundedGrid(lb=lb,ub=ub,weight=[1,1,1,1,1,1000])
    tt = bg.sample(9)
    pyplot.scatter(tt[1,:],tt[5,:])
