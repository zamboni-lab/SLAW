import numpy as np
from sklearn.linear_model import Lasso,LinearRegression
from sklearn.preprocessing import StandardScaler
from scipy.optimize import minimize,Bounds

class maxOptimizer:
  def __init__(self):
    pass

  def get_maximum(self, points, values, removed = -1.0):
    bool = [val!=removed for val in values]
    points = points[bool]
    nvar = points.shape[1]
    values = [values[idx] for idx in range(len(values)) if bool[idx]]
    pindex = values.index(max(values))
    return points[pindex, :], values[pindex] , [True]*nvar

def compute_interactions_terms(points):
  inter_term = []
  for i in range(points.shape[1] - 1):
    for j in range(i + 1, points.shape[1]):
      inter_term.append(points[:, i] * points[:, j])
  return np.stack(inter_term, axis=1)

def compute_square_terms(points):
  inter_term = []
  for i in range(points.shape[1]):
    inter_term.append(np.square(points[:, i]))
  return np.stack(inter_term, axis=1)

def make_interaction_table(points):
  npoints = points.copy()
  return np.concatenate([npoints, compute_interactions_terms(points),
                         compute_square_terms(points)], axis=1)

def compute_polynom(x, supp_args):
  norm, coef, inter = supp_args
  arr = np.array(x)
  arr = arr.reshape(1, len(arr))
  nval = make_interaction_table(arr)
  nval = norm.transform(nval)
  return inter + np.sum(nval * coef)

def fit_surface(points, values):
  mftable = make_interaction_table(points)
  scaler = StandardScaler()
  stable = scaler.fit_transform(mftable)
  lr = LinearRegression(fit_intercept=True,normalize=False)
  vlr = lr.fit(stable, values)
  coef = vlr.coef_
  inter = vlr.intercept_
  return scaler, coef, inter

def get_range(points):
  bmin = np.min(points, axis=0)
  bmax = np.max(points, axis=0)
  return bmin, bmax

def find_approximate_maximum(points, values,impacting=0.1):
  nvalues = [-v for v in values]
  scaler, coef, inter = fit_surface(points, nvalues)
  ###We sum the linear and square coefficient for each variable.
  nvar = points.shape[1]
  ntotal = len(coef)
  lb, ub = get_range(points)
  bounds = Bounds(lb=lb, ub=ub)
  pindex = values.index(max(values))
  x0 = points[pindex, :]
  vmax = minimize(compute_polynom, x0, args=[scaler, coef, inter], method='L-BFGS-B', bounds=bounds)
  ###We calculate the one for which the coefficient is
  coef_linear = abs(coef[0:nvar])
  coef_squared = abs(coef[(ntotal-nvar):ntotal])
  max_contrib = max(max(coef_linear),max(coef_squared))
  threshold = max_contrib*0.1
  valid = [(coef_linear[idx]>threshold or coef_squared[idx]>threshold) for idx in range(len(coef_linear))]
  return vmax.x,-vmax.fun-inter,valid

class rsmOptimizer:
  def __init__(self):
    pass
  def get_maximum(self, points, values, removed=-1.0):
    ##We removed the forbbiden values
    bool = [val!=removed for val in values]
    points = points[bool]
    values = [values[idx] for idx in range(len(values)) if bool[idx]]

    return find_approximate_maximum(points, values)

if __name__=="__main__":
  n_test = 200
  fake_data = np.stack([np.random.uniform(-1,1,n_test),np.random.uniform(-1,1,n_test),np.random.uniform(-1,1,n_test)],axis=0)
  fake_data = fake_data.T
  values = 0.5*fake_data[:,0]**2+1.5*fake_data[:,1]**2+3*fake_data[:,0]+fake_data[:,1]*fake_data[:,0]*2.3+0.05*fake_data[:,2]
  values = values+np.random.normal(len(values))
  ropt = rsmOptimizer()
  ropt.get_maximum(fake_data,values)

  mopt = maxOptimizer()
  mopt.get_maximum(fake_data,values)