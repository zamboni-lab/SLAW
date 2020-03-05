import math
###We impote the minimization function
import numpy as np
import inspect
import concurrent.futures


###We consider that

def upper_bound(x,supp):
  points, val_points, k = supp
  '''
  :param x: The point of evaluation of the fuction
  :param points: The prviosuly evalautated points
  :param val_points: The function values of the prviously evaluated point
  :param k: The lipizitchian constant k
  :return:  The valu eof the upper bounds in k
  '''
  ###We strat by substracting the value of X to all elements of points
  diff_points = [p-x for p in points]
  m_diff_points = np.array(diff_points)

  ###We tranform k in a squre matrix
  pow_k = [kk**2 for kk in k]
  pow_k = np.diag(pow_k)

  ###We compute the k derived part
  k_term = np.sqrt(np.apply_along_axis(lambda x, k: np.matmul(np.matmul(x, k), x.T), 1, m_diff_points, k=pow_k))
  fun_values = k_term+val_points

  ###This return the value of the upper bound
  return -np.min(fun_values)

####The upper bound which is piecewise linear can be approximated easly.
def maximize_upper_bound(points,val_points,k,lbi,ubi):
  ndim = points.shape[1]
  max_point = [0]*ndim
  ##We find the best coo
  points = np.asarray(points)
  ndim = points.shape[1]
  npoints = points.shape[0]
  val_points=np.asarray(val_points)
  xorder = np.apply_along_axis(np.argsort,0,points).T
  for d in range(ndim):
    dxorder = xorder[d,:]
    dxorder = dxorder.astype("int_")
    ub = ubi[d]
    lb = lbi[d]
    ##Al the dimesnion values
    xval = points[dxorder,d]
    yval = val_points[dxorder]
    kval = k[d]
    sx = xval[0:(npoints-1)]+xval[1:npoints]
    dy = np.diff(yval)
    middle_x = 0.5*(dy/kval+sx)
    middle_y = yval[0:(npoints-1)]+(middle_x-xval[0:(npoints-1)])*kval

    ###We add the bound if necessary
    ##Lower bound
    if np.abs(xval[0]-lb) > 0.00001:
      middle_x = np.insert(middle_x,0,lb)
      nval = yval[0]+(xval[0]-lb)*kval

      ###We insert the indices at the first positon
      middle_y = np.insert(middle_y,0,nval)
    if (np.abs(xval[npoints-1] - ub) > 0.00001):
      middle_x = np.append(middle_x, ub)
      nval = yval[npoints-1]+(ub-xval[npoints-1])*kval
      middle_y = np.append(middle_y,nval)
    ###We fidn the index of the maximum middle point
    index = np.where(middle_y == np.amax(middle_y))
    ###We return the best candidate
    try:
      max_point[d]=middle_x[index].item()
    except ValueError:
      ##In this case we return None
      return None
  return max_point


def calculateLipschitzConstant(points,val_points):
  ###We calculate the highest slope sbetween neibouring point
  mpoints = np.asarray(points)
  ndim = mpoints.shape[1]
  xorder = np.apply_along_axis(np.argsort,0,mpoints).T
  val_points=np.asarray(val_points)
  k = [None]*ndim
  for d in range(ndim):
    dxorder = xorder[d,:]
    dxorder = dxorder.astype("int_")
    dy = np.diff(val_points[dxorder])
    dx = np.diff(mpoints[dxorder,d])
    slopes = np.abs(dy/dx)
    ###We take the maximu slopes as the possible slopes
    k[d]=np.max(slopes)
  return k

def LIPO(lb,ub,func,max_call=50,initial_points = 4,fixed_arguments=None):
  '''
  :param constraints: a scipy.optimize.Bounds object
  :param func: The function ot be optimized
  :param max_call: The maximum number of call to the function allowed INCLUDING the initial points
  :param initial_points: The number of initial point sused for bounds estimation
  :return: The optimized paramters
  '''
  ###We count the number of arguments
  if fixed_arguments is None:
    fixed_arguments = {}
  args_name = inspect.getfullargspec(func)[0]

  to_optimize = [name for name in args_name if name not in fixed_arguments]
  # print(to_optimize)

  ##We read the dimension from the contraints
  ndim = len(lb)
  if len(to_optimize) != ndim:
    raise Exception("Arguments don't match.")

  counter = 0
  vpar = [None]*ndim
  ###We sample each paraneters across the different dimension
  for idx, (ilb,iub) in enumerate(zip(lb,ub)):
    par_seq = np.linspace(ilb, iub, 2*initial_points+1)
    par_seq=par_seq[range(1,2*initial_points+1,2)]
    vpar[idx] = par_seq

  ###We initialize the points matrix
  points = np.array(vpar).T
  counter += initial_points
  points = np.append(points, [[0.0]*ndim]*(max_call-counter), axis=0)

  val_points=[0.0]*initial_points

  def wrap_func(x,min_names,fixed):
    dict_min = dict(zip(min_names,x))
    dict_call = {**dict_min,**fixed}
    ###We compute the function value
    return func(**dict_call)

  for ip in range(len(val_points)):
    # dict_call = dict(zip(nn,points[ip,:]))
    ##We add the fixed arugments.
    val_points[ip] = wrap_func(points[ip,:],to_optimize,fixed_arguments)
    ###We compute the function value

  ###Computing the intiali Lipshictz constant function
  k =calculateLipschitzConstant(points[range(counter),:],val_points)

  while counter < max_call:
    ###we start by finding a first minimization point
    x0 = [0.0] * ndim
    ###We sample the inital k values
    for ik in range(ndim):
      ##The intial guess of the data.
      x0[ik]= np.random.uniform(lb[ik],ub[ik],1)[0]
    new_point = maximize_upper_bound(points[range(counter),:], val_points, k, lb,ub)
    ###if the returned point is None it s that we already found an upper boud.
    if new_point is None:
      break
    ###If the point is already present
    val_new_point = wrap_func(new_point,to_optimize,fixed_arguments)
    val_points.append(val_new_point)
    points[counter,:] = new_point
    ##We recompute the Lipschitz constant evnetually
    counter = counter+1
    k = calculateLipschitzConstant(points[range(counter),:], val_points)
  ###We just have oto return the best values now
  mpos = val_points.index(max(val_points))
  return(points[mpos,:])

###Stupid random grid search
def random_grid(lb,ub,func,max_call=1000,num_points=1000,fixed_arguments=None):
  '''
  :param constraints: a scipy.optimize.Bounds object
  :param func: The function ot be optimized
  :param max_call: The maximum number of call to the function allowed INCLUDING the initial points
  :param initial_points: The number of initial point sused for bounds estimation
  :return: The optimized paramters
  '''
  ###We count the number of arguments
  if fixed_arguments is None:
    fixed_arguments = {}
  args_name = inspect.getfullargspec(func)[0]

  to_optimize = [name for name in args_name if name not in fixed_arguments]
  # print(to_optimize)

  ##We read the dimension from the contraints
  ndim = len(lb)
  if len(to_optimize) != ndim:
    raise Exception("Arguments don't match.")

  counter = 0
  vpar = [None]*ndim
  ##We first inestigate hos may point

  ###We just sample every parameters across the different
  for idx, (ilb,iub) in enumerate(zip(lb,ub)):
    par_seq = np.random.uniform(ilb, iub, num_points)
    vpar[idx] = par_seq

  ###We initialize the points matrix
  points = np.array(vpar).T


  def wrap_func(x,min_names,fixed):
    dict_min = dict(zip(min_names,x))
    dict_call = {**dict_min,**fixed}
    ###We compute the function value
    return func(**dict_call)

  with concurrent.futures.ProcessPoolExecutor(self.max_jobs) as executor:
    executor.map(run_cl, largs)


  while counter < max_call:
    ###we start by finding a first minimization point
    x0 = [0.0] * ndim
    ###We sample the inital k values
    for ik in range(ndim):
      ##The intial guess of the data.
      x0[ik]= np.random.uniform(lb[ik],ub[ik],1)[0]
    new_point = maximize_upper_bound(points[range(counter),:], val_points, k, lb,ub)
    ###if the returned point is None it s that we already found an upper boud.
    if new_point is None:
      break
    ###If the point is already present

    ###This is done in parallel for every paramters
    val_new_point = wrap_func(new_point,to_optimize,fixed_arguments)
    val_points.append(val_new_point)
    points[counter,:] = new_point
    ##We recompute the Lipschitz constant evnetually
    counter = counter+1
    k = calculateLipschitzConstant(points[range(counter),:], val_points)
  ###We just have oto return the best values now
  mpos = val_points.index(max(val_points))
  return(points[mpos,:])




if __name__=="__main__":
  lb=[-1.5,-2.5]
  ub=[3.5,2.5]
  def test(x,y,k):
    return -5*(np.sin(10*(x+y))+2*(x+y)**2)+k
  voptim = LIPO(lb,ub,test,10,4,fixed_arguments={"k":2})

  lb2 = [-5]
  ub2 = [5]
  def test2(x):
    return -(5*np.sin(10*x)+2*(x**2))
  voptim2 = LIPO(lb2,ub2,test2,50,4)


# # Make data.
# X = np.arange(-2.5, 2.5, 0.25)
# Y = np.arange(-2.5, 2.5, 0.25)
# X, Y = np.meshgrid(X, Y)
# Z = test(X,Y,k=2)
#
#
# fig = plt.figure()
# ax = Axes3D(fig)
# surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# ax.set_zlim(-10, 150)
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# fig.colorbar(surf, shrink=0.5, aspect=5)
#
# plt.show()