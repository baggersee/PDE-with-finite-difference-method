import numpy as np
import matplotlib.pyplot as plt

# Specifications of the problem
fx_tStart = lambda x: np.sin(2*np.pi*x)*np.sin(2*np.pi*x)  # Initial condition at t = t_start
gt_xLeft = lambda t: 0.   # Boundary condition Dirichlet type in x = x_left
gt_xRight = lambda t: 0.  # Boundary condition Dirichlet type in x = x_right

x_left = 0
x_right = 1
t_start = 0
t_finish = 1

# Discretization data
m = 3 # number of nodes in the x-axis
n = 3 # number of nodes in the t-axis --> number of time instants

h = (x_right - x_left)/(m-1) # space for the x-axis discretization
k = (t_finish - t_start)/(n-1) # space for the y-axis discretization
D = 1 # diffusion coeficient
sigma = D*k/(pow(h,2))

X = np.array([ x_left + i*h for i in range(m)])
T = np.array([ t_start + i*k for i in range(n)])

# matrix construction
A = np.diag([1-2*sigma]*(m-2))
for i in range((m-2)-1):
    A[i,i+1] = sigma
    A[i+1,i] = sigma

U = np.array([fx_tStart(X[i]) for i in range (1,m-1)])
# aqui usar recursividad!!!
for t in range(n):
    U = np.matmul(A,U)
    U[0] = sigma*gt_xLeft(T[t])
    U[-1] = sigma*gt_xRight(T[t])

U = np.concatenate((np.array([gt_xLeft(T[-1])]),U))
U = np.concatenate((np.array([gt_xRight(T[-1])]),U))
