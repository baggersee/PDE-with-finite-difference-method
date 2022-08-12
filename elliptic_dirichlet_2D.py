import numpy as np
import matplotlib.pyplot as plt

# Specifications of the problem
f = lambda x,y: 0 # function for the Poisson equation
gx_yBottom = lambda x: np.sin(np.pi*x)  # Boundary condition Dirichlet type in y = y_bottom
gx_yTop = lambda x: np.sin(np.pi*x)     # Boundary condition Dirichlet type in y = y_right
gy_xLeft = lambda y: 0                  # Boundary condition Dirichlet type in x = x_left
gy_xRight = lambda y: 0                 # Boundary condition Dirichlet type in x = x_right

x_left = 0
x_right = 1
y_bottom = 0
y_top = 1

# Discretization data
m = 12 # number of nodes in the x-axis
n = 12 # number of nodes in the y-axis

h = (x_right - x_left)/(m-1) # space for the x-axis discretization
k = (y_top - y_bottom)/(n-1) # space for the y-axis discretization
t = -2*(1/(h*h) + 1/(k*k))

X = np.array([ x_left + t*h for t in range(m)])
Y = np.array([ y_bottom + t*k for t in range(m)])

# Matrix calculation
A = np.diag([1]*m*n)
for i in range(1,m-1):
    for j in range(1,n-1):
        p = i + (j)*m 
        A[p, p] = t
        A[p, i + 1 + (j)*m] = 1/(h*h)
        A[p, i - 1 + (j)*m] = 1/(h*h)
        A[p, i + (j + 1)*m] = 1/(k*k)
        A[p, i + (j - 1)*m] = 1/(k*k)
A_ = np.linalg.inv(A)

b = np.array([gx_yBottom(X[j]) for j in range(m)])
for j in range(1,n-1):
    x = np.array([f(X[i],Y[j]) for i in range(m)])
    x[0] = gy_xLeft(Y[i])
    x[-1] = gy_xRight(Y[i])
    b = np.concatenate((b,x))
b = np.concatenate((b,np.array([gx_yTop(X[j]) for j in range(m)]))) 

# Solving
U = np.matmul(A_,b)
U = U.reshape(m,n,order='F').copy()
U = U.transpose()
X, Y = np.meshgrid(X, Y)

# Plotting 3D

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, U)
plt.grid()
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_zlim(-0.5,1)
plt.show()




