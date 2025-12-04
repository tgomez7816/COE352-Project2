import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt

#domain and n nodes
N = 11  #n (can be changed)
x_start, x_end = 0.0, 1.0  #domain in x
t_start, t_end = 0.0, 1.0  #time interval
dx = (x_end - x_start) / (N - 1)  #element size
x_nodes = np.linspace(x_start, x_end, N)  #node position

#time step parameters
dt = 0.01  #step size
num_steps = int((t_end - t_start) / dt)

#define the r.h.s f(x, t)
def f(x, t):
    return (np.pi**2 - 1) * np.exp(-t) * np.sin(np.pi * x)

#initial condition u(x, 0) = sin(pi x)
def initial_condition(x):
    return np.sin(np.pi * x)

#dirichlet boundary conditions
def boundary_left(t):
    return 0.0

def boundary_right(t):
    return 0.0

#2-point Gaussian quadrature points + weights [-1, 1]
gauss_pts = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
gauss_wts = np.array([1.0, 1.0])

#basis functions and their derivatives for linear elements
def lagrange_basis(xi):
    # xi: local coordinate in [-1, 1]
    N = np.array([(1 - xi) / 2, (1 + xi) / 2])
    dN_dxi = np.array([-0.5, 0.5])
    return N, dN_dxi

#global mass and stiffness matrices
M = np.zeros((N, N))  # Mass 
K = np.zeros((N, N))  # Stiffness 

for elem in range(N - 1):
    #node index for element
    n1, n2 = elem, elem + 1
    #n positions
    x1, x2 = x_nodes[n1], x_nodes[n2]
    #jacobian for mapping [-1, 1] to [x1, x2]
    J = (x2 - x1) / 2

    #local mass and stiffness
    M_loc = np.zeros((2, 2))
    K_loc = np.zeros((2, 2))

    #gaussian quadrature integration
    for k in range(2):
        xi = gauss_pts[k]
        w = gauss_wts[k]
        N_vals, dN_dxi = lagrange_basis(xi)
        #map derivative to physical space
        dN_dx = dN_dxi / J

        #compute local mass and stiffness contributions
        M_loc += np.outer(N_vals, N_vals) * J * w
        K_loc += np.outer(dN_dx, dN_dx) * J * w

    #global matrices
    M[n1:n2+1, n1:n2+1] += M_loc
    K[n1:n2+1, n1:n2+1] += K_loc

#sol vector
u = initial_condition(x_nodes)

#time step loop
for step in range(num_steps):
    t = t_start + step * dt

    #rhs vector (load vector)
    F = np.zeros(N)
    for elem in range(N - 1):
        n1, n2 = elem, elem + 1
        x1, x2 = x_nodes[n1], x_nodes[n2]
        J = (x2 - x1) / 2
        F_loc = np.zeros(2)
        for k in range(2):
            xi = gauss_pts[k]
            w = gauss_wts[k]
            N_vals, _ = lagrange_basis(xi)
            #map local xi to physical x
            x_phys = ((1 - xi) * x1 + (1 + xi) * x2) / 2
            f_val = f(x_phys, t)
            F_loc += N_vals * f_val * J * w
        F[n1:n2+1] += F_loc

    #lhs matrix for implicit Euler: M + dt*K
    A = M + dt * K
    #rhs vector
    b = M @ u + dt * F

    #dirichlet boundary conditions
    A[0, :] = 0
    A[0, 0] = 1
    b[0] = boundary_left(t + dt)
    A[-1, :] = 0
    A[-1, -1] = 1
    b[-1] = boundary_right(t + dt)

    #solve linear system
    u = solve(A, b)

#plot
plt.plot(x_nodes, u, label='Numerical')
plt.plot(x_nodes, np.exp(-t_end) * np.sin(np.pi * x_nodes), '--', label='Analytic')
plt.xlabel('x')
plt.ylabel('u(x, t)')
plt.legend()
plt.title('1D Heat Equation Solution at t = {:.2f}'.format(t_end))
plt.show()
