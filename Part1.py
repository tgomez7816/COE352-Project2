import numpy as np

#problem set up
# spatial domain (0,1) + time domain (0,1)
x_s, x_n = 0.0, 1.0 
t_s, t_n = 0.0, 1.0
# interval divided into N-1 sections
N = 11 #easy variable to change
x = np.linspace(x_s, x_n)
h = (x_n - x_s)/(N-1) # equally spaced nodes
# inital conditon -> u(x,0)=sin(pix)
def init_con(x):
    return #enter sin(pix)
u_0 = init_con(x)

#   create sol vector u_i^0=sin(pix_i) [starting pt for time integration]
# dirichlet boundary conditions
def bc_l(t):

def bc_r(t):

#   set sol to zero at boundary nodes for all time steps -> u(0,t)=u(1,t)=0
#     change global matrix to zero except for diagional (set = to 1, w/r.h.s = 0 for those nodes)
# define source funct.
#   eval source @ each time step -> contribute to r.h.s of system
#   allow for general f(x,t) -> make general funct. def
def f(x, t): #source function easily changable when defined separately 


#basic functions + mesh

# create mesh/grid for domain using N, nodes (starting case N=11)
#   discretize domain into elements
#   Uniform mesh nodes -> x_i = x_0+i(h) w/ h = (x_end - x_0)/(N-1)
#     element defined by 2 nodes (e_j is bounded by x_i and x_i+1)
# construct 1D lagrange basis funct.s for nodes
# map each element from physical space to the parent space for integration
#   linear mapping -> parent coordinate squiggle included on [-1,1]
#   jacobian
#   isoparametric mapping


#element + matrices creation
# local mass entry
# locak stiffness entry
# 2nd order Gaussian quadrature
#   reference element mapping -> apply quadrature to matrix entries
#   global mass + stiffness matrix


#time integration
# backward euler -> better for stiff eqs
# at each time step
#   apply inital condition t=0
#     set sol vector u^0 
#   apply dirichlet boundary conditions
#     mod matrix -> all entries @row k to 0 except A_kk = 1
#     set corresponding entry in r.h.s to boundary value
#   integrate r.h.s -> load vector + gaussian quadrature
#   update solution using backward euler 


#generalization


#output/visual


#test/validation

