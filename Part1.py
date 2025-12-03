import numpy

#problem set up
# spatial domain (0,1) + time domain (0,1)
# interval divided into N-1 sections
#   nodes @ x_i = i(h) -> h = 1/(N-1)
# inital conditon -> u(x,0)=sin(pix)
#   create sol vector u_i^0=sin(pix_i) [starting pt for time integration]
# dirichlet boundary conditions
#   set sol to zero at boundary nodes for all time steps -> u(0,t)=u(1,t)=0
#     change global matrix to zero except for diagional (set = to 1, w/r.h.s = 0 for those nodes)
# define source funct.
#   eval source @ each time step -> contribute to r.h.s of system
#   allow for general f(x,t) -> make general funct. def


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


#time integration
#backward vs. forward euler


#generalization


#output/visual


#test/validation

