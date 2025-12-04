import numpy as np

#elements
num_elements = 10

#nodes (linear elements: nodes = elements + 1)
num_nodes = num_elements + 1

#length
L = 1.0

#coordinates in x-space
x = np.linspace(0, L, num_nodes)

#element (each row: [node1, node2])
elements = np.array([[i, i+1] for i in range(num_elements)])

#global mass and stiffnesss
M = np.zeros((num_nodes, num_nodes))  #mass
K = np.zeros((num_nodes, num_nodes))  #xtiffness

#gauss quadrature points + weights for 2-point rule in parent space 
gauss_pts = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
gauss_wts = np.array([1.0, 1.0])

#loop over each element to local matrices and add to global
for e in range(num_elements):
    # Get global node indices for this element
    n1, n2 = elements[e]
    
    #x-coordinates of the nodes
    x1, x2 = x[n1], x[n2]
    
    #jacobian of the mapping from parent space to x-space
    J = (x2 - x1) / 2.0  # dx/dξ for linear mapping
    
    #local mass and stiffness matrices (2x2 for linear elements)
    Me = np.zeros((2, 2))
    Ke = np.zeros((2, 2))
    
    #loop over Gauss points for numerical integration
    for i in range(len(gauss_pts)):
        xi = gauss_pts[i]      #gp in parent space
        w = gauss_wts[i]       #weight
        
        #shape functions eval
        N = np.array([0.5 * (1 - xi), 0.5 * (1 + xi)])
        
        #derivatives of shape functions w.r.t. parent coordinate
        dN_dxi = np.array([-0.5, 0.5])
        
        #derivatives w.r.t. x: dN/dx = dN/ds * ds/dx = dN/ds / J
        dN_dx = dN_dxi / J
        
        #local mass matrix fauss point
        Me += np.outer(N, N) * J * w
        
        #tiffness matrix gauss point
        Ke += np.outer(dN_dx, dN_dx) * J * w
    
    #local into global matrices
    #d Me and Ke to the correct in M and K
    indices = [n1, n2]
    for a in range(2):
        for b in range(2):
            M[indices[a], indices[b]] += Me[a, b]
            K[indices[a], indices[b]] += Ke[a, b]

#global mass and stiffness 
print("Global Mass Matrix M:\n", M)
print("Global Stiffness Matrix K:\n", K)
