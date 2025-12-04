import numpy as np

# Number of elements
num_elements = 10

# Number of nodes (linear elements: nodes = elements + 1)
num_nodes = num_elements + 1

# Domain length
L = 1.0

# Node coordinates in x-space
x = np.linspace(0, L, num_nodes)

# Element connectivity (each row: [node1, node2])
elements = np.array([[i, i+1] for i in range(num_elements)])

# Initialize global mass and stiffness matrices
M = np.zeros((num_nodes, num_nodes))  # Mass matrix
K = np.zeros((num_nodes, num_nodes))  # Stiffness matrix

# Gauss quadrature points and weights for 2-point rule in parent space [-1, 1]
gauss_pts = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
gauss_wts = np.array([1.0, 1.0])

# Loop over each element to assemble local matrices and add to global
for e in range(num_elements):
    # Get global node indices for this element
    n1, n2 = elements[e]
    
    # Get x-coordinates of the nodes
    x1, x2 = x[n1], x[n2]
    
    # Compute the Jacobian of the mapping from parent space to x-space
    J = (x2 - x1) / 2.0  # dx/dξ for linear mapping
    
    # Initialize local mass and stiffness matrices (2x2 for linear elements)
    Me = np.zeros((2, 2))
    Ke = np.zeros((2, 2))
    
    # Loop over Gauss points for numerical integration
    for i in range(len(gauss_pts)):
        xi = gauss_pts[i]      # Gauss point in parent space
        w = gauss_wts[i]       # Corresponding weight
        
        # Shape functions evaluated at xi
        N = np.array([0.5 * (1 - xi), 0.5 * (1 + xi)])
        
        # Derivatives of shape functions w.r.t. parent coordinate ξ
        dN_dxi = np.array([-0.5, 0.5])
        
        # Derivatives w.r.t. x: dN/dx = dN/dξ * dξ/dx = dN/dξ / J
        dN_dx = dN_dxi / J
        
        # Compute local mass matrix contribution at this Gauss point
        Me += np.outer(N, N) * J * w
        
        # Compute local stiffness matrix contribution at this Gauss point
        Ke += np.outer(dN_dx, dN_dx) * J * w
    
    # Assemble local matrices into global matrices
    # Add Me and Ke to the correct positions in M and K
    indices = [n1, n2]
    for a in range(2):
        for b in range(2):
            M[indices[a], indices[b]] += Me[a, b]
            K[indices[a], indices[b]] += Ke[a, b]

# Print the assembled global mass and stiffness matrices
print("Global Mass Matrix M:\n", M)
print("Global Stiffness Matrix K:\n", K)
