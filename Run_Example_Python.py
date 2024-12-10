import DynamicRelaxation_Python as Drpy
import numpy as np
import matplotlib.pyplot as plt

'''
This is an example script that shows all the different functions being used in order: 

1) generate a flat, square mesh (with dummy options)
2) deform the position of some fixed nodes (with dummy options)
3) run dynamic relaxation 
4) optionally plot or animate and save as a gif
'''

# ------------------- Generate a mesh ----------------------
node_count= 11
node_length= 100
mesh_type = 'quad'
fixed_node_option = 'all'

verts, edges, fixed, fixed_indices_OG, midlines, diagonals = Drpy.generate_mesh_square(node_count,node_length,mesh_type,fixed_node_option)

# get rest lengths here - assuming that equilibrium is the undeformed mesh
# ATTENTION: this may not be appropriate for your case
edge_vectors = verts[edges[:, 1]] - verts[edges[:, 0]]
rest_lengths = np.linalg.norm(edge_vectors, axis=1, keepdims=True)  # rest lengths of edges


# ------------------- Deform the mesh ----------------------

shape_function = 'user-defined'#'sinusoidal' # can be 'pyramid', 'sinusoidal', 'gaussian', 'user-defined': height list in max_height
PointsDeformed = [102,74] # midlines # can be 'midlines', 'diagonals' or 'user-defined': given as a list of indices
max_height = [70,50] # can be single val for funcitons. single val or list with height per index if user defined
verts, fixed_indices_New, fixed_indices_Total, free_indices_Total = Drpy.apply_shape_function(verts, PointsDeformed, fixed_indices_OG, shape_function, max_height, plot_option=True)

# -------------------- applying any nodewise loading, like gravity -----------------------------

vertexCount = len(verts)
AppliedNodalForces = np.zeros((vertexCount, 3)) #zero force in 3 dimensions
# could optionally add mass * gravity as m*[0,0,-9.81], but given mass is ficticious and assumed lumped
# during the analysis, this may not be meaningful without due care 

# ------------------- Use Dynamic Relaxation to find minimum energy config, post deformation ----------------------

# DR params
params = {
    'MaximumIterations': 150,
    'tolerance': 5e-1,
    'dt': 0.49,  # Time interval
    'm': 1,  # Fictitious mass
    'Elastic': 1000,  # Elastic modulus (MPa)
    'Area': (0.1**2) * np.pi,  # Cross-sectional area (note: use ** for exponentiation)
    'link_pretension': [0],  # link pre-tension in N (between nodes, not at them)
    'restlengthOG': rest_lengths,
}

LinkTension, xyz, KEAllStore, edges, history = Drpy.dynamic_relaxation(params,verts, edges, free_indices_Total, fixed_indices_Total, AppliedNodalForces)

# ------------------------- visualise
fig = plt.figure()
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

# plot Initial geoemtry
Drpy.plot_geometry_nofig(ax1, history['xyz'][0], edges, fixed_indices_Total,'Initial')
# plot final geoemtry
Drpy.plot_geometry_nofig(ax2,xyz, edges, fixed_indices_Total,'Final')

plt.show()

#------------ plot kinetic energy convergence --------
figKE = plt.figure()
plt.plot(history['KE'])
plt.title('KE vs iterations')
plt.show()

# ---------- animate the time history and save a gif -----------

# Example usage
animation = Drpy.animate_geometry(history, edges, fixed_indices_Total)
animation.save("geometry_animation.gif", writer="pillow", fps=3)
print("Animation saved as 'geometry_animation.gif'")


