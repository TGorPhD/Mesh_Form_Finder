import numpy as np

def dynamic_relaxation(params, xyz, edges, free_indices, fixed, applied_nodal_forces):
    """
    Perform dynamic relaxation to find the equilibrium state of a structure.
    Returns:
        LinkTension, xyz, KEAllStore, edges, history
    """

    # Extract parameters
    dt = params['dt']
    area = params['Area']
    E_elastic = params['Elastic']
    restlengthOG = params['restlengthOG']
    tolerance = params['tolerance']
    max_iterations = params['MaximumIterations']
    applied_link_pretension = params['link_pretension']

    # Initialize variables
    vertex_count = xyz.shape[0]
    edge_count = edges.shape[0]
    K_S = (E_elastic * area) / restlengthOG  # Initial geometric stiffness
    #print('KS is:', K_S)
    #print('restlengthOG is:', restlengthOG)

    v = np.zeros((vertex_count, 3))  # Initializing velocities
    S = np.zeros((vertex_count, 3))  # Forces on nodes
    R = np.zeros((vertex_count, 3))  # Residual forces
    KE = 0  # Initial kinetic energy
    reset_flag = False

    # Initialize history
    history = {
        "xyz": [],
        "v": [],
        "KE": [],
        "m": [],
        "iters": 0,
    }

    # Main loop variables
    t = 1
    difference = np.inf
    KE_store = []

    while difference > tolerance and t <= max_iterations:
        # Step 1: Initialization
        vp = np.copy(v)  # Store previous velocities
        xyz0 = np.copy(xyz)  # Store previous positions
        KE_0 = KE  # Store previous kinetic energy

        # Step 2: Compute Link Lengths and Forces
        link_vectors_t = xyz[edges[:, 1]] - xyz[edges[:, 0]]
        link_lengths_t = np.linalg.norm(link_vectors_t, axis=1, keepdims=True)
        #print(link_lengths_t)

        link_tension = np.full_like(link_lengths_t, applied_link_pretension)  # Initialize constant link tension
        link_tension += K_S * (link_lengths_t - restlengthOG)
        link_tension[link_tension < 0] = 0  # No compression allowed

        # Step 3: Node-Specific Mass Computation #
        m = np.zeros(vertex_count)
        geom_stiffness = link_tension / link_lengths_t
        #print(geom_stiffness)
        for node in range(vertex_count):
            adjoining_edges = np.sum(edges == node, axis=1)
            #print(adjoining_edges)
            adjoining_edges = adjoining_edges[:, np.newaxis] 
            stiffness = np.sum((K_S / restlengthOG) * adjoining_edges + geom_stiffness * adjoining_edges)
            #print('stiffness', stiffness)
            #print(K_S.shape, restlengthOG.shape, geom_stiffness.shape, adjoining_edges.shape)
            m[node] = stiffness
        

        # Step 4/5: Force Resolution
        edge_force_s = np.zeros((edge_count, 3))
        for i in range(edge_count):
            edge_force_s[i] = (link_vectors_t[i] / link_lengths_t[i]) * link_tension[i]
            S[edges[i, 0]] += edge_force_s[i]
            S[edges[i, 1]] -= edge_force_s[i]

        R[free_indices] = applied_nodal_forces[free_indices] + S[free_indices]

        # Step 6: Update Velocities and Positions
        A = 1  # Damping constant
        B = (dt / m[:, None])

        v[free_indices] = A * vp[free_indices] + B[free_indices] * R[free_indices]
        xyz[free_indices] = xyz0[free_indices] + dt * v[free_indices]

        KE = np.sum(0.5 * m * np.linalg.norm(v, axis=1)**2)
        KE_store.append(KE)

        # Step 7/8: Kinetic Energy Correction
        if t > 2 and not reset_flag and KE > KE_0:
            E = KE_store[-2] - KE_store[-1]
            D = KE_store[-3] - KE_store[-2]
            q = max(0, min(E / (E - D), 1))  # Clamp q between 0 and 1

            R_t = m[:, None] * (v - vp) / dt
            xyz -= (dt * (1 + q) * vp) + (dt**2 / 2) * q * (R_t / m[:, None])

            # Reset variables
            KE = 0
            v = np.zeros((vertex_count, 3))
            R = np.zeros((vertex_count, 3))
            S = np.zeros((vertex_count, 3))
            reset_flag = True
        else:
            reset_flag = False

        # Step 9: Convergence Check
        difference = np.linalg.norm(xyz - xyz0, ord='fro')
        history["xyz"].append(np.copy(xyz))
        history["v"].append(np.copy(v))
        history["KE"].append(KE)
        history["m"].append(m)

        t += 1

    history["iters"] = t - 1
    return link_tension, xyz, KE_store, edges, history