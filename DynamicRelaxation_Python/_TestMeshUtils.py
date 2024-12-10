import numpy as np
import matplotlib.pyplot as plt

def generate_mesh_square(node_count, length, mesh_type='quad', fixed_node_option='all'):
    """
    Generate a square mesh with specified properties.

    Parameters:
        node_count (int): Number of nodes per side (>= 2).
        length (float): Total length of the square side (> 0).
        mesh_type (str): 'quad' for quadrilateral, 'tri' for triangular, or 'tri_both' for both diagonals in each quad.
        fixed_node_option (str): Controls fixed nodes:
            'all' - All edge nodes fixed.
            'everyOther' - Every other edge node fixed.
            'none' - No nodes fixed.

    Returns:
        verts (numpy.ndarray): Nodal coordinates (n x 3 array).
        edges (numpy.ndarray): Connectivity matrix (m x 2 array).
        fixed (numpy.ndarray): Fixed nodes (p x 3 array).
        fixed_indices (numpy.ndarray): Indices of fixed nodes.
        midlines (dict): Indices for horizontal and vertical midlines.
        diagonals (dict): Indices for main and anti-diagonals.
    """
    # Input validation
    if node_count < 2:
        raise ValueError("node_count must be at least 2.")
    if length <= 0:
        raise ValueError("length must be greater than 0.")
    if mesh_type not in {'quad', 'tri', 'tri_both'}:
        raise ValueError('mesh_type must be "quad", "tri", or "tri_both".')
    if fixed_node_option not in {'all', 'everyOther', 'none'}:
        raise ValueError('fixed_node_option must be "all", "everyOther", or "none".')

    # Step 1: Create nodes
    spacing = length / (node_count - 1)
    x, y = np.meshgrid(np.linspace(0, length, node_count), np.linspace(0, length, node_count))
    verts = np.column_stack((x.ravel(), y.ravel(), np.zeros(node_count**2)))  # All nodes in the z=0 plane initially

    # Step 2: Define connectivity
    edges = []
    for i in range(node_count - 1):
        for j in range(node_count - 1):
            # Get indices of the current square's corners
            n1 = i * node_count + j
            n2 = n1 + 1
            n3 = n1 + node_count
            n4 = n3 + 1

            if mesh_type == 'quad':
                # Quadrilateral connectivity
                edges.extend([(n1, n2), (n2, n4), (n4, n3), (n3, n1)])
            elif mesh_type == 'tri':
                # Triangular connectivity
                edges.extend([(n1, n2), (n2, n4), (n4, n1), (n1, n3), (n3, n4)])
            elif mesh_type == 'tri_both':
                # Triangular connectivity (both diagonals per quad)
                edges.extend([(n1, n2), (n2, n4), (n4, n1)])  # Triangle 1
                edges.extend([(n2, n3), (n3, n4), (n4, n2)])  # Triangle 2

    # Remove duplicate edges
    edges = np.unique(np.sort(edges, axis=1), axis=0)

    # Step 3: Determine fixed nodes
    edge_indices = np.unique(
        np.concatenate([
            np.arange(node_count),  # Bottom edge
            np.arange(node_count - 1, node_count**2, node_count),  # Right edge
            np.arange(node_count**2 - 1, node_count**2 - node_count - 1, -1),  # Top edge
            np.arange(node_count**2 - node_count, -1, -node_count)  # Left edge
        ])
    )

    if fixed_node_option == 'all':
        fixed_indices = edge_indices
    elif fixed_node_option == 'everyOther':
        fixed_indices = edge_indices[::2]
    else:  # 'none'
        fixed_indices = np.array([], dtype=int)
    
    fixed = verts[fixed_indices]

    # Step 4: Identify midlines and diagonals
    midlines = {
        'horizontal': np.where(np.abs(verts[:, 1] - length / 2) < spacing / 2)[0],
        'vertical': np.where(np.abs(verts[:, 0] - length / 2) < spacing / 2)[0]
    }

    diagonals = {
        'main': np.where(np.abs(verts[:, 0] - verts[:, 1]) < spacing / 2)[0],
        'anti': np.where(np.abs(verts[:, 0] + verts[:, 1] - length) < spacing / 2)[0]
    }

    return verts, edges, fixed, fixed_indices, midlines, diagonals

def apply_shape_function(verts, lines, fixed_indices_OG, shape_function, max_height, plot_option=False):
    """
    Applies a shape function to the vertices along specified lines.

    Parameters:
        verts (numpy.ndarray): Vertex coordinates (n x 3 array).
        lines (dict): Dictionary with diagonal indices or other line definitions.
        shape_function (str): Shape function ('pyramid', 'gaussian', 'sinusoidal', 'user-defined').
        max_height (float): Maximum height for the shape function.
        plot_option (bool): If True, plot before and after the transformation.

    Returns:
        updated_verts (numpy.ndarray): Modified vertex coordinates (n x 3 array).
        new_fixed_indices (numpy.ndarray): Indices of updated vertices.
    """
    # Validate shape function
    valid_shapes = {'pyramid', 'gaussian', 'sinusoidal','user-defined'}
    if shape_function not in valid_shapes:
        raise ValueError(f"Invalid shape function. Choose from {valid_shapes}.")

    # Combine all line indices into a single list
    if shape_function == 'user-defined':
        line_indices = lines
    else:
        line_indices = np.unique(np.concatenate([lines[key] for key in lines]))

    fixed_indices_New = line_indices

    # Extract vertices corresponding to the line indices
    line_verts = verts[line_indices]

    # Calculate distances from the line's center for each vertex
    line_center = np.mean(line_verts[:, :2], axis=0)  # Center in the x-y plane
    distances = np.linalg.norm(line_verts[:, :2] - line_center, axis=1)

    # Normalize distances to [0, 1]
    max_distance = np.max(distances)
    normalized_distances = distances / max_distance

    # Compute new z values based on the selected shape function
    if shape_function == 'pyramid':
        # Linear height distribution (inverted V shape)
        new_z = max_height * (1 - normalized_distances)
    elif shape_function == 'gaussian':
        # Gaussian curve
        new_z = max_height * np.exp(-4 * (normalized_distances ** 2))
    elif shape_function == 'sinusoidal':
        # Half sinusoidal wave (arch-like shape)
        new_z = max_height * (np.sin(np.pi * (1 - normalized_distances) / 2) ** 2)
    elif shape_function == 'user-defined':
        # max height defined manually
        new_z = max_height 
        pass
    else:
        raise ValueError("Unsupported shape function.")

    # Assign new z values to the vertices
    updated_verts = verts.copy()
    updated_verts[line_indices, 2] = new_z

    #get new total fixed indices
    fixed_indices_Total = np.unique(np.hstack((fixed_indices_OG, fixed_indices_New)))
    AllIndices = (np.arange(len(verts)))
    free_indices_Total = np.array(list((set(AllIndices) | set(fixed_indices_Total)) - (set(AllIndices) & set(fixed_indices_Total))))
    #print(AllFreeIndices)

    # Plot before and after if plot_option is enabled
    if plot_option:
        fig = plt.figure(figsize=(12, 6))

        # Original vertices
        ax1 = fig.add_subplot(121, projection='3d')
        ax1.scatter(verts[:, 0], verts[:, 1], verts[:, 2], c='b', s=20, label="Original free nodes")
        ax1.scatter(line_verts[:, 0], line_verts[:, 1], line_verts[:, 2], c='r', s=50, label="Original fixed nodes")
        ax1.set_title("Original Mesh")
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        ax1.set_zlabel("Z")
        ax1.grid(True)
        ax1.legend()
        ax1.view_init(elev=20, azim=30)

        # Updated vertices
        ax2 = fig.add_subplot(122, projection='3d')
        ax2.scatter(updated_verts[:, 0], updated_verts[:, 1], updated_verts[:, 2], c='b', s=20, label="Updated free nodes (unchanged)")
        ax2.scatter(updated_verts[line_indices, 0], updated_verts[line_indices, 1], updated_verts[line_indices, 2],
                    c='r', s=50, label="Updated fixed nodes (changed)")
        ax2.set_title(f"Mesh After Applying {shape_function.capitalize()} Shape")
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        ax2.set_zlabel("Z")
        ax2.grid(True)
        ax2.legend()
        ax2.view_init(elev=20, azim=30)

        plt.tight_layout()
        plt.show()

    return updated_verts, fixed_indices_New, fixed_indices_Total, free_indices_Total