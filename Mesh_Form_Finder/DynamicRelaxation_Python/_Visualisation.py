import matplotlib.pyplot as plt
import numpy as np

def plot_geometry_genfig(xyz, edges, fixed_indices, titlePlot = str):
    """
    Plot geometry in 3D.

    Parameters:
    - xyz: np.array of shape (N, 3), node coordinates.
    - edges: np.array of shape (M, 2), indices of edges connecting nodes.
    - fixed: np.array of shape (P, 3), coordinates of fixed nodes.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    fixed = xyz[fixed_indices]
    
    # Plot edges
    for edge in edges:
        x_vals = [xyz[edge[0], 0], xyz[edge[1], 0]]
        y_vals = [xyz[edge[0], 1], xyz[edge[1], 1]]
        z_vals = [xyz[edge[0], 2], xyz[edge[1], 2]]
        ax.plot(x_vals, y_vals, z_vals, 'b')  # Blue lines for edges
    
    
    # Plot free nodes
    #ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c='k', marker='o', label='Nodes')
    
    # Plot fixed nodes
    ax.scatter(fixed[:, 0], fixed[:, 1], fixed[:, 2], c='r', marker='*', label='Fixed Nodes')

    # Set equal scaling by adjusting the limits
    max_range = np.max(np.ptp(fixed, axis=0))  # Find the max range in x, y, z
    mid_x = np.mean(fixed[:, 0])
    mid_y = np.mean(fixed[:, 1])
    mid_z = np.mean(fixed[:, 2])

    # Set the same limit for all axes
    ax.set_xlim(mid_x - max_range / 2, mid_x + max_range / 2)
    ax.set_ylim(mid_y - max_range / 2, mid_y + max_range / 2)
    ax.set_zlim(mid_z - max_range / 2, mid_z + max_range / 2)

    
    # Labels, title, and grid
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(titlePlot)
    ax.legend()
    ax.grid(True)
    ax.set_box_aspect([1, 1, 1])  # Equal scaling
    
    plt.show()

def plot_geometry_nofig(ax, xyz, edges, fixed_indices, titlePlot=""):
    """
    Plot geometry in 3D on the given axes.

    Parameters:
    - ax: Matplotlib 3D axes object where the plot will be drawn.
    - xyz: np.array of shape (N, 3), node coordinates.
    - edges: np.array of shape (M, 2), indices of edges connecting nodes.
    - fixed_indices: np.array of indices of fixed nodes in `xyz`.
    - titlePlot: str, title of the plot.
    """
    ax.clear()  # Clear the previous frame for animation
    
    fixed = xyz[fixed_indices]
    
    # Plot edges
    for edge in edges:
        x_vals = [xyz[edge[0], 0], xyz[edge[1], 0]]
        y_vals = [xyz[edge[0], 1], xyz[edge[1], 1]]
        z_vals = [xyz[edge[0], 2], xyz[edge[1], 2]]
        ax.plot(x_vals, y_vals, z_vals, 'b')  # Blue lines for edges
    
    # Plot fixed nodes
    ax.scatter(fixed[:, 0], fixed[:, 1], fixed[:, 2], c='r', marker='*', label='Fixed Nodes')

    # Set equal scaling by adjusting the limits
    max_range = np.max(np.ptp(fixed, axis=0))  # Find the max range in x, y, z
    mid_x = np.mean(fixed[:, 0])
    mid_y = np.mean(fixed[:, 1])
    mid_z = np.mean(fixed[:, 2])

    # Set the same limit for all axes
    ax.set_xlim(mid_x - max_range / 2, mid_x + max_range / 2)
    ax.set_ylim(mid_y - max_range / 2, mid_y + max_range / 2)
    ax.set_zlim(mid_z - max_range / 2, mid_z + max_range / 2)

    # Labels, title, and grid
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(titlePlot)
    ax.legend()
    ax.grid(True)
    ax.set_box_aspect([1, 1, 1])  # Equal scaling

def animate_geometry(history, edges, fixed_indices_Total):
    '''
    This function creates an animated visualization of geometry in 3D. 
    It uses a history of node coordinates (`history['xyz']`) and edges to animate changes over iterations.
    
    Parameters:
    - history: dict, contains the history of iterations with key 'xyz' for coordinates and 'iters' for frame count.
    - edges: np.array, shape (M, 2), defines the connectivity of nodes for edges.
    - fixed_indices_Total: np.array or list of int, indices of fixed nodes.
    
    Returns:
    - anim: FuncAnimation object, the generated animation.
    '''
    from matplotlib.animation import FuncAnimation
    from functools import partial
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    frames = history['iters']
      # Total number of frames in the animation

    # Ensure fixed_indices_Total is a NumPy array of integers
    fixed_indices_Total = np.array(fixed_indices_Total, dtype=int)

    # Data generator
    def data_gen():
        for i in range(frames):
            yield history['xyz'][i]  # Yield coordinates for each frame

    # Update function
    def update(frame, ax, edges):
        ax.clear()  # Clear previous frame
        plot_geometry_nofig(ax, frame, edges, fixed_indices_Total, titlePlot="Animated Geometry")

    # Create the animation using `partial` to pass additional arguments
    anim = FuncAnimation(fig, partial(update, ax=ax, edges=edges), frames=data_gen, interval=200, save_count=frames)
    
    return anim
