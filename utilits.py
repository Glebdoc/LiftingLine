import numpy as np
from geometry import Point, Vortex

def spatial_transform(theta, vector, flatten=True):
    vector.transpose()
    Rx = np.array([[1, 0, 0], 
                   [0, np.cos(theta), -np.sin(theta)], 
                   [0, np.sin(theta), np.cos(theta)]])
    if flatten:
        return Rx@vector.flatten()
    else:
        return Rx@vector
    



def assembleLeftVortex(xw, yw, zw, shift, translation):
    left = []
    for i in range(len(xw)-1):
        vector1 = spatial_transform(shift, np.array([xw[i+1], yw[i+1], zw[i+1]])) + translation
        vector2 = spatial_transform(shift,np.array([xw[i], yw[i], zw[i]]))        + translation
        left.append(Vortex(Point(*vector1), Point(*vector2), 1))
    return left

def assembleRightVortex(xw, yw, zw, shift, translation):
    right = []
    for i in range(len(xw)-1):
        vector1 = spatial_transform(shift, np.array([xw[i], yw[i], zw[i]]))       + translation
        vector2 = spatial_transform(shift, np.array([xw[i+1], yw[i+1], zw[i+1]])) + translation
        right.append(Vortex(Point(*vector1), Point(*vector2), 1))
    return right

def spanwise_discretization(N, ROOT_LOCATION, TIP_LOCATION, discretization_type='uniform'):
    if discretization_type == 'uniform':
        return np.linspace(ROOT_LOCATION, TIP_LOCATION, N+1)
    if discretization_type == 'cosine':
        # Generate cosine-spaced points between 0 and π
        theta = np.linspace(np.pi, 0, N+1)
        
        # Cosine-spaced points in the range [-1, 1]
        cos_spaced = np.cos(theta)
        
        # Map these points to the spanwise locations
        half_span = (TIP_LOCATION - ROOT_LOCATION) / 2
        spanwise_locations = ROOT_LOCATION + half_span + (cos_spaced) * half_span
        
        return spanwise_locations
    
def set_axes_equal(ax):
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    """

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])