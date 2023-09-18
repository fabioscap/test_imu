import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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

# Load poses_pred from the text file
poses_pred = []
with open("../examples/output_pred.txt", 'r') as file:
    lines = file.readlines()
    for i in range(0, len(lines), 4):
        T = np.empty((4, 4), dtype=float)
        string_list = [lines[i],lines[i+1],lines[i+2],lines[i+3]]
        for i, line in enumerate(string_list):
          values = line.strip().split()

          T[i, :] = [float(value) for value in values]
        R = T[:3, :3]
        t = T[:3, 3]
        poses_pred.append((R,t))

poses_gt = []
with open("../examples/output_gt.txt", 'r') as file:
    lines = file.readlines()
    for i in range(0, len(lines), 4):
        T = np.empty((4, 4), dtype=float)
        string_list = [lines[i],lines[i+1],lines[i+2],lines[i+3]]
        for i, line in enumerate(string_list):
          values = line.strip().split()

          T[i, :] = [float(value) for value in values]
        R = T[:3, :3]
        t = T[:3, 3]
        poses_gt.append((R,t))


# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Plot poses_pred
""" for p in range(len(poses_pred)):
    R,t = poses_pred[p]
    R_gt, _ = poses_gt[p]
    # Define axes colors
    colors = ['r', 'g', 'b']

    # Plot the X, Y, and Z axes of the reference frame
    i=1
    axis = R.T[1]
    ax.quiver(
        t[0], t[1], t[2], axis[0], axis[1], axis[2],
        color="blue", label=f'Poses Pred Axis {i + 1}', length=0.1
    )
    axis = R_gt.T[1]
    ax.quiver(
        t[0], t[1], t[2], axis[0], axis[1], axis[2],
        color="red", label=f'Poses Pred Axis {i + 1}', length=0.1
    ) """
    

# Extract trajectory points (translation vectors) for poses_pred
trajectory_pred = np.array([t for _, t in poses_pred])

# Plot the whole trajectory for poses_pred
ax.plot(
    trajectory_pred[:, 0], trajectory_pred[:, 1], trajectory_pred[:, 2],
    label='Trajectory Pred', color='blue'  # Choose a color for poses_pred
)

# Plot poses_gt
""" for R, t in poses_gt:
    # Plot the X, Y, and Z axes of the reference frame for poses_gt
    for i, axis in enumerate(R.T):
        colors = ['r', 'g', 'b']
        ax.quiver(
            t[0], t[1], t[2], axis[0], axis[1], axis[2],
            color=colors[i], label=f'Poses GT Axis {i + 1}', length=0.5
        ) """

# Extract trajectory points (translation vectors) for poses_gt
trajectory_gt = np.array([t for _, t in poses_gt])

# Plot the whole trajectory for poses_gt
ax.plot(
    trajectory_gt[:, 0], trajectory_gt[:, 1], trajectory_gt[:, 2],
    label='Trajectory GT', color='red'  # Red color for poses_gt
)

# Set equal axis scales and labels
ax.set_box_aspect([1, 1, 1])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
set_axes_equal(ax)


# Show the plot
plt.show()