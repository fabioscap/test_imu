# takes input csv files
# time x y z r p y
# produces trajectory plot
# different files have different colors

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import argparse

def read_csv(file_path):
    col_names = ['t', 'x', 'y', 'z', 'a', 'b', 'c']
    return pd.read_csv(file_path, sep=",", names=col_names)

def rotation_matrix(a, b, c):
    Rx = np.array([[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]])
    Ry = np.array([[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]])
    Rz = np.array([[math.cos(c), -math.sin(c), 0], [math.sin(c), math.cos(c), 0], [0, 0, 1]])
    return np.dot(Rz, np.dot(Ry, Rx))

def plot_axes(ax, pos, R, scale=1.0):
    axes = np.eye(3)
    origin = np.array([pos]).T

    transformed_axes = np.dot(R, axes) * scale

    # Plotting each axis
    for i, color in zip(range(3), ['k', 'g', 'b']):
        ax.plot([origin[0], origin[0] + scale*R[0,i]], 
                [origin[1], origin[1] + scale*R[1,i]], 
                #[origin[2], origin[2] + transformed_axes[2, i]], 
                color=color)
        break

def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def plot_poses(df, ax, plot_poses=False,color="r",label="trajectory"):

    positions = np.zeros((2,len(df)))
    #positions = np.zeros((3,len(df)))
    for i in range(len(df)):
        # pos = df.loc[i, ['x', 'y', 'z']]
        pos = df.loc[i, ['x', 'y']]
        positions[:,i] = np.array(pos)
        a, b, c = df.loc[i, ['a', 'b', 'c']]
        R = rotation_matrix(a, b, c)
        if plot_poses:
          plot_axes(ax, pos, R, scale=3.0) 

    ax.plot(positions[0,:], 
            positions[1, :], 
            #positions[2,:],
            label=label, color=color,linewidth=2)




def main():
    parser = argparse.ArgumentParser(description="Plot 3D poses from CSV files.")
    parser.add_argument('files', nargs='+', help="Paths to CSV files")
    args = parser.parse_args()

    colors = ['green', 'red', 'blue', 'purple', 'orange']

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    fig, ax = plt.subplots(figsize=[16,9])

    for file_path, color in zip(args.files, colors):
        df = read_csv(file_path)
        plot_poses(df, ax, color=color, plot_poses=True, label=file_path)

    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    # ax.set_zlabel('Z [m]')
    # set_axes_equal(ax)
    ax.axis('equal')
    # ax.view_init(90, 0)
    plt.legend(fontsize=25)
    plt.grid(True)


    plt.show()
    

if __name__ == "__main__":
    main()