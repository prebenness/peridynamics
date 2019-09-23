from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plot_coors(X):
    max, min = np.amax(X), np.amin(X)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X[:,0], X[:,1], X[:,2])
    ax.set_xlim(max, min)
    ax.set_ylim(max, min)
    ax.set_zlim(max, min)
    plt.show()
