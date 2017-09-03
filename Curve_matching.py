from __future__ import division
import matplotlib as mpl
import numpy as np
from scipy.misc import comb
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import math
import icp

if __name__ == "__main__":
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #ax = fig.add_subplot(111, projection='3d')
    x = np.linspace(0.5,10.0,50)
    y = [1/t for t in x]
    z = [i for i in range(0,50)]
    ax.plot(x,y,z,label='target')
    B = np.array([x,y,z]).T

    x_1 = x + 0.05 * np.random.normal(0, 1, len(x))
    p = [t+1 for t in x_1]
    y_1 = y + 0.05 * np.random.normal(0, 1, len(y))
    q = [t+1 for t in y_1]
    z_1 = z + 0.05 * np.random.normal(0, 1, len(z))
    s = [t+1 for t in z_1]

    #ax.plot(p,q,s,label='trial')

    comp = np.array([p,q,s]).T
    print comp.shape
    theta1 = ( 220.0 / 360) * 2 * np.pi
    rot1 = np.array([[math.cos(theta1), -math.sin(theta1), 0],
                   [math.sin(theta1),  math.cos(theta1), 0], [0, 0, 1]])
    x = []
    y = []
    z = []
    for t in comp:
        f = np.dot(rot1, t)
        x.append(f[0])
        y.append(f[1])
        z.append(f[2])
    
    ax.plot(x,y,z,label='source')

    A = np.array([x,y,z]).T
    print 'A:', A
    T,d = icp.icp(A,B)
    print 'T: ',T
    R_f = T[0:3, 0:3]
    t_f = T[0:3, 3]

    A_f = np.dot(R_f, A.T).T + t_f
    ax.plot(A_f.T[0],A_f.T[1],A_f.T[2],label='registered')
    ax.legend()
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()
