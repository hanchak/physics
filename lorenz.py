'''
Lorenz strange attractor

Mike Hanchak
04OCT2017
'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  
import matplotlib.animation as animation

# Lorenz parameters
sigma = 10
beta = 8/3
rho = 28

N = 10000  # number of points to simulate
L = 1000 # length of tail
pts = np.zeros((3, N+1))  # use numpy array to store all of the data
dt = 0.01  # pseudo-time step
xlim = (-15,15)
ylim = xlim
zlim = (15,45)

# initial conditions
pts[0,0] = 1
pts[2,0] = 10

# set up the figure
fig = plt.figure(figsize=(6,6))
fig.canvas.set_window_title('Lorenz Attractor: left click and drag to rotate')
ax = fig.add_subplot(111, projection='3d', facecolor='k')
ax.set_xlim3d(xlim)
ax.set_ylim3d(ylim)
ax.set_zlim3d(zlim)
plt.axis('off')
ax.set_position([0,0,1,1])

# this is a blank plot; we will change its data in the draw() function
scatter, = ax.plot([], [], [], 'c')


def draw(j):
    
    # get the current point:
    x = pts[0,j]
    y = pts[1,j]
    z = pts[2,j]
    
    # get derivatives (Lorenz)
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    
    # Euler integration step
    x += dx * dt
    y += dy * dt
    z += dz * dt
    
    # update array of points
    pts[:,j+1] = np.array([x,y,z])
    
    # update plot data
    if j < L:
        scatter.set_data(pts[0,:j], pts[1,:j])
        scatter.set_3d_properties(pts[2,:j])
    else:
        scatter.set_data(pts[0,j-L:j], pts[1,j-L:j])
        scatter.set_3d_properties(pts[2,j-L:j])  
    
    scatter.set_color([np.cos(j/L*2)**2, np.cos(j/L*3)**2, np.cos(j/L*5)**2, 1])
    #return [scatter]  # only needed if blit=True


ani = animation.FuncAnimation(fig, draw, np.arange(0,N), interval=0)

plt.show()
