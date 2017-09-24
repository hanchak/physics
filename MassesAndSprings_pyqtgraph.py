'''
Bouncing Ball code using object oriented programming. (3D Version!)

This code uses numpy arrays for all vector math

Code originally from here:
https://github.com/cpbotha/bwtl-python-tutorials/blob/master/
  part%205%20-%20object%20oriented%20programming%20and%20bouncing%20balls.ipynb
'''
#%% IMPORTS
import numpy as np
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
from random import random

# gravitational acceleration on Earth in m*s^-2
g = 0  #9.81
Ks = 50.  # spring constant
Kd = 0.3   # dashpot constant (viscous)
Kv = 0.2  # ball damping (viscous)
eqdist = 0.25  # equilibrium distance of spring
m = 1  # mass of ball

Nballs = 8

# acceleration vector due to g
ag = np.array([0., 0., -g]).reshape((1,3))
# coefficient of restitution (ratio of velocity after and before bounce)
cor = 0.55

# bounds of the room
xlim = [-0.5, 0.5]
ylim = xlim
zlim = xlim  #[0., 2.]

#  delta t
delta_t = 0.002

# create figure
app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.setBackgroundColor('k')
w.opts['distance'] = 1
w.show()
w.setWindowTitle('Masses and Springs')

#g = gl.GLGridItem()
#w.addItem(g)
#g = gl.GLAxisItem()
#w.addItem(g)

#%% DEFINE CLASS
class Ball():
    #global w
    '''This class creates a ball and provides for physics of motion'''
    def __init__(self, xy, v):
        #global w
        """:param xy: Initial position, :param v: Initial velocity. """
        # cast inputs as 2-d numpy arrays
        self.xy = xy.reshape((1,3)) #np.array(xy)
        self.v = v.reshape((1,3))  #np.array(v)
        
        # generate sphere data and colors.
        md = gl.MeshData.sphere(rows=10, cols=20, radius=0.05)
        color = np.random.randint(1,5, (1,4)) /4
        #color[0,3] = 1
        colors = np.repeat(color, md.faces().shape[0], axis=0)
        md.setFaceColors(colors)
        
        g = gl.GLMeshItem(meshdata=md, smooth=True)
        g.setShader('shaded')
        g.translate(self.xy[0,0],self.xy[0,1],self.xy[0,2])
        
        w.addItem(g)
        self.scatter = g

    def update(self, others):
        '''update the position from equations of motion & update plot data'''
        
        if self.xy[0,0] <= xlim[0] or self.xy[0,0] >= xlim[1]:
            # hit the left or right wall, reflect x component
            self.v[0,0] *= - cor

        if self.xy[0,1] <= ylim[0] or self.xy[0,1] >= ylim[1]:
            # hit the bottom or top wall, reflect y component
            self.v[0,1] *= - cor
            
        if self.xy[0,2] <= zlim[0] or self.xy[0,2] >= zlim[1]:
            # hit the bottom or top wall, reflect y component
            self.v[0,2] *= - cor
        
        
        # loop through all other balls to calcuate attractive force
        force = np.zeros((1,3))
        
        for other in others:

            vec = (other.xy - self.xy).reshape((1,3))  # position vec of other w.r.t. selfff
            
            dist = np.sqrt(vec[0,0]**2 + vec[0,1]**2 + vec[0,2]**2)  # scalar distance
                      
            if dist > 0.0001:
                unitvec = (vec / dist).reshape((1,3))  # unit vector from self to other

                vec_rel = (other.v - self.v).reshape(3).dot(unitvec.reshape(3))  #rel vel along unit vec
                
                # spring force between balls
                force +=  (dist - eqdist) * Ks * unitvec
                # damping force between balls
                force +=  vec_rel * Kd * unitvec
                
            
        # add damping force to the ball proportional to absolute velocity
        dforce = - Kv * self.v
        
        # Euler integration: v(i+1) = v(i) + acc*dt
        self.v += (dforce/m + force/m + ag) * delta_t

        # Euler integration: x(i+1) = x(i) + vel*dt
        dxy = self.v * delta_t
        self.xy += dxy
        
        # keep xs and ys inside box:
        #self.xy[0] = np.clip(self.xy[0], xlim[0], xlim[1])
        #self.xy[1] = np.clip(self.xy[1], ylim[0], ylim[1])
        #self.xy[2] = np.clip(self.xy[2], zlim[0], zlim[1])

        # update the scatter plot with the new x and y positions:
        self.scatter.translate(dxy[0,0],dxy[0,1],dxy[0,2])
        


#%% CONSTRUCT OBJECTS
balls = []
for i in range(Nballs):
    startxy = np.array([random()*(xlim[1]-xlim[0])+ xlim[0], \
                      random()*(ylim[1]-ylim[0])+ ylim[0], \
                      random()*(zlim[1]-zlim[0])+ zlim[0]])
    
    balls.append( Ball( startxy, np.zeros(3) ) )


#%% ANIMATION ROUTINES
def animate():
    for ball in balls:  # loop through every ball
        ball.update(balls)  # call each ball's update method

		
#%% CALL ANIMATION
t = QtCore.QTimer()
t.timeout.connect(animate)
t.start(2)  # frame delay in ms

import sys
if sys.flags.interactive != 1 or not hasattr(QtCore, 'PYQT_VERSION'):
    QtGui.QApplication.instance().exec_()
