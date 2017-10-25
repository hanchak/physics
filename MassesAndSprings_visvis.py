'''
Bouncing Ball code using object oriented programming. (3D Version!)

This code uses numpy arrays for all vector math

Code originally from here:
https://github.com/cpbotha/bwtl-python-tutorials/blob/master/
  part%205%20-%20object%20oriented%20programming%20and%20bouncing%20balls.ipynb
'''
#%% IMPORTS
import numpy as np
from random import random
import visvis as vv

# gravitational acceleration on Earth in m*s^-2
g = 0  #9.81
Ks = 50.  # spring constant
Kd = 0.3   # dashpot constant (viscous)
Kv = 0.2  # ball damping (viscous)
eqdist = 0.25  # equilibrium distance of spring
m = 1  # mass of ball

Nballs = 10

# acceleration vector due to g
ag = np.array([0., 0., -g]).reshape((1,3))
# coefficient of restitution (ratio of velocity after and before bounce)
cor = 0.55

# bounds of the room
xlim = [-0.5, 0.5]
ylim = xlim
zlim = xlim  #[0., 2.]

#  delta t
delta_t = 0.003

# create figure

axes = vv.gca()
axes.bgcolor = 'k'
axes.axis.visible=False


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
        self.object = vv.solidSphere(vv.Point(self.xy), scaling=(0.05, 0.05, 0.05))
        self.object.faceColor = [random(), random(), random()]
        self.trans = self.object.transformations[0]


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
        #self.xy[0,0] = np.clip(self.xy[0,0], xlim[0], xlim[1])
        #self.xy[0,1] = np.clip(self.xy[0,1], ylim[0], ylim[1])
        #self.xy[0,2] = np.clip(self.xy[0,2], zlim[0], zlim[1])

        # update the scatter plot with the new x and y positions:
        self.trans.dx = self.xy[0,0]
        self.trans.dy = self.xy[0,1]
        self.trans.dz = self.xy[0,2]



#%% CONSTRUCT OBJECTS
balls = []
for i in range(Nballs):
    startxy = np.array([random()*(xlim[1]-xlim[0])+ xlim[0], \
                      random()*(ylim[1]-ylim[0])+ ylim[0], \
                      random()*(zlim[1]-zlim[0])+ zlim[0]])
    
    balls.append( Ball( startxy, np.zeros(3) ) )

    axes.SetLimits(rangeX=xlim, rangeY=ylim, rangeZ=zlim)
    axes.camera.SetViewParams({'zoom':1})
    axes.camera.fov = 10
    

#%% ANIMATION ROUTINES
def animate(event):
    for ball in balls:  # loop through every ball
        ball.update(balls)  # call each ball's update method
    axes.Draw()


#%% CALL ANIMATION
timer = vv.Timer(axes, 0, False) 
timer.Bind(animate) 
timer.Start() 
app = vv.use() 
app.Run() 



