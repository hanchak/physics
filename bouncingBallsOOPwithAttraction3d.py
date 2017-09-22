'''
Bouncing Ball code using object oriented programming. (3D Version!)

This code uses numpy arrays for all vector math

Code originally from here:
https://github.com/cpbotha/bwtl-python-tutorials/blob/master/
  part%205%20-%20object%20oriented%20programming%20and%20bouncing%20balls.ipynb
'''
#%% IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  
#import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

# gravitational acceleration on Earth in m*s^-2
g = 0  #9.81
Ks = 50.  # spring constant
Kd = 0.5   # dashpot constant (viscous)
Kv = 0.5  # ball damping (viscous)
eqdist = 0.5  # equilibrium distance of spring
m = 1  # mass of ball

Nballs = 6

# acceleration vector due to g
ag = np.array([0., 0., -g])
# coefficient of restitution (ratio of velocity after and before bounce)
# see http://en.wikipedia.org/wiki/Coefficient_of_restitution
cor = 0.50

# bounds of the room
xlim = (0.,2.)
ylim = xlim
zlim = xlim

#  delta t
delta_t = 0.005

# create figure
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')  #ax = p3.Axes3D(fig)
ax.set_xlim3d(xlim)
ax.set_ylim3d(ylim)
ax.set_zlim3d(zlim)
#plt.axis('off')
#ax.grid('off')


#%% DEFINE CLASS
class Ball():
    '''This class creates a ball and provides for physics of motion'''
    def __init__(self, xy, v):
        """
        :param xy: Initial position.
        :param v: Initial velocity.
        """
        # cast inputs as numpy arrays
        self.xy = np.array(xy)
        self.v = np.array(v)

        # set up a plot object, which we will reference later.
        self.scatter, = ax.plot([], [], [], 'o', markersize=20)

    def update(self, others):
        '''update the position from equations of motion & update plot data'''
        if self.xy[0] <= xlim[0] or self.xy[0] >= xlim[1]:
            # hit the left or right wall, reflect x component
            self.v[0] *= - cor

        if self.xy[1] <= ylim[0] or self.xy[1] >= ylim[1]:
            # hit the bottom or top wall, reflect y component
            self.v[1] *= - cor
            
        if self.xy[2] <= zlim[0] or self.xy[2] >= zlim[1]:
            # hit the bottom or top wall, reflect y component
            self.v[2] *= - cor
        
        # loop through all other balls to calcuate attractive force
        force = np.zeros(3)
        
        for other in others:
            # skip self
            if self.scatter != other.scatter:  # dont check self
                vec = other.xy - self.xy  # position vec of other w.r.t. selfff
                
                dist = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)  # scalar distance
                
                unitvec = vec / dist  # unit vector from self to other
                
                vec_rel = (other.v - self.v).dot(unitvec)  #rel vel along unit vec
                
                if dist > 0.0001:
                    # spring force between balls
                    force +=  (dist - eqdist) * Ks * unitvec
                    # damping force between balls
                    force +=  vec_rel * Kd * unitvec
                
        # add damping force to the ball proportional to absolute velocity
        dforce = - Kv * self.v
        
        # Euler integration: v(i+1) = v(i) + acc*dt
        self.v += (dforce/m + force/m + ag) * delta_t

        # Euler integration: x(i+1) = x(i) + vel*dt
        self.xy += self.v * delta_t

        # keep xs and ys inside box:
        self.xy[0] = np.clip(self.xy[0], xlim[0], xlim[1])
        self.xy[1] = np.clip(self.xy[1], ylim[0], ylim[1])
        self.xy[2] = np.clip(self.xy[2], zlim[0], zlim[1])

        # update the scatter plot with the new x and y positions:
        self.scatter.set_data(self.xy[0], self.xy[1])
        self.scatter.set_3d_properties(self.xy[2])
        


#%% CONSTRUCT OBJECTS


balls = [Ball( (np.random.random(3)) * ylim[1], np.random.randn(3)*0.0 ) \
         for i in range(Nballs)]
#balls2 = [Ball( (np.random.random(2)) * ylim[1], np.random.randn(2)*0.0 ) \
#         for i in range(Nballs)]

#%% ANIMATION ROUTINES
def init():
    return []

def animate(t):
    # t is time in seconds
    for ball in balls:
        ball.update(balls)
        
    #for ball in balls2:
    #    ball.update(balls2)
    #ax.view_init(elev=10., azim=t*100)
    # have to return an iterable
    #return [ball.scatter for ball in balls] #+balls2]
    return [ball.scatter for ball in balls]

#%% CALL ANIMATION

# interval in milliseconds
# we're watching in slow motion (delta t is shorter than interval)
ani = animation.FuncAnimation(fig, animate, np.arange(0,1,delta_t) \
                              , init_func=init, interval=1, blit=0)

plt.show()
