'''
Bouncing Ball code using object oriented programming.

This code uses numpy arrays for all vector math

Code originally from here:
https://github.com/cpbotha/bwtl-python-tutorials/blob/master/
  part%205%20-%20object%20oriented%20programming%20and%20bouncing%20balls.ipynb
'''
#%% IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# gravitational acceleration on Earth in m*s^-2
g = 0  #9.81
Ks = 50.  # spring constant
Kd = 0.2   # dashpot constant (viscous)
eqdist = 0.5  # equilibrium distance of spring
m = 1  # mass of ball

Nballs = 15

# acceleration vector due to g
ag = np.array([0.,-g])
# coefficient of restitution (ratio of velocity after and before bounce)
# see http://en.wikipedia.org/wiki/Coefficient_of_restitution
cor = 0.50

# bounds of the room
xlim = (0.,2.)
ylim = (0.,2.)

#  delta t
delta_t = 0.005

# create figure
fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=xlim, ylim=ylim)
ax.grid('off')


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
        self.scatter, = ax.plot([], [], 'o', markersize=20)

    def update(self, others):
        '''update the position from equations of motion & update plot data'''
        if self.xy[0] <= xlim[0] or self.xy[0] >= xlim[1]:
            # hit the left or right wall, reflect x component
            self.v[0] *= - cor

        if self.xy[1] <= ylim[0] or self.xy[1] >= ylim[1]:
            # hit the bottom or top wall, reflect y component
            self.v[1] *= - cor
        
        # loop through all other balls to calcuate attractive force
        force = np.zeros(2)
        
        for other in others:
            # skip self
            if self.scatter != other.scatter:  # dont check self
                vec = other.xy - self.xy
                
                dist = np.sqrt(vec[0]**2 + vec[1]**2)
                
                unitvec = vec / dist
                
                vec_rel = (other.v - self.v).dot(unitvec)
                
                if dist > 0.0001:
                    # spring force between balls
                    force +=  (dist - eqdist) * Ks * unitvec
                    # damping force between balls
                    force +=  vec_rel * Kd * unitvec
                
        # add damping force to the ball
        dforce = - 0.5 * self.v
        
        # Euler integration: v(i+1) = v(i) + acc*dt
        self.v += (dforce/m + force/m + ag) * delta_t

        # Euler integration: x(i+1) = x(i) + vel*dt
        self.xy += self.v * delta_t

        # keep xs and ys inside box:
        self.xy[0] = np.clip(self.xy[0], xlim[0], xlim[1])
        self.xy[1] = np.clip(self.xy[1], ylim[0], ylim[1])

        # update the scatter plot with the new x and y positions:
        self.scatter.set_data(self.xy)


#%% CONSTRUCT OBJECTS


balls = [Ball( (np.random.random(2)) * ylim[1], np.random.randn(2)*0.0 ) \
         for i in range(Nballs)]

#%% ANIMATION ROUTINES
def init():
    return []

def animate(t):
    # t is time in seconds
    for ball in balls:
        ball.update(balls)

    # have to return an iterable
    return [ball.scatter for ball in balls]

#%% CALL ANIMATION

# interval in milliseconds
# we're watching in slow motion (delta t is shorter than interval)
ani = animation.FuncAnimation(fig, animate, np.arange(0,1,delta_t) \
                              , init_func=init, interval=0.1, blit=True)

plt.show()
