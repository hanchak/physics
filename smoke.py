'''
particle system to emulate rising smoke

Mike Hanchak
12OCT17
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

############# SETUP ######################
# bounds of the room
xlim = (-1.5,1.5)
ylim = (-0.,3.)

# create figure
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111, autoscale_on=False, xlim=xlim, ylim=ylim, facecolor='k')
ax.set_position([0,0,1,1])

########### DEFINE CLASS ######################
# particle class
class particle():
    '''simple particle class'''
    def __init__(self, x, y, vx, vy):
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.alpha = 1
        self.markersize = 10
        self.scatter, = ax.plot([], [], 'wo', ms=self.markersize, markeredgewidth=0)
        
    def update(self):
        self.x += self.vx
        self.y += self.vy
        self.alpha -= 0.015
        self.markersize += 0.7
        if self.alpha < 0:
            self.alpha = 0
        
    def show(self):
        self.scatter.set_data(self.x, self.y)
        self.scatter.set_alpha(self.alpha)
        self.scatter.set_markersize(self.markersize)


############## RUN ##########################
# run
Particles = []
def run(t):
    # make a new particle for every frame
    vx = np.random.randn(1)/400
    vy = 0.04 + np.random.randn(1)/300
    Particles.append(particle(0., 0., vx , vy))
    
    # loop through all particles to update and show them
    # delete ones for which the alpha is 0
    for i,item in enumerate(Particles):
       item.update()
       if item.alpha == 0:
           Particles.pop(i)
           item.scatter.remove()
       item.show()
    
    # return a list of iterables for the animation routine
    return [item.scatter for item in Particles]

################## ANIMATE ########################
# call the animation:
ani = animation.FuncAnimation(fig, run, range(1000), interval=0, blit=True)
plt.show()