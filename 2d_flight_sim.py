# 2-d flight simulator
# M. Hanchak, 07MAR18
# roughly a Cessna 172
# inputs: AOA and thrust
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

rad2deg = lambda ang:ang * 180 / np.pi

def C_L(aoa):
    # roughly NACA 2412 
    return 2.25/20 * rad2deg(aoa) + 0.25
    
def C_D(aoa):
    # roughly NACA 2412 
    return 0.02/100 * rad2deg(aoa)**2 + 0.02

# constants
g = 9.81  # m/s2
m = 1136  # kg
rho = 1.0  # kg/m3
A_L = 16  # m2
A_D = 1  # m2
dt = 0.50  # s
aoa = 0.0  # rad
F_T = 80  # N

# equations of motion
def accels(aoa, F_T, vx, vy):
    v_mag = np.sqrt(vx**2 + vy**2)
    theta = np.arctan2(vy,vx)
    F_D_mag = 0.5*C_D(aoa)*A_D*rho*v_mag**2
    F_L_mag = 0.5*C_L(aoa)*A_L*rho*v_mag**2
    F_G_mag = m*g
    ax = 1/m*(F_T*np.cos(theta + aoa) - F_L_mag*np.sin(theta) - F_D_mag*np.cos(theta) )
    ay = 1/m*(F_T*np.sin(theta + aoa) + F_L_mag*np.cos(theta) - F_D_mag*np.sin(theta) -m*g)
    return ax, ay
    
#%% set up plotting
fig, ax = plt.subplots()

pos_x = [0]
pos_y = [1000]
vel_x = [80]  # m/s
vel_y = [0]

#myline = plt.Line2D(pos_x, pos_y)
myline, = plt.plot(pos_x, pos_y, animated=True)                                                                                                                                                                                                                                                                                                                              

#ax.add_line(myline)



#%% ANIMATION ROUTINES
def init():
    ax.set_xlim(0,10000)
    ax.set_ylim(pos_y[0]-500, pos_y[0]+500)
    return myline,

def animate(t):
    # t is time in seconds

    ax, ay = accels(aoa, F_T, vel_x[-1], vel_y[-1])
    
    vel_x.append(vel_x[-1] + ax*dt)
    vel_y.append(vel_y[-1] + ay*dt)
    
    pos_x.append(pos_x[-1] + vel_x[-1]*dt)
    pos_y.append(pos_y[-1] + vel_y[-1]*dt)
    
    myline.set_data(pos_x, pos_y)
    
    #time.append(t)
    
    # have to return an iterable?
    return myline,

#%% CALL ANIMATION
# interval in milliseconds
# we're watching in slow motion (delta t is shorter than interval)
ani = animation.FuncAnimation(fig, animate, np.arange(0,100,dt) \
                              , init_func=init, interval=1, blit=True, repeat=False)

plt.show()
