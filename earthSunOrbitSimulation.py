###############
# Simple orbit simulator: numnerical integration method
#
# Mike Hanchak, 19MAR2018
#
###############
import numpy as np
import matplotlib.pyplot as plt
import datetime
from scipy.optimize import fsolve

deg2rad = np.pi / 180
rad2deg = 1 / deg2rad

'''
From Wikipedia, Earth's orbit:
epoch J2000.0[nb 3] 
aphelion 152.10×106 km (94.51×106 mi)
 1.0167 AU[nb 4] 
perihelion 147.10×106 km (91.40×106 mi)
0.98329 AU[nb 4] 
semimajor axis 149.60×106 km (92.96×106 mi)
 1.000001018 AU[11] 
eccentricity 0.0167086[11] 
inclination 7.155° to Sun's equator
 1.578690°[12] to invariable plane 
longitude of the ascending node 174.9°[11] 
longitude of perihelion 102.9°[11] 
argument of periapsis 288.1°[11][nb 5] 
period 365.256363004 days[13] 
average speed 29.78 km/s (18.50 mi/s)[3]
 107,200 km/h (66,600 mph) 
 
angle of perihelion relative to vernal equinox = 102.9 deg

'''
# parameters
# Earth Sun system
mu = 1.32712440018e20  # m3/s2
au = 149597870700  # m/au

perihelion = 0.9832899 * au  # m
a = 149.60e9  # m
ecc = 0.0167086
speed_perihelion = np.sqrt( mu * (2/perihelion - 1/a))

earth_tilt = 23.43693 * deg2rad
long_peri = 102.94719 * deg2rad  # rad
earth_ang_vel = 7.292115090e-5  # rad/s
earth_period = 2*np.pi*np.sqrt(a**3/mu) / 3600  # hr

# time between perihelion and vernal equinox
time_perihelion = datetime.datetime(2018, 1, 3, 5, 35, tzinfo=datetime.timezone.utc)
time_vernal = datetime.datetime(2018, 3, 20, 16, 15, tzinfo=datetime.timezone.utc)
delta = time_vernal - time_perihelion
time_perihelion_vernal = delta.total_seconds() / 3600  # hr

# Dayton, OH
lat, long = 39.7589 * deg2rad, -84.1916 * deg2rad

def EarthPosAtTime(t):
    # position of earth relative to perihelion by direct numerical integration    
    r_p = np.array([0.9832899 ,0])  # au
    v_p = np.array([0, speed_perihelion ])  # m/s
    
    dt = 1  # hr
    
    # initialize vectors
    r = r_p.copy()
    v = v_p.copy()
    
    #N = np.round(t/dt)
    T = 0
    #for i in range(N):
    while T <= t:

        # SI units
        a = -mu* r / (np.sqrt(r[0]**2 + r[1]**2))**3 / au**2  # m/s/s
        
        v = v + a*dt*3600
        
        r = r + v*dt*3600/au
        
        T += dt
        
    return r

def EarthPosAtTime2(t):
    # position of earth relative to perihelion by Kepler's Method
    # https://en.wikipedia.org/wiki/True_anomaly
    
    mean_anomaly = 2 * np.pi * t / earth_period
    eccentric_anomaly = fsolve(lambda E:E - ecc*np.sin(E) - mean_anomaly, 0)[0]
    
    #true_anomaly = np.arccos((np.cos(eccentric_anomaly) - ecc)/(1 - ecc*np.cos(eccentric_anomaly)))
    true_anomaly = 2*np.arctan2(np.sqrt(1+ecc)*np.sin(eccentric_anomaly/2),np.sqrt(1-ecc)*np.cos(eccentric_anomaly/2))
    
    r_mag = a*(1 - ecc**2) / (1 + np.cos(true_anomaly)) / au
    
    return np.array( ( r_mag*np.cos(true_anomaly) , r_mag*np.sin(true_anomaly) ) )

def rotZ(ang):
    # return the rotation matrix (ang in radians)
    return np.array([[np.cos(ang),-np.sin(ang),0,0],[np.sin(ang),np.cos(ang),0,0],[0,0,1,0],[0,0,0,1]])

def rotX(ang):
    # return the rotation matrix (ang in radians)
    return np.array([[1,0,0,0],[0, np.cos(ang),-np.sin(ang),0],[0, np.sin(ang),np.cos(ang),0],[0,0,0,1]])
    
def rotY(ang):
    # return the rotation matrix (ang in radians)
    return np.array([[np.cos(ang),0,np.sin(ang),0],[0,1,0,0],[-np.sin(ang),0,np.cos(ang),0],[0,0,0,1]])
    
def trans(r):
    if len(r) == 2:
        z = 0
    else:
        z = r[2]
    # return the translation matrix
    return np.array([[1,0,0,r[0]],[0,1,0,r[1]],[0,0,1,z],[0,0,0,1]])

# sun frame to local frame
def sunPosition(given_time):
    # take a given time and output location of sun in local frame (alt and azi)

    #for given_time in range(12):
    #for given_time in np.arange(5.5, 9000, 24*10):

    #given_time = 0  # hours since vernal equinox UTC

    r = EarthPosAtTime2(given_time + time_perihelion_vernal)
    
    greenwich_ang = given_time * earth_ang_vel*3600
    
    mat = rotZ(long_peri)@trans(r)@rotZ(-long_peri)@rotX(earth_tilt)@rotZ(greenwich_ang+long)@rotY(-lat)
    
    # get the sun's position in local frame
    sun_fixed = np.array([[0],[0],[0],[1]])
    sun_local = np.linalg.inv(mat) @ sun_fixed
    
    sun_alt = np.arctan2(sun_local[0], np.sqrt(sun_local[1]**2 + sun_local[2]**2))*rad2deg
    sun_azi = 90 + np.arctan2(sun_local[2], sun_local[1])*rad2deg

    return sun_alt, sun_azi

for i in range(0,365,75):
    alts = []
    azis = []
    for given_time in np.arange(-6+i*24, 18+i*24, 0.25):
    #for given_time in np.arange(5.62, 9000, 24):    
        sun_alt, sun_azi = sunPosition(given_time)
        alts.append(sun_alt)
        azis.append(sun_azi)
        
    plt.plot(azis, alts)
    
#print(sun_local, sun_alt, sun_azi)

#plt.plot(r_out[0,:], r_out[1,:])
#plt.plot(0,0,'ro')
plt.show()