################################################################################
# Simple Earth orbit simulator: Kepler's method with coordinate transformations
#
# Uses Vernal equinox (UTC) of 2018 as reference time.  However, this does not
#   coorespond to high no0n at Greenwich... there is an offset +4.25 hours.
#
# Uses standard time... convert to DST at own risk.
#
# Mike Hanchak, 19MAR2018
#
################################################################################
import numpy as np
import matplotlib.pyplot as plt
import datetime, pytz
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
earth_ang_vel = 2*np.pi/86164.091   #7.2921159e-5  # rad/s
earth_period = 2*np.pi*np.sqrt(a**3/mu) / 3600  # hr

# time between perihelion and vernal equinox
# TIME AT VERNAL EQUINOX IS EPOCH FOR THIS CODE!
time_perihelion = datetime.datetime(2018, 1, 3, 5, 35, tzinfo=datetime.timezone.utc)
time_vernal = datetime.datetime(2018, 3, 20, 16, 15, tzinfo=datetime.timezone.utc)
vernal_offset = time_vernal - datetime.datetime(2018, 3, 20, 12, tzinfo=datetime.timezone.utc)
delta = time_vernal - time_perihelion
time_perihelion_vernal = delta.total_seconds() / 3600  # hr

# Dayton, OH
lat, long = 39.7589 * deg2rad, -84.1916 * deg2rad
local_tz = pytz.timezone('US/Eastern')
# # topocentric
# lat, long = 0,0
# local_tz = pytz.timezone('UTC')
# # Greenwich
# lat, long = 51.476852* deg2rad, -0.000500* deg2rad
# local_tz = pytz.timezone('UTC')


def EarthPosAtTime(t):
    # position of earth relative to perihelion by Kepler's Method
    # https://en.wikipedia.org/wiki/True_anomaly
    
    mean_anomaly = 2 * np.pi * t / earth_period
    eccentric_anomaly = fsolve(lambda E:E - ecc*np.sin(E) - mean_anomaly, 0)[0]
    
    #true_anomaly = np.arccos((np.cos(eccentric_anomaly) - ecc)/(1 - ecc*np.cos(eccentric_anomaly)))
    num = np.sqrt(1+ecc)*np.sin(eccentric_anomaly/2)
    den = np.sqrt(1-ecc)*np.cos(eccentric_anomaly/2)
    true_anomaly = 2*np.arctan2(num, den)
    
    r_mag = a*(1 - ecc**2) / (1 + ecc*np.cos(true_anomaly)) / au
    
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
    # take a given time relative to vernal equinox and output location of sun in local frame (alt and azi)

    r = EarthPosAtTime(given_time + time_perihelion_vernal)
    
    greenwich_ang = vernal_offset.total_seconds()/3600*15*deg2rad  # angle of greenwich meridian from vernal equinox
    time_ang = given_time * earth_ang_vel*3600
    
    mat = rotZ(long_peri)@trans(r)@rotZ(-long_peri)@rotX(-earth_tilt)@rotZ(greenwich_ang+long+time_ang)@rotY(-lat)
    
    # get the sun's position in local frame
    sun_fixed = np.array([[0],[0],[0],[1]])
    sun_local = np.linalg.inv(mat) @ sun_fixed
    
    sun_alt = np.arctan2(sun_local[0], np.sqrt(sun_local[1]**2 + sun_local[2]**2))*rad2deg
    
    sun_azi = 90 + np.arctan2(sun_local[2], sun_local[1])*rad2deg
    if sun_azi > 180:
        sun_azi -= 360  # adjust angle extents for prettier plot
    
    return sun_alt, sun_azi, sun_local
   
suns = []
for day in range(0,365,365):
    alts = []
    azis = []
    
    for hour in np.arange(-6,6,0.25):  #np.arange(-6+i*24, 17+i*24, 0.25):
    #for given_time in np.arange(5.62, 9000, 24):
        given_time = datetime.datetime(2018, 3, 20, 12, 0, tzinfo=local_tz) + datetime.timedelta(hours=day*24+hour)
        offset = given_time - time_vernal
        sun_alt, sun_azi, sun_local = sunPosition(offset.total_seconds()/3600)
        alts.append(sun_alt)
        azis.append(sun_azi)
        suns.append(sun_local[:])
    
        plt.plot(sun_azi,sun_alt,'bo')
        plt.text(sun_azi,sun_alt,given_time.strftime('%d%b%y %H:%M'))
    
    plt.plot(azis, alts,'b-')    
    
plt.grid(True)
#given_time = datetime.datetime.now(tz=pytz.timezone('EST'))
# given_time = datetime.datetime(2018, 3, 20, 12, 36, tzinfo=pytz.timezone('US/Eastern'))
# offset = given_time - time_vernal + vernal_offset
# print(sunPosition(offset.total_seconds()/3600))

############## 3D plot of sun in local frame #################    
# from mpl_toolkits.mplot3d import Axes3D
# suns = np.concatenate(suns, axis=1)   
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot(suns[0,:],suns[1,:],suns[2,:])


plt.grid(True)
plt.show()
