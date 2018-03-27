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
import os

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

earth_tilt = 23.43693 * deg2rad  # rad, earth's tilt
long_peri = 102.94719 * deg2rad  # rad, longitude of the earth's perihelion
earth_ang_vel = 2*np.pi/86164.091   #7.2921159e-5  # rad/s
earth_period = 2*np.pi*np.sqrt(a**3/mu) / 3600  # hr

# time between perihelion and vernal equinox
# UTC TIME AT VERNAL EQUINOX IS EPOCH FOR THIS CODE!
time_perihelion = datetime.datetime(2018, 1, 3, 5, 35, tzinfo=datetime.timezone.utc)
time_vernal = datetime.datetime(2018, 3, 20, 15, 0, tzinfo=datetime.timezone.utc)
vernal_offset = time_vernal - datetime.datetime(2018, 3, 20, 12, tzinfo=datetime.timezone.utc)
greenwich_ang = vernal_offset.total_seconds()/3600*15*deg2rad  # rad

#delta = time_vernal - time_perihelion
#time_perihelion_vernal = delta.total_seconds() / 3600  # hr

# Dayton, OH
lat, long = 39.7589 * deg2rad, -84.1916 * deg2rad
local_tz = pytz.timezone('US/Eastern')
# # topocentric
#lat, long = 0,0
#local_tz = pytz.timezone('UTC')
# # Greenwich
# lat, long = 51.476852* deg2rad, -0.000500* deg2rad
# local_tz = pytz.timezone('UTC')



def earthPosInFixed(given_time):
    # give the vector of the earth in heliocentric fixed frame, Kepler's Method
    
    # calculate time since perihelion (reference for Kepler's Method)
    t = (given_time - time_perihelion).total_seconds()/3600  # hr
    
    # mean anomaly = proportion of circular orbit area covered in the time
    mean_anomaly = 2 * np.pi * t / earth_period
    
    # Eccentric anomaly: Kepler's equation
    eccentric_anomaly = fsolve(lambda E:E - ecc*np.sin(E) - mean_anomaly, 0)[0]
    
    # true_anomaly: elliptical orbit angle of satellite
    num = np.sqrt(1+ecc)*np.sin(eccentric_anomaly/2)
    den = np.sqrt(1-ecc)*np.cos(eccentric_anomaly/2)
    true_anomaly = 2*np.arctan2(num, den)
    
    # elliptic orbit: distance to satellite
    r_mag = a*(1 - ecc**2) / (1 + ecc*np.cos(true_anomaly)) / au
    
    # return the x and y positions of satellite
    return rotZ(long_peri + true_anomaly) @ np.array([[r_mag],[0],[0],[1]])
    
def earthRotAng(given_time):
    # returns the angle of Greenwich meridian from the vernal equinox at the given time
    
    # how much time has passed since vernal equinox?
    t = given_time - time_vernal
    # angle is time * spin velocity
    time_ang = t.total_seconds() * earth_ang_vel
    
    return time_ang + greenwich_ang

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

def celestial2fixed(ra, dec):
    # function to convert rigt acension and declination data on the celestial
    # sphere to x,y,z points in the fixed, heliocentric frame
    xyz = np.array([[np.cos(dec*deg2rad)*np.cos(ra*deg2rad)],[np.cos(dec*deg2rad)*np.sin(ra*deg2rad)],[np.sin(dec*deg2rad)],[1]])
    return rotX(-earth_tilt) @ xyz

def fixed2local(r_obj, given_time, long, lat, trans_flag):
    # coordinate transform from fixed, heliocentric reference frame to local
    # frame at viewers lat and long
    # last argument = True: includes a translation to earth center, not needed
    #     for star positions, only for planets, sun
    r_earth = earthPosInFixed(given_time)
    ang_earth = earthRotAng(given_time)
    
    # transformation matrix:
    if trans_flag:
        mat = trans(r_earth) @ rotX(-earth_tilt) @ rotZ(ang_earth + long)@rotY(-lat)
    else:
        mat = rotX(-earth_tilt) @ rotZ(ang_earth + long)@rotY(-lat)
    return np.linalg.inv(mat) @ r_obj

def xyz2AltAzi(r):
    # given x,y,z in the local frame, return altitude and azimuth
    alt = np.arctan2(r[0], np.sqrt(r[1]**2 + r[2]**2) ) * rad2deg
    
    azi = 90 + np.arctan2(r[2], r[1]) * rad2deg
    if azi > 180:
        azi -= 360  # adjust angle extents for prettier plot
    
    return alt, azi

def altAzi2SkyView(alt, azi):
    # convert alt and azi (in degrees) to x,y positions on a circular sky view map
    r = (90 - alt) / 90
    th = 270 - azi
    x = r*np.cos(th*deg2rad)
    y = r*np.sin(th*deg2rad)
    
    return x,y

# here is a loop to plot suns position in local long and lat   
# suns = []
# for day in range(0,365,365):
#     alts = []
#     azis = []
#     
#     for hour in np.arange(-6,6,0.25):
# 
#         given_time = datetime.datetime(2018, 3, 23, 12, 0, tzinfo=local_tz) + datetime.timedelta(hours=day*24+hour)
#         
# 
#         sun_fixed = np.array([[0],[0],[0],[1]])
#         sun_alt, sun_azi = xyz2AltAzi( fixed2local(sun_fixed, given_time, long, lat) )
#         
#         alts.append(sun_alt)
#         azis.append(sun_azi)
# 
#         plt.plot(sun_azi,sun_alt,'bo')
#         plt.text(sun_azi,sun_alt,given_time.strftime('%d%b%y %H:%M'))
#     
#     plt.plot(azis, alts,'b-')    
#     
# plt.grid(True)

# plot a skymap
given_time = datetime.datetime(2018, 3, 27, 19, 0, tzinfo=local_tz)
#given_time = datetime.datetime.now(tz=local_tz)
plt.axis('equal')

sun_fixed = np.array([[0],[0],[0],[1]])

# first plot sun
plt.plot(*altAzi2SkyView( *xyz2AltAzi( fixed2local(sun_fixed, given_time, long, lat, True) ) ),'yo',markersize = 10)

# plot horizon, 45 alt, and zenith dot
th = np.linspace(0,2*np.pi,100)
plt.plot(np.cos(th), np.sin(th), 'g-')
plt.plot(0.5*np.cos(th), 0.5*np.sin(th), 'r:')
plt.plot(0,0,'r+')

# plot direction markers
plt.text(1,0,'West', color='w')
plt.text(-1,0,'East', color='w') 
plt.text(0,1,'North', color='w')
plt.text(0,-1,'South', color='w')

#plt.figure()

# plot stars
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.gca(projection='3d')
plt.gca().set_facecolor('black')

os.chdir(r'C:\Users\HanchaMS\Documents\Python Scripts\physics')
with open('brightstar_2018.txt') as file:
    for i in range(5):
        file.readline()
        
    for i in range(1469):
        s = file.readline()
        
        ratext = s[27:37].split()
        ra = (float(ratext[0]) +  float(ratext[1])/60 + float(ratext[2])/3600) * 15 
        
        dectext = s[41:49].split()
        dec = float(dectext[0]) +  float(dectext[1])/60 + float(dectext[2])/3600
        mag = float(s[59:64])
        if s[40] == '-':
            dec = -1*dec
        if mag > 4.5:
            #print(ra,dec)
            alt, azi = xyz2AltAzi( fixed2local( celestial2fixed(ra, dec) , given_time, long, lat, False))
            #r = fixed2local( celestial2fixed(ra, dec) , given_time, long, lat, False)
            #r = celestial2fixed(ra, dec)
            
            if alt > 0:
            #    plt.plot(r[1],r[2],'w.' , markersize = round(mag) )
            #    plt.text(-r[1], r[2] , s[16:19], fontsize=8, color='w')
                x,y = altAzi2SkyView(alt,azi )
                plt.plot(x,y ,'w.' , markersize = round(mag) )
                plt.text(x,y , s[16:19], fontsize=8, color='w')
        
        #plt.plot(azi,alt  ,'k.' )
        #plt.plot(ra,dec  ,'k.' )
        #ax.plot(r[0],r[1],r[2],'k.')
        
plt.show()
