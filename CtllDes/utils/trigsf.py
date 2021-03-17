import math
import numpy as np 
from .. import units as u

def SSA(a,b,A):
    #returns c, B,C where capital letters indicate angles, and small ones great circle arcs
    B = math.asin((math.sin(b)*math.sin(A))/math.sin(a))%(2*u.pi)
    c = MSL(a,b,A,B)
    #C = MAL(A,B,a,b) not implemented because i dont need it (yet)
    
    return c#, B #,C

def AAS(A,B,a):
    b = math.asin((math.sin(B)*math.sin(a))/math.sin(A))
    c = MSL(a,b,A,B)
    #C = MAL(A,B,a,b) not implemenet because i dont need it (yet)
    return b, c#,C

def MSL(a,b,A,B):
    sinc_n = math.sin(a)*math.cos(b)*math.cos(B)+math.sin(b)*math.cos(a)*math.cos(A)
    sinc_d = 1 - math.sin(a)*math.sin(b)*math.sin(A)*math.sin(B)
    sinc = sinc_n/sinc_d
    
    cosc_n = math.cos(a)*math.cos(b)-math.sin(a)*math.sin(b)*math.cos(A)*math.cos(B)
    cosc_d = 1 - math.sin(a)*math.sin(b)*math.sin(A)*math.sin(B)
    cosc = cosc_n/cosc_d
    
    return math.atan2(sinc,cosc)

def pitolat(rad):
    """Util solamente por ahora para traducir resultados de AAS a lat y long"""
    sgn = np.sign(rad)
    if sgn == 1:
        if  rad > u.pi/2:
            return u.pi/2 - rad%(u.pi/2)
        else:
            return rad
    else:
        if rad < -u.pi/2:
            return -u.pi/2 + rad%(u.pi/2)
        else: 
            return rad

def c2s(x,y,z):
    """Transformacion de coordenadas cartesianas a esfericas.
    
    Recibe:
        x,y,z
    Devuelve: r, lat, lon
    """
    r = (x**2+y**2+z**2)**(1/2)
    lat = np.arccos( (x**2+y**2)**(1/2) /r)
    if z<0:
        lat *= -1
    lon = np.arctan2(y,x)

    return r,lat,lon

def rad2deg(rad):
    return rad*u.rad2deg

def deg2rad(deg):
    return deg/u.rad2deg
