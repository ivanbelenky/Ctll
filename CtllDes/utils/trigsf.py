import math
import numpy as np 


import astropy.units as u


def SSA(a,b,A):
    #returns c, B,C where capital letters indicate angles, and small ones great circle arcs
    B = math.asin((math.sin(b)*math.sin(A))/math.sin(a))%(2*np.pi)
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
        if  rad > np.pi/2:
            return np.pi/2 - rad%(np.pi/2)
        else:
            return rad
    else:
        if rad < -np.pi/2:
            return -np.pi/2 + rad%(np.pi/2)
        else: 
            return rad

def c2s(x,y,z):
    """Coordinates transformation. Cartesian to Spherical.
    
    Parameters
    ----------
    x : ~astropy.units.quantity.Quantity
        x coordinate
    y : ~astropy.units.quantity.Quantity
        y coordinate
    z : ~astropy.units.quantity.Quantity
        z coordinate

    Returns
    -------
    r : ~astropy.units.quantity.Quantity
        radii 
    lat : ~numpy.ndarray
        latitudes            
    lon : ~numpy.ndarray
        longitudes 
    """
    
    r = (x**2+y**2+z**2)**(1/2)
    lat = np.sign(z)*np.arccos( ((x**2+y**2)**(1/2)) /r)
    lon = np.arctan2(y,x)

    return r,lat,lon


def get_lam(r,FOV,R):
    """Maximum attractor's central angle.
    
    Parameters
    ----------
    r : ~astropy.units.quantity.Quantity
        positions in kilometers or any distance metric
    FOV : ~astropy.units.quantity.Quantity
        field of fiew in radians
    R : ~astropy.units.quantity.Quantity
        attractor mean radius

    Returns
    -------
    lams :  ~astropy.units.quantity.Quantity
        maximum(s) attractor's central angle(s).

    """

    radiis = np.sqrt(np.sum(r**2,axis=1))
    rho = np.arcsin(R/radiis)
    eps = np.arccos((np.sin(FOV))/(np.sin(rho)))

    lams = (np.pi/2)*u.rad - FOV - eps

    return lams

def get_angles(lons,lats,t_lon,t_lat):
    """Angles between two points, first and second,
     in the surface of a sphere.

    Parameters
    ----------
    lons : ~astropy.units.quantity.Quantity
        first longitudes in radians
    lats : ~astropy.units.quantity.Quantity 
        first latitudes in radians
    t_lons : ~astropy.units.quantity.Quantity
        second longitudes in radians
    t_lats : ~astropy.units.quantity.Quantity 
        second latitudes in radians

    Returns
    -------
    angles : ~numpy.ndarray
        angle between first and second points

    """
    
    s_lat_tgt = np.sin(t_lat)
    c_lat_tgt = np.cos(t_lat)
    s_lat_ssps = np.sin(lats)
    c_lat_ssps = np.cos(lats)
    c_lon_r = np.cos(t_lon-lons)
    angles = np.arccos(s_lat_tgt*s_lat_ssps+
        c_lat_ssps*c_lat_tgt*c_lon_r)    

    return angles