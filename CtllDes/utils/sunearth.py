from poliastro.bodies import Earth, Sun
from poliastro.frames import Planes
from poliastro.constants import J2000
from poliastro.ephem import Ephem

import astropy.units as u
from astropy.time import Time

import numpy as np

from ..core import ctll,satellite


def bodies_vector(T, dt, body1, body2, epoch=J2000, timedelta = None):
    """Returns vector for every dt in [0,T] from bodie1 to bodie2
    
    Parameters
    ----------
    T : float
        time of propagation
    dt : float
        time interval
    body1 : ~poliastro.bodies.Body
        first body
    body2 : ~poliastro.bodies.Body
        second body
    epoch : ~astropy.time.Time, optional
        Epoch offset
    timedelta : ~astropy.time.Time, optional
        desired time of propagation
    
    Returns
    -------
    diff : ~astropy.units.quantity.Quantity
        vectors for timedelta propagation from body1 to body2
    """
    
    DT_TIME = 500
    if not timedelta:
        time_delta = Time([epoch + j*u.s for j in np.linspace(0,3600*24*1.01*T, int(3600*24*T/DT_TIME)) ] )
    else:
        time_delta = timedelta
    
    ephem_b1 = Ephem.from_body(body1, time_delta.tdb)
    ephem_b2 = Ephem.from_body(body2, time_delta.tdb)
    
    tofs = Time([epoch + j*dt*u.s for j in range(10,int(3600*24*T/dt)+10)])
    r_b1, _ = ephem_b1.rv(tofs)
    r_b2, _ = ephem_b2.rv(tofs)
    
    diff = (r_b2-r_b1).to(u.km)
    return diff