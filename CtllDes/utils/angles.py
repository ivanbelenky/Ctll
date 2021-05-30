import numpy as np

import astropy.units as u

from poliastro.constants import J2000
from poliastro.bodies import Sun

from . import sunearth

def get_passage_raan(sat, target_lon, target_lat, lon_offset = 0):
	"""Get right ascension of the ascending node in
	order to fly over target at the closest 
	achievable distance.

	Parameters
	----------
	sat : ~CtllDes.core.satellite.Sat
	    Satellite with desired orbit
	target_lon : ~astropy.units.quantity.Quantity
	    target longitude in degrees
	target_lat : ~astropy.units.quantity.Quantity
	    target latitude in degrees
	    
	Returns
	-------
	raan : ~astropy.units.quantity.Quantity
	    right ascension of the ascending node
	    such that, closest passage to target
	    is achieved by the satellite.
	"""

	t_lon = target_lon.to(u.rad)
	t_lat = target_lat.to(u.rad)

	t_lon = t_lon+np.pi*u.rad


	period = sat.orbit.period
	lon,lat = sat.ssps(period.to(u.d).value, lon_offset=lon_offset, dt=1, J2=True)
	raan = t_lon.value-lon[np.abs(lat.value-t_lat.value).argmin()].value
	
	if np.abs(raan) > np.pi:
	    raan = -1*np.sign(raan)*(2*np.pi-np.abs(raan))
	return raan*u.rad


def get_passage_LTDN(sat, epoch=J2000, LTDN=12*u.h):
	"""Get longitude offset for propagation. If
	right ascension of the ascending node correction
	has already been calculated using get_passage_raan,
	in order to preserve passage over target, lon_offset
	must be added to raan.

	Parameters
	----------
	sat : ~CtllDes.core.satellite.Sat
	    Satellite with desired orbit
	epoch : ~astropy.time.Time, optional
	    Epoch, default to J2000. 
	LTDN : ~astropy.time.Time, optional
	    Desired Local Time of Descending node

	Returns
	-------
	lon_offset : ~astropy.units.quantity.Quantity
	    longitude offset
	"""

	t_ltdn = (180*u.deg-sat.orbit.nu.to(u.deg))/(360*u.deg)*sat.orbit.period
	t_ltdn = t_ltdn.to(u.d).value

	attr = sat.attractor

	_dt = t_ltdn*3600*24/1000 

	r,v = sat.rv(t_ltdn, dt=_dt, J2=True)
	r_tldn = r[-1]

	diff = sunearth.bodies_vector(t_ltdn, _dt, attr, Sun, epoch=epoch, timedelta=None)
	diff_ltdn = diff[-1]

	diff_x = diff_ltdn[0]
	diff_y = diff_ltdn[1]

	r_tldn_x = r_tldn[0]
	r_tldn_y = r_tldn[1]

	theta_r = np.arctan2(r_tldn_y,r_tldn_x)
	theta_diff = np.arctan2(diff_y,diff_x)

	delta_theta = (LTDN/(12*u.h)-1)*np.pi*u.rad

	t_sgn = np.sign(theta_r-theta_diff)
	t_abs = np.abs(theta_r-theta_diff)

	if t_abs > np.pi*u.rad:
	    t_ = t_sgn*(2*np.pi*u.rad - t_abs)
	else:
	    t_ = theta_r - theta_diff

	lon_offset = t_ + delta_theta
	return lon_offset
    
    


