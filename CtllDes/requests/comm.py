import astropy.units as u

import numpy as np
import pandas as pd

from . import coverage 

from ..utils import trigsf




def default_comm_data(sat, target, T, dt=1., **kwargs):
	"""Useful information for coverage data. That is
	positions, velocities, latitudes and longitudes, elevations,
	and times of view of the target. Just the covered data is displayed.

	Parameters
	----------
	sat : CtllDes.core.satellite.Sat
		satellite
	target : CtllDes.targets.targets.Target
		desired target
	T : int | float
		desired time of analysis
	dt : int | float
		desired time interval of sampling
	kwargs : dict
		keywrd arguments for the propagation module
	"""

	r,v = sat.rv(T,dt,**kwargs)
	lons,lats = sat.ssps(T,dt,**kwargs)

	view = coverage.get_view(lons, lats, r, target, sat.attractor.R_mean)
	
	lons_ = []
	lats_ = []
	r_ = []
	v_ = []
	times = []

	for i in range(len(view)):
		if view[i] == 1:
			lons_.append(lons[i])
			lats_.append(lats[i])
			r_.append(r[i])
			v_.append(v[i])
			times.append(dt*i)

	t_lon = (target.x*u.deg).to(u.rad)
	t_lat = (target.y*u.deg).to(u.rad)

	elevations = trigsf.get_elevations(u.Quantity(r_),
		u.Quantity(lons_),
		u.Quantity(lats_),
		t_lon,
		t_lat,
		sat.attractor.R_mean)


	data = [{'r':r_[i],
		'v':v_[i],
		'lon ssp':lons_[i],
		'lat ssp':lats_[i],
		'elevation':elevations[i],
		'time':times[i]
	} for i in range(len(lons_))]

	df = pd.DataFrame(data = data, columns = ['r', 'v', 'lon ssp', 'lat ssp', 'elevation', 'time'])

	return df



	


