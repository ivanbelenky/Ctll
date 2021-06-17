import astropy.units as u

import numpy as np
import pandas as pd

from . import coverage 

from ..utils import trigsf




def default_comm_data(sat, target, T, dt=1., lon_offset=0, **kwargs):
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
	lons,lats = sat.ssps(T,dt,lon_offset=lon_offset,**kwargs)

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

	elevations = trigsf.get_elevations(u.Quantity(r_,unit=u.km),
		u.Quantity(lons_,unit=u.rad),
		u.Quantity(lats_,unit=u.rad),
		t_lon,
		t_lat,
		sat.attractor.R_mean)


	data = [{'r':r_[i],
		'v':v_[i],
		'lon ssp':lons_[i].value,
		'lat ssp':lats_[i].value,
		'elevation':elevations[i].value,
		'time':times[i]
	} for i in range(len(lons_))]

	df = pd.DataFrame(data = data, columns = ['r', 'v', 'lon ssp', 'lat ssp', 'elevation', 'time'])

	return df



def time_over_mask(comm_df, mask, dt):
	"""Returns time of communication considering mask as minimum elevation

	Parameters
	----------
	comm_df : ~pandas.DataFrame
		communication dataframe
	mask : float
		mask value, minimum eleveation angle
	dt : float
		time interval 
	"""


	if mask.unit == u.deg:
		elevation = comm_df['elevation'].to_numpy(dtype=np.float64)
	else:
		elevation = comm_df['elevation'].to_numpy(dtype=np.float64)*180/np.pi

	elevation_mask = [e for e in elevation if e > mask.value]
	time_over = len(elevation_mask)*(dt)*3600 

	return time_over


