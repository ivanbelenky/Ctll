import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point
import collections.abc
from collections.abc import Iterable


from astropy import units as u

import os
import sys

fp_country = sys.path[0] + '/CtllDes/targets/borders-simple/TM_WORLD_BORDERS_SIMPL-0.3.shp'
fp_state = sys.path[0] + '/CtllDes/targets/borders-states/ne_10m_admin_1_states_provinces.shp'

#TODO: update docstrings. Redo the functions adding the filepath in order to eliminate the database for borders

TOL = 10**-5 

class Targets(collections.abc.Set):
	""""""

	def __init__(self,tgts,tag=None):
		""""""
		self._targets = lst = list()
		self._tag = tag if tag else "No Tag"

		if isinstance(tgts,Target):
			lst.append(tgts)
		elif isinstance(tgts,Point):
			lst.append(Target(tgts.x,tgts.y))
		elif not isinstance(tgts,Iterable):
			raise TypeError("Targets must be Target or Point, or iterable collection of Target or Point objects")
		else: 	
			for target in tgts:
				if not isinstance(target,(Point,Target)):
					raise TypeError("Targets must be Target or iterable collection of Target objects") 
				if target not in lst:
					if isinstance(target,Point):
						lst.append(Target(target.x,target.y))
					else:
						lst.append(target)

	@property
	def targets(self):
		return self._targets
	

	def __iter__(self):
		return iter(self.targets)

	def __contains__(self,value):
		return value in self.targets

	def __len__(self):
		return len(self.targets)

	@property
	def coords(self):
		lst = np.array([[target.x,target.y] for target in self.targets])
		return lst
	

	@classmethod
	def from_country(cls,country,N=50):
		""""""
		if not isinstance(country,str):
			raise TypeError("country argument must be string")
		if not isinstance(N,(float,int)):
			raise TypeError("N must be float or int")
		elif N<1:
			raise ValueError("N must be at least 1")
		

		fp = fp_country
		data = gpd.read_file(fp)
		
		country_row = data.loc[data['NAME'] == country]
		country_row['geometry']
		cg = country_row['geometry'] 

		lon_min,lat_min,lon_max,lat_max = cg.total_bounds
		lats = np.linspace(lat_min,lat_max, N)

		points = []
		for lat in lats:
		    new_N = int(np.abs(N*np.cos(lat*np.pi/180)))
		    lons = np.linspace(lon_min,lon_max,new_N)
		    for lon in lons:
		        points.append(Point(lon,lat))
		
        
		candidate_points = gpd.GeoSeries(points)
		inside_points = gpd.GeoSeries([point for point in candidate_points if point.within(cg.values[0])])
		return cls(inside_points,tag=country)

	@classmethod
	def from_state(cls,state,N=50):
		""""""
		if not isinstance(state,str):
			raise TypeError("state argument must be string")
		if not isinstance(N,(float,int)):
			raise TypeError("N must be float or int")
		elif N<1:
			raise ValueError("N must be >= 1")

		
		fp = fp_state
		data = gpd.read_file(fp)
		data = gpd.read_file(fp)


		#TODO: Check all possible column keys where the string may be found

		state_row = data.loc[data['name'] == state]
		state_row['geometry']
		cg = state_row['geometry'] 

		lon_min,lat_min,lon_max,lat_max = cg.total_bounds
		
		lats = np.linspace(lat_min,lat_max, N)
		points = []
		for lat in lats:
		    new_N = int(np.abs(N*np.cos(lat*np.pi/180)))
		    lons = np.linspace(lon_min,lon_max,new_N)
		    for lon in lons:
		        points.append(Point(lon,lat))

		        
		candidate_points = gpd.GeoSeries(points)
		inside_points = gpd.GeoSeries([point for point in candidate_points if point.within(cg.values[0])])

		return cls(inside_points,tag=state)

	def plot(self, use_3d=False):
		figure = plt.figure()
		x = [target.x for target in self.targets]
		y = [target.y for target in self.targets]
		if use_3d:
			centroid = [sum(x)/len(x), sum(y)/len(y)]
			#TODO
			raise NotImplemented
		else: 	
			plt.scatter(x,y,c='k',s=0.5)

		return figure


class Target(object):
	""""""
	def __init__(self,lon,lat):
		"""Targets represents a point in a sphere surface, longitude
		and latitude respectively. This values must be provided in 
		degrees.

		Parameters
		----------
		lon : float | int
			longitude 
		lat : float | int
			latitude 
		 
		"""

		#TODO: fix lat lon in coverage or here.

		self._lat = lat 
		self._lon = lon 

	@property
	def lat(self):
		return self._lat


	@property
	def lon(self):
		return self._lon

	@property
	def x(self):
		return self._lon

	@property
	def y(self):
		return self._lat
	
	def __eq__(self, other):
		if np.abs(self.x-other.x) < TOL and np.abs(self.y-other.y) < TOL:
			return True
		else:
			return False

	def __ne__(self,other):
		if np.abs(self.x-other.x) > TOL or np.abs(self.y-other.y) < TOL:
			return True
		else:
			return False
