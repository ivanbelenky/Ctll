import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point
from collections.abc import Iterable


class Targets(collections.abc.Set):
	def __init__(self,targets,tag=None):
		self._targets = lst = list()
		self._tag = tag if tag else "No Tag"
		if isinstance(targets,Target):
			lst.append(targets)
		elif isinstance(targets,Point):
			lst.append(Target(targets.x,targets.y))
		elif not isinstance(targets,Iterable):
			raise TypeError("Targets must be Target or Point, or iterable collection of Target or Point objects")
		else: 	
			for target in targets:
				if not isinstance(target,(Point,Target)):
					raise TypeError("Targets must be Target or iterable collection of Target objects") 
				if target not in lst:
					if isinstance(target,Point)
					lst.append(Target(point.x,point.y))

	@property
	def targets(self):
		return self._targets
	

	def __iter__(self):
		return iter(self.targets)

	def __contains__(self,value):
		return value in self.targets

	def __len__(self):
		len(self.elements)


	@classmethod
	def from_country(cls,country,N=50):
		if not isintance(country,str):
			raise TypeError("country argument must be string")
		if not isinstance(N,(float,int)):
			raise TypeError("N must be float or int")
		elif N<1:
			raise ValueError("N must be at least 1")
		
		fp = 'borders-simple/TM_WORLD_BORDERS_SIMPL-0.3.shp'
		data = gpd.read_file(fp)
		data = gpd.read_file(fp)

		country_row = data.loc[data['NAME'] == country]
		country_row['geometry']
		cg = country_row['geometry'] 

		lon_min,lat_min,lon_max,lat_max = cg.total_bounds
		lats = np.linspace(lat_min,lat_max, N)

		points = []
		for lat in lats:
		    N = int(np.abs(NUM*np.cos(lat*np.pi/180)))
		    lons = np.linspace(lon_min,lon_max,N)
		    for lon in lons:
		        points.append(Point(lon,lat))

		        
		candidate_points = gpd.GeoSeries(points)
		inside_points = gpd.GeoSeries([point for point in candidate_points if point.within(cg.values[0])])

		return cls(inside_points,tag=country)

	@classmethod
	def from_state(cls,state,N=50):
		if not isintance(state,str):
			raise TypeError("state argument must be string")
		if not isinstance(N,(float,int)):
			raise TypeError("N must be float or int")
		elif N<1:
			raise ValueError("N must be at least 1")

		
		fp = 'borders-states/ne_10m_admin_1_states_provinces.shp'
		data = gpd.read_file(fp)
		data = gpd.read_file(fp)

		#TODO: Check all possible column keys where the string may be found

		country_row = data.loc[data['name'] == country]
		country_row['geometry']
		cg = country_row['geometry'] 

		lon_min,lat_min,lon_max,lat_max = cg.total_bounds
		
		lats = np.linspace(lat_min,lat_max, N)

		points = []
		for lat in lats:
		    N = int(np.abs(NUM*np.cos(lat*np.pi/180)))
		    lons = np.linspace(lon_min,lon_max,N)
		    for lon in lons:
		        points.append(Point(lon,lat))

		        
		candidate_points = gpd.GeoSeries(points)
		inside_points = gpd.GeoSeries([point for point in candidate_points if point.within(cg.values[0])])

		return cls(inside_points,tag=country)

	def plot(self, use_3d=False):
		figure = plt.figure()
		x = [target.x for target in self.targets]
		y = [target.y for target in self.targets]
		if use_3d:
			centroid = [sum(x)/len(x), sum(y)/len(y)]

		else: 	
			plt.scatter(x,y,c='k',s='0.1')

		return figure


class Target(object):
	"""docstring for Target"""
	def __init__(self,lon,lat):
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
		return self._y
	



		