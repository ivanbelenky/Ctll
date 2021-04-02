import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point



class Targets(collections.abc.Set):
	def __init__(self,targets,tag=None):
		self._targets = lst = list()
		self._tag = tag if tag else "No Tag"
		for target in targets:
			if not isinstance(target,(Point,Target)):
				raise TypeError("targets must contain Point or Target object ") 
			if target not in lst:
				lst.append(target)

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
			raise TypeError("Country argument must be string")
		
		fp = 'borders-simple/TM_WORLD_BORDERS_SIMPL-0.3.shp'
		data = gpd.read_file(fp)
		data = gpd.read_file(fp)
		#data.loc[data['NAME'].where('Argentina')]

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
	def from_city(cls,city,N=50):
		if not isintance(country,str):
			raise TypeError("Country argument must be string")
		
		#TODO: download files for city borders
		fp = 'borders-simple/TM_WORLD_BORDERS_SIMPL-0.3.shp'
		data = gpd.read_file(fp)
		data = gpd.read_file(fp)
		#data.loc[data['NAME'].where('Argentina')]

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
	



		