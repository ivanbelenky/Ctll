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
cities_fp = sys.path[0] + '/CtllDes/targets/location-cities/worldcities.csv'


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
	
	@property
	def lons(self):
		return [t.lon for t in self.targets]

	@property	
	def lats(self):
		return [t.lat for t in self.targets]


	@classmethod
	def from_country(cls,country,N=50):
		"""Returns Targets built from country.

		Parameters 
		----------
		country : string
			Name, FIPS, ISO2 and ISO3 codes are valid entries	
		N : int, optional
			proportional to square root of #targets.

		If not found a ValueError exception will be thrown.  

		"""

		if not isinstance(country,str):
			raise TypeError("country argument must be string")
		if not isinstance(N,(float,int)):
			raise TypeError("N must be float or int")
		elif N<1:
			raise ValueError("N must be at least 1")
		

		fp = fp_country
		data = gpd.read_file(fp)
		
		names = data['NAME'].values
		fips = data['FIPS'].values
		iso2 = data['ISO2'].values
		iso3 = data['ISO3'].values 

		if country in names:
			country_row = data.loc[data['NAME'] == country]
			cg = country_row['geometry']
		elif country in fips:
			country_row = data.loc[data['FIPS'] == country]
			cg = country_row['geometry']
		elif country in fips:
			country_row = data.loc[data['ISO2'] == country]
			cg = country_row['geometry']
		elif country in fips:
			country_row = data.loc[data['ISO3'] == country]
			cg = country_row['geometry']
		else:
			raise ValueError(f"{country}: invalid Name, FIPS, ISO2 or ISO3 code")


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
		"""Returns Targets built from state.

		Parameters 
		----------
		state : string
			Name, English Name, ISO or ADM1 codes are valid entries	
		N : int, optional
			proportional to square root of #targets.

		If not found a ValueError exception will be thrown.  

		"""	

		if not isinstance(state,str):
			raise TypeError("state argument must be string")
		if not isinstance(N,(float,int)):
			raise TypeError("N must be float or int")
		elif N<1:
			raise ValueError("N must be >= 1")

		
		fp = fp_state
		data = gpd.read_file(fp)

		adms = data['adm1_code'].values
		isos = data['iso_3166_2'].values
		names = data['name'].values
		alt_names = data['name_alt'].values
		en_names = data['name_en'].values


		if state in adms:
			state_row = data.loc[data['adm1_code'] == state]
			cg = state_row['geometry']
		elif state in isos:
			state_row = data.loc[data['iso_3166_2'] == state]
			cg = state_row['geometry']
		elif state in names:
			state_row = data.loc[data['name'] == state]
			cg = state_row['geometry']
		elif state in en_names:
			state_row = data.loc[data['name_en'] == state]
			cg = state_row['geometry']
		elif state in alt_names:	
			state_row = data.loc[data['name_alt'] == state]
			cg = state_row['geometry']
		else:
			raise ValueError(f"{state}: invalid Name, English Name, ISO or ADM1 code")


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
	"""
	"""

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
	
	@classmethod
	def from_city(cls, city, country = None):
		"""Return Target from city name and iso code if specified
		
		Parameters
		----------
		city : string 
			City Name, Name_ASCII or admin_name
		country : string, optional
			Country Name, or  ISO2/ISO3 code

		If country is not specified and multiple entries are found with
		specified city name, the first value will be taken into account.
		
		"""


		cities_df = pd.read_csv(cities_fp)
		
		names = cities_df['city'].values	
		ascii_names = cities_df['city_ascii'].values
		admin_names = cities_df['admin_name'].values
		ids = cities_df['id'].values

		c_names = cities_df['country'].values
		c_iso2  = cities_df['iso2'].values
		c_iso3  = cities_df['iso3'].values


		if city in names:
			_city = cities_df.loc[cities_df['city'] == city]
		elif city in ascii_names:
			_city = cities_df.loc[cities_df['city_ascii'] == city]
		elif city in admin_names:
			_city = cities_df.loc[cities_df['admin_name'] == city]
		elif city in ids:
			_city = cities_df.loc[cities_df['id'] == city]
		else:
			raise ValueError(f"{city}: invalid city Name, ASCII Name, id or admin Name")

		if not country:
			lon,lat = _city[['lng','lat']].values[0]
		else: 
			if country in c_names:
				lon, lat = _city.loc[_city['country'] == country ][['lng','lat']].values[0]
			elif country in c_iso2:
				lon, lat = _city.loc[_city['iso2'] == country ][['lng','lat']].values[0]
			elif country in c_iso3:
				lon, lat = _city.loc[_city['iso3'] == country ][['lng','lat']].values[0]
			else:
				raise ValueError(f"{country}: invalid country Name or ISO2/ISO3 Code for {city} city")

		return cls(lon,lat)


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
