from ..core import ctll, satellite, instrument
from ..targets.targets import Target
from ..utils import trigsf


import collections.abc
from collections.abc import Iterable
from astropy import units as u



import numpy as np
import pandas as pd

try:
    from functools import cached_property  # type: ignore
except ImportError:
    from cached_property import cached_property  # type: ignore

TOL = 1**(-10)



def get_view(lons, lats, r, target, R):
	"""Get string indicating in view considering maximum 
	field of view possible, i.e. elevation = 0 

	Parameters
	----------
	lons : ~astropy.units.quantity.Quantity
		longitudes in radians
	lats : ~astropy.units.quantity.Quantity
			latitudes in radians
	r : ~astropy.units.quantity.Quantity
		position, distance
	target : CtllDes.targets.targets.Target
		desired target
	R : ~astropy.units.quantity.Quantity
		attractors mean distance
	"""
	t_lon = (target.x*u.deg).to(u.rad)
	t_lat = (target.y*u.deg).to(u.rad)

	lams_0 = trigsf.get_lam_0(r,R)

	angles = trigsf.get_angles(lons, lats, t_lon, t_lat)

	cov = []
	for a,l in zip(angles,lams_0):
		if a < l: 
			cov.append(1)
		else:
			cov.append(0)

	return cov


def push_broom(FOV_vert, lons, lats, r, target, R):
	view = get_view(lons, lats, r, target, R)
	
	cov = []
	for i in range(1,len(view)):
		if view[i] == 1:
			radii = np.sqrt(np.sum(r[i]**2))

			t_lon = (target.x*u.deg).to(u.rad)-lons[i]
			t_lat = (target.y*u.deg).to(u.rad)-lats[i]
							
			inc = np.arctan2(np.sin(lons[i]-lons[i-1])*np.sin(lats[i]-lats[i-1]),np.cos(lats[i-1]))
			cos_inc = np.cos(i)
			sin_inc = np.sin(i)
			angle = trigsf.get_angles(0,0,lons[i]-lons[i-1],lats[i]-lats[i-1])

			t_lon = np.arctan2(cos_inc*np.sin(t_lon)*np.sin(t_lat)+sin_inc*np.cos(t_lat),
			 np.cos(t_lon)*np.sin(t_lat))
			t_lat = np.arccos(-sin_inc*np.sin(t_lon)*np.sin(t_lat)+cos_inc*np.cos(t_lat)) 

			lam = trigsf.get_lam(r[i],FOV_vert,R)

			if t_lon < angle/2 and t_lat < lam:
				cov.append(1)
			else:
				cov.append(0)
		else:
			cov.append(0)
			pass

	return cov



#TODO: implement
def symmetric_with_roll():
	pass



def symmetric_disk(FOV_min,FOV_max,lons,lats,r,target,R):
	"""coverage method.

	Disk of coverage centered on subsatellite point.
	
	Parameters
	----------
	FOV_min : ~Astropy.units.quantity.Quantity
		minimum field of view in radians
	FOV_max : ~Astropy.units.quantity.Quantity
		maximum field of view in radians

	* : default coverage parameters
		help(CtllDes.request.coverage.Instrument.coverage) for more
		info.

	"""
		
	#HACK: there were 2 options, do this nasty thing here 
	#or in targets.py I had to make that call.

	if FOV_max < FOV_min:
		raise ValueError("Wrong FOV ordering")


	lams_min = trigsf.get_lam(r,FOV_min,R)
	lams_max = trigsf.get_lam(r,FOV_max,R)

	angles = trigsf.get_angles(lons,lats,(target.x*u.deg).to(u.rad),
		(target.y*u.deg).to(u.rad))

	cov = []
	for angle,lam_min,lam_max in zip(angles,lams_min,lams_max):
		if  lam_min <= angle <= lam_max:
			cov.append(1)
		else:
			cov.append(0)
	
	return cov




def symmetric(FOV,lons,lats,r,v,target,R):
	"""Circle of coverage centered on ssp"""
		
	#HACK: there were 2 options, do this nasty thing here 
	#or in targets.py I had to make that call.



	lams = trigsf.get_lam(r,FOV,R)
	angles = trigsf.get_angles(lons,lats,(target.x*u.deg).to(u.rad),
		(target.y*u.deg).to(u.rad))

	cov = []
	for angle,lam in zip(angles,lams):
		if angle <= lam:
			cov.append(1)
		else:
			cov.append(0)
	
	return cov



COV_METHODS = [symmetric, symmetric_disk, symmetric_with_roll]



class Coverages(collections.abc.Set):
	"""Container for Coverage objects. Frame where you can 
	get all the data from the coverage analysis. It is meant
	to be created with the classmethods. 
	"""
	def __init__(self,covs,tag=None):
		self._covs = lst = list()
		self._tag = tag if tag else "No Tag"

		if isinstance(covs,Coverage):
			lst.append(covs)
		elif not isinstance(covs,Iterable):
			raise TypeError("covs must be Coverage, or iterable collection of Coverage objects")
		else: 	
			for cov in covs:
				if not isinstance(cov,Coverage):
					raise TypeError("covs must be a collection of Coverage objects") 
				if cov not in lst:
					lst.append(cov)

		self._targets = {(covv.target.lon,covv.target.lat) for covv in self.covs}
		self._sats_id = {covv.sat_id for covv in self.covs }


	@property
	def covs(self):
		return self._covs
	
	@property
	def targets(self):
		tgts = [Target(ll[0],ll[1]) for ll in self._targets]
		return tgts
	
	@property
	def sats_id(self):
		return self._sats_id
	


	def __iter__(self):
		return iter(self.covs)

	def __contains__(self,value):
		return value in self.covs

	def __len__(self):
		len(self.covs)

	def __str__(self):
		return self.tag

	#TODO: how to create dataframes in python.
	def to_df(self):		
		"""Returns dataframe with all the merit figures
		form the Coverage contained.

		"""

		df = pd.DataFrame(columns=['T',
									'dt',
									'Satellite ID',
									'Target',
									'accumulated',
									'mean gap light',
									'mean gap dark',
									'response time',
									'average time gap',
                          			'max gap'])

		data = [{'T': cov.T,
				'dt': cov.dt,
				'Satellite ID': cov.sat_id, 
	        	'Target':(cov.target.lon,cov.target.lat),
    	    	'accumulated': cov.accumulated,
        		'mean gap light': cov.mean_gap_light,
       			'mean gap dark': cov.mean_gap_dark,
        		'response time': cov.response_time,
        		'average time gap': cov.avg_time_gap,
        		'max gap': cov.max_gap} for cov in self.covs]

		df = df.append(data,ignore_index=True)

		return df
		

	@classmethod
	def from_ctll(cls,ctll,targets,T,dt=1.,**kwargs):
		"""Get coverages from constellation.
		
		Parameters
		----------
		ctll : CtllDes.core.ctll.Ctll
			CtllDes constellation object
		targets : CtllDes.targets.targets.Targets
			Desired targets to build coverages
		T : float | int
			Desired time of coverage in days
		dt : float | int, optional
			time of sampling in seconds

		Returns
		-------
		Coverages object containing only the coverage from
		instruments that have _coverage function implemented.
		
		TODO: provide duck typing or not for instruments in ctll?
		"""

		ctll_covs = []
		for sat in ctll.sats:
			print(f'Satellite {ctll.sats.index(sat)+1} of {ctll.N}')
			try:
				sat_covs = Coverages.from_sat(sat,targets,T,dt,f=True,**kwargs)
			except Exception as e:
				print(e)
				pass
			else:
				ctll_covs += sat_covs

		if not len(ctll_covs): 
			raise Exception("Constellation has no Coverage Instruments")

		return Coverages(ctll_covs,tag=ctll.__str__())


	@classmethod
	def from_sat(cls,sat, targets,T, dt=1.,f=False,**kwargs):
		"""Build list of coverage objects from satellite
	
		Parameters
		----------
		sat : ~CtllDes.core.sat
			sat object
		targets : ~CtllDes.targets.targets.Targets 
			Desired targets of coverage
		T : float
			Desired Time of analysis in days.
		dt : float | int, optional
			time of sampling in seconds
		f :  boolean
			Not to be modfied from default state
		Returns
		-------
		Coverages : list
			List of Coverage objects, one for each target and  
			coverage instrument.
		
		Not the best practice to return two different types,
		depending on internal or external use.
		but since the use of this library is quite reduced. This 
		will be fixed later.

		"""
		if isinstance(targets,Target):
			targets = [targets]

		cov_instruments = [instr for instr in sat.cov_instruments]
		if not len(cov_instruments):
			raise Exception("No coverage instruments found on" +
				f" satID : {sat.id}")	

		lons,lats = sat.ssps(T,dt,**kwargs)
		r,v = sat.rv(T,dt,**kwargs)


		total = len(cov_instruments)*len(targets)
		count = 1
		sat_coverages = []
			
		for instr in cov_instruments:
			for target in targets:
				#cov = isCovered(lons,lats,r,target,sat.attractor.R_mean,instr.coverage())
				cov = instr.coverage(lons,lats,r,v,target,sat.attractor.R_mean)
				sat_coverages.append(Coverage(cov,target,T,dt,sat.id))
				print(f'target {target.x}° {target.y}°. {count} of {total}')
				count += 1

		if f:
			return sat_coverages
		else:
			return Coverages(sat_coverages,tag=sat.__str__())



	def collapse_sats(self,id_lst=None):
		""" Collapse Coverages object and create new one from specified
		satellites. If non specified, it will colapse all sats.
		
		Parameters
		----------
		id_lst : list, optional
			Satellites ID list
		"""		
		lst = id_lst if id_lst else self.sats_id

		sat_tgt = []

		for target in self.targets:

			sat_tgt.append([covv for covv in self.covs 
				if (covv.sat_id in lst and covv.target == target)])

		Covs = []
		for sats in sat_tgt:
			cov = sats[0]
			newcov = np.max(np.array([sat.cov for sat in sats]),axis=0)
			print(newcov)
			Covs.append(Coverage(newcov,cov.target,cov.T,cov.dt))

		return Coverages(Covs,tag=f'collapsed satellites: {lst}')


	

		

class Coverage(object):
	"""Coverage is the main class of coverages.py, composed 
	of the data needed to make all the desired analysis in order to get the 
	merit figures for a specific target, during a time of flight.
	"""
	
	def __init__(self,
		cov,
		target,
		T,
		dt,
		sat_id=None,
	):
		"""Constructor for Coverage
		
		Parameters
		----------
		cov : Iterable
			Iterable collection of 1 or 0 values representing
			view and not in view.
		target : CtllDes.targets.targets.Target
			Target related to cov collection.
		T : float | int
			Time of flight in days
		dt : float | int
			time interval in seconds
		sat_id : uuid.UUID, optional
			Satellite id related to the cov collection.
		"""

		self._cov = cov
		self._target = target
		self._sat_id = sat_id if sat_id else ''
		self._T = T
		self._dt = dt

		self._set_merit_figures()



	@property
	def cov(self):
		return self._cov

	@property
	def sat_id(self):
		return self._sat_id
	
	
	@property
	def target(self):
		return self._target
	
	@property
	def T(self):
		return self._T
	
	@property
	def dt(self):
		return self._dt
	
		
	def _set_merit_figures(self):
		self._accumulated = self._accum()
		self._mean_gap_light = self._mean_gap(1)
		self._mean_gap_dark = self._mean_gap(0)
		
		resp,avg,max_gap = self._resp_avg_max()
		
		self._response_time = resp
		self._avg_time_gap = avg
		self._max_gap = max_gap

	@property
	def accumulated(self):
		"""Accumulated time of view [s]"""
		return self._accumulated

	@property
	def mean_gap_light(self):
		"""Mean gap time in view [s]"""
		return self._mean_gap_light*self.dt

	@property 
	def mean_gap_dark(self):
		"""Mean gap time not in view [s]"""
		return self._mean_gap_dark*self.dt

	@property
	def response_time(self):	
		"""Average time to view target. This is
		obtained by sampling the covs attribute
		and calculating the time interval to next
		view.
	 	"""

		return self._response_time*self.dt

	@property
	def avg_time_gap(self):
		"""Average time of the dark gap. Is obtained
		sampling randomly the covs, if it is a gap
		it adds that gap_time. It repeats this a couple
		of times and calculates the the mean [s].
		"""

		return self._avg_time_gap*self.dt

	@property
	def max_gap(self):
		"""Maximum dark gap [s]"""
		return self._max_gap*self.dt
	

	def _accum(self):
		"""Integrated time of light (in view) [s]"""
		return self.dt*sum(self.cov)


	def _mean_gap(self,view):
		"""Returns the mean gap in view."""

		gap = 0 
		switch = 0 
		c_gap = 0
		
		for c in self.cov:
			if switch == 0 and c == view:
				switch = 1
				c_gap += 1
				gap += 1
			elif switch == 1 and c != view:
				switch = 0
			elif switch == 1 and c == view:
				gap += 1
			elif switch == 0 and c != view:
				continue

		if c_gap == 0:
			return 0
		else:
			return gap/c_gap	

	
	
	def _resp_avg_max(self):
		"""Returns response time, average time gap, 
		and max gap in seconds"""
		
		idx = []
		resp_time = 0
		time_gap = 0
		gap = 0
		switch = 0
		max_gap = 0
		N = len(self.cov)

		for c in self.cov:
			if switch == 0 and c == 0:
				switch = 1
				gap += 1
			elif switch == 1 and c != 0:
				switch = 0
				resp_time += (gap+1)*gap/2
				time_gap += gap*gap
				if gap > max_gap:
					max_gap = gap
				gap = 0
			elif switch == 1 and c == 0:
				gap += 1
			elif switch == 0 and c != 0:
				continue

		return resp_time/N,time_gap/N,max_gap

	

	def __add__(self,other):
		if self.target != other.target:
			raise ValueError("Targets must be equal")
		elif self.T-other.T<TOL or self._dt-other.dt<TOL:
			raise ValueError("T and dt must be equal")
		elif len(self.cov) != len(other.cov):
			raise ValueError("Tolerance may be ill-defined")
		else:
			new_cov = np.array(self.cov) | np.array(other.cov)
			return Coverage(list(new_cov),self.target,self.T,self.dt)
		

	def plot_lapida(self, **kwargs):
		"""Grafico de lapidas. 

		Funcion escalon en los intervalos donde el objetivo es cubierto.
		En esencia es la misma funcion step de pyplot. En caso que quiera ser 

		Recibe: 
			**kwargs: argumentos propios de la funcion matplotlib.pyplot.step
			(para mas informacion recurrir a documentacion oficial de matplotlib
			https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.pyplot.step.html)
		Devuelve:
			 
		"""
		
		merged = self.merge()
		x = np.array( [self.dt*i/3600 for i in range(len(merged))] ) 
		y = np.array(merged)

		fig = plt.figure()
		plt.step(x,y,**kwargs)
		plt.title(str(self.acumulado()))

		return fig







