from poliastro import ephem

from poliastro.bodies import Earth
from poliastro.twobody import Orbit, propagation,states
from poliastro.constants import J2000
from poliastro.frames import Planes

from astropy.time import Time, TimeDelta
from astropy import time, units as u


from . import ctll
from .specs import Specifications
from .instrument import Instrument
from ..utils import trigsf

import uuid
import numpy as np
from scipy.interpolate import interp1d


SAT_ST = {
	"On":'Online',
	"Off":'Offline'
}


PROPAGATOR_DT = 50

class Sat(object):

	def __init__(self,
		state,	
		status = SAT_ST["On"],
		spec = None,
		instruments=None,
		epoch=J2000
	):

		"""Constructor.

		Parameters
		----------
		state : ~poliastro.twobody.states.ClassicalState
			State for satellite orbit
		status : string, optional
			SAT_STATUS string
		spec : ~Ctlldes.core.spec
			specifications
		epoch : ~astropy.time.Time, optional
		    Epoch, default to J2000.
		instruments : list, optional
			List of Instrument objects

		"""
		self.state = state
		self.status = status
		self.spec = spec if spec else Specifications()

		self._epoch = epoch

		if not instruments:
			self.instruments = []
		elif not isinstance(isntrumentss,list):
			self.instruments = [instruments]
		else:
			self.instruments = instruments


		self._id = uuid.uuid4()
		self._orbit = Orbit(state,epoch)
		self._Propagator = Propagator(self.orbit)

	@property
	def state(self):
		return self._state

	@state.setter
	def state(self,state):
		if not isinstance(state,states.BaseState): 
			raise Exception("Invalid state")
		else:
			self._state = state
	
	@property
	def status(self):
		return self._status
	
	@status.setter
	def status(self,status):
		if status not in SAT_ST.values():
			print("Warning")
			self._status = SAT_ST["On"]
		else:
			self._status = status
	
	@property
	def spec(self):
		return self._spec
	
	@spec.setter
	def spec(self,spec):
		if not isinstance(spec,Specifications): 
			raise Exception("Invalid satellite specifications")
		self._spec = spec

	@property
	def instruments(self):
		return self._isntruments

	@instruments.setter
	def instruments(self,instruments):
		for instr in instruments:
			if not isinstance(instr, Instrument): 
				raise Exception(f"{type(instr)} is not an Instrument object")
		self._instruments = instruments

	@property
	def id(self):
		return self._id

	@property
	def Propagator(self):
		return self._Propagator
	

	@property
	def attractor(self):
		return self.state.attractor
	
	@property
	def m(self):
		return self._spec.m
	
	@property
	def Cd(self):
		return self._spec.Cd

	@property
	def A_over_m(self):
		return self._spec._A_over_m
	
	@property
	def orbit(self):
		return self._orbit
	
	@property
	def cov_instruments(self):
		covInstr = [instr for instr in self.instruments
		if "_coverage" in set(dir(instr))]
		return self.covInstr
	
	
	
	@classmethod
	def from_orbit(
			cls,
			orbit,
			status=SAT_ST["On"],
			spec=None,
			instruments=None
	):
		"""Returns Satellite from specified orbit.

		Parameters
		----------
		orbit : poliastro.twobody.orbit.Orbit
			Specific Orbit
		status : string, optional
			SAT_STATUS string
		spec : ~Ctlldes.core.spec
			specifications
		instruments : list, optional
			List of Instrument objects

		"""

		return cls(orbit._state,status,spec,instruments,orbit.epoch)

	@classmethod
	def from_vectors(
			cls,
			r,
			v,
			attractor=Earth,
			epoch=J2000,
			plane=Planes.EARTH_EQUATOR,
			status=SAT_ST["On"],
			spec=None,
			instruments=None
	):
		"""Return Satellite from position and velocity vecotrs
		
		Parameters
		----------
		r : ~astropy.units.Quantity
			Position vector wrt attractor center.
		v : ~astropy.units.Quantity
			Velocity vector.	
		attractor : Body
			Main attractor
		epoch : ~astropy.time.Time, optional
		    Epoch, default to J2000.
		plane : ~poliastro.frames.Planes
		    Fundamental plane of the frame.
		status : string, optional
			SAT_STATUS string
		spec : ~Ctlldes.core.spec
			specifications
		instruments : list, optional
			List of Instrument objects

		"""

		orbit = Orbit.from_vectors(attractor, r, v, epoch, planes)
		return Sat.from_orbit(orbit,status,spec,instruments)


	def rv(self,T,dt=1.,method=propagation.cowell,**kwargs):
		""" Propagates orbit, T days of flight.
		
		Parameters
		----------
		T : int
			Desired time of propagation in days
		dt : float
			Size of intervals in seconds
		
		Returns
		-------
		rrvv : list
			list of tuples (rr,vv), rr and vv are Quantity 
			objects array, size=floor(T*24*3600/dt)
		
		"""
		rr,vv = self.Propagator.get_rv(T,dt,method,**kwargs)

		return rr,vv

	def ssps(self,T,dt=1.,method=propagation.cowell,**kwargs):
		""" Get subsatellite points, T days of flight.

		Parameters
		----------
		T : int
			Desired time of propagation
		dt : float
			Size of intervals in seconds
		
		Returns
		-------
		sspss : list
			list of tuples (lat,long), lat and lon are Quantity 
			objects array, size = floor(T*24*3600/dt)
		
		Exactly and very near polar orbits interpolation results in domain error.

		"""
		
		w = self.attractor.angular_velocity
		lats,lons = self.Propagator.get_ssps(T,dt,w,method=propagation.cowell,**kwargs)
		
		return lats,lons

	def update_instruments(self,instruments):
		"""Updates satellite's insruments

		
		"""
		self.instruments = instruments

	def update_status(self,newStatus):
		"""Updates satellite's status
		
		Parameters
		----------
		newStatus : string
			SAT_STATUS string
		"""

		self.status = newStatus


	def info(self):
		print(f"id : {self.id}\t[a ecc inc raan argp nu] : "+
			f"[{self.orbit.a:.1f} {self.orbit.ecc:.1f} {self.orbit.inc:.1f}"+
			f" {self.orbit.raan:.1f} {self.orbit.argp:.1f} {self.orbit.nu:.1f}]"+
			f"\t\tstatus : {self.status}")
		 





class Propagator(object):

	def __init__(
		self,
		orbit,
	 	T = 1,
	 	method = propagation.cowell,	   
		**kwargs
	):
		
		"""Propagator object builder. Build coordinates"""
		self._DT = PROPAGATOR_DT 


		self._orbit = orbit 
		self._T = T 
		self._kwargs = kwargs
		self._method = method
		#TODO: check if this is ok unpacking.

		self._tofs = None
		self._setcoords()
		

		#hardsetted to optimize performance on interpolation



	@property
	def orbit(self):
		return self._orbit

	@orbit.setter
	def orbit(self,orbit):
		if not isinstance(orbit, Orbit):
		#TODO, check if orbits are recognized by Orbit	
			raise Exception(f"orbit must be {type(Orbit)}")
	
	@property
	def T(self):
		return self._T

	@T.setter
	def T(self,T):
		if not isinstance(T,float) and not isinstance(T,int): 	
			raise Exception("T has to be float or int ")
		if T <= 0: raise Exception("T must be positive and != 0")
		self._T = T
	
	@property
	def kwargs(self):
		return self._kwargs
	


	@property
	def method(self):
		return self._method
	
	@method.setter
	def method(self,method):
		if method.func_name not in propagation.ALL_PROPAGATORS:
			raise Exception("Unknown propagation method") 


	@property
	def dt(self):
		return self._dt
	
	@dt.setter
	def dt(self,dt):
		if not isinstance(dt,float): raise Exception("dt has to be float")
		if dt <= 0: raise Exception("dt must be positive and !=0")	
		self._dt = dt

	@property
	def coords(self):
		return self._coords



	def _setcoords(self):
		"""Update tofs and coords. 
	
		Note
		----
		After updating T or kwargs this method is going to 
		be called in order to update the times of flight of the
		propagator as well as the resulting coordinate for the new
		timeframe or perturbative force.
		"""


		self.tofs = TimeDelta(np.linspace(0,self.T*24*3600*u.s,
		 num=int(self.T*24*3600/self._DT)))
		
		self._coords = propagation.propagate(self.orbit,
			self.tofs,method=self.method,**self.kwargs)


	def get_rv(self,T,dt=1.,method=propagation.cowell,**kwargs):

		"""Get position and velocity for specified Time of flight
		
		Parameters
		----------
		T : float
			Time of flight to propagate.
		dt : float, optional
			time interval between interpolation. 
		method : callable, optional
			Propagation method, default to cowell.

		Returns
		-------
		r, v : tuple 
			Tuple of Quantity arrays. 

		Note
		----
		The keyword arguments correspond to the perturbative force
		arguments. In the submodule perturb, combination of 
		usual perturbation forces are defined. 

		"""



		# if not isinstance(T,float): 
		# 	raise Exception("T has to be float")
		# if T <= 0: raise Exception("T must be positive and != 0")


		# if not isinstance(dt,float): raise Exception("dt has to be float")
		# if dt <= 0: raise Exception("dt must be positive and !=0")	



		flag = False
		if kwargs:
			if kwargs != self.kwargs:
				flag = True
				self.kwargs = kwargs
		if T > self.T:
			flag = True 
			self.T = T
		if method != self.method:
			flag = True
			self.method = method

		if flag:
			self._setcoords()

		tofs = TimeDelta(np.linspace(0, T*24*3600*u.s, num=int(T*24*3600/dt)))


		ephemerides = ephem.Ephem(self.coords,self.tofs, Planes.EARTH_EQUATOR)
		
		return ephemerides.rv(tofs)
	

	def get_ssps(self,T,dt,w,method=propagation.cowell,**kwargs):
		"""Return subsatellite points"""

		flag = False
		if kwargs:
			if kwargs != self.kwargs:
				flag = True
				self.kwargs = kwargs
		if T > self.T:
			flag = True 
			self.T = T
		if method != self.method:
			flag = True
			self.method = method

		if flag:
			self._setcoords()

		ephemerides = ephem.Ephem(self.coords,self.tofs, Planes.EARTH_EQUATOR)
		rr,vv = ephemerides.rv()
		print("le shape is ", rr.shape)
		rs,lats,lons = trigsf.c2s(rr[:,0].value, rr[:,1].value, rr[:,2].value)

		tofs = np.linspace(0,self.T*24*3600*u.s,num=int(self.T*24*3600/self._DT))
		
		#de-rotation
		lons = np.array([((lon*u.rad)-(t*w)).value%(2*np.pi) for lon,t in zip(lons,tofs)])*u.rad
		lats = lats*u.rad

		#this manipulation is needed to fold the interpolation 
		adder = 2*np.pi * u.rad
		for i in range(1,len(lons)):
			if np.abs(lons[i].value-lons[i-1].value) > 1:
				sgn=np.sign(lons[i].value-lons[i-1].value)
				lons[i:] += -sgn*adder
		
		#interpolation of manipulated and folded points 
		flats = interp1d(tofs,lats)
		flons = interp1d(tofs,lons)
		tofs = np.linspace(0,T*3600*24*u.s,int(T*3600*24/dt))
		lats = flats(tofs)
		lons = flons(tofs) % (2*np.pi)
		
		return lats,lons


