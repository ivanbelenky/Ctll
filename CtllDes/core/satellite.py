from poliastro import ephem

from poliastro.bodies import Earth
from poliastro.twobody import Orbit, propagation,states
from poliastro.constants import J2000
from poliastro.frames import Planes


from poliastro.constants import rho0_earth, H0_earth
from poliastro.core.perturbations import atmospheric_drag_exponential, third_body, J2_perturbation

from astropy.time import Time, TimeDelta
from astropy import time, units as u


from .specs import Specifications
from .instrument import Instrument
from ..utils import trigsf

import uuid
import numpy as np
from scipy.interpolate import interp1d
import time 
import copy

from poliastro.core.util import jit


SAT_ST = {
	"On":'Online',
	"Off":'Offline'
}




PROPAGATOR_DT = 300

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
		elif not isinstance(instruments,list):
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
			print("Warning, notImplemented")
			self._status = SAT_ST["On"]
		else:
			self._status = status
	
	@property
	def spec(self):
		return self._spec
	
	@spec.setter
	def spec(self,specif):
		if not isinstance(specif,Specifications): 
			raise Exception("Invalid satellite specifications")
		self._spec = specif

	@property
	def instruments(self):
		return self._instruments

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
		if 'coverage' in set(dir(instr))]
		return covInstr

	@property
	def a(self):
		return self.orbit.a

	@property
	def ecc(self):
		return self.orbit.ecc
	
	@property
	def inc(self):
		return self.orbit.inc
		
	@property
	def raan(self):
		return self.orbit.raan
	
	@property
	def argp(self):
		return self.orbit.argp

	@property
	def nu(self):
		return self.orbit.nu
	
	
	
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

		orbit = Orbit.from_vectors(attractor, r, v, epoch, plane)
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

		kw = self._parse_kwargs(**kwargs)
		
		rr,vv = self.Propagator.get_rv(T,dt,method,**kw)
		

		return rr,vv


	def ssps(self,T,dt=1.,method=propagation.cowell,lon_offset=0, **kwargs):
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
			list of tuples (lon,lat), lat and lon are Quantity 
			objects array, size = floor(T*24*3600/dt)
		
		Exactly and very near polar orbits interpolation results in domain error.

		"""
		

		kw = self._parse_kwargs(**kwargs)		
		w = self.attractor.angular_velocity
		lons,lats = self.Propagator.get_ssps(T,dt,w,method,lon_offset,**kw)
		
		return lons,lats


	def ssps_from_r(self,r,T,dt=1.,method=propagation.cowell,lon_offset=0,**kwargs):
		"""Same as ssps but builds from rr
		"""
		
		kw = self._parse_kwargs(**kwargs)
		w = self.attractor.angular_velocity
		lons,lats = self.Propagator.get_ssps_from_rr(T,dt,w,r,method,lon_offset,**kw)

		return lons,lats
	
	def __str__(self):
		return f'{self.orbit}'


	def _parse_kwargs(self,**kw):
		
		drag = 0
		J2 = 0
		fad = False
		args = {} 

		for k in kw:
			if k == 'drag':
				if kw[k]:
					drag = 1
				continue
			if k == 'J2':
				if kw[k]:
					J2 = 1
				continue
			if k == 'ad':
				fad = True
				continue

			args[k]=kw[k]

		if fad:
			ad_args = copy.deepcopy(args)

		if J2:
			args['J2'] = self.attractor.J2.value
			args['R'] = self.attractor.R.to(u.km).value
		if drag:
			if self.attractor != Earth:
				drag = False
				print(f"No atmospheric drag model for {self.attractor} ")
				pass
			else:
				args['R'] = self.attractor.R.to(u.km).value
				args['H0'] = H0
				args['rho0'] = rho0
				args['C_D'] = self.spec.Cd.value
				args['A_over_m'] = self.spec.A_over_m.to_value(u.km** 2/u.kg)

		if J2 and drag:
			args['ad'] = J2_and_drag
			return args 
		elif J2:
			args['ad'] = J2_perturbation
			return args
		elif drag:
			args['ad'] = atmospheric_drag_exponential
			return args
		else:
			return {}

		#TODO: interpret ad if given, alongside J2 and drag.


	def update_spec(self,spec):
		"""Updates satellite's specs
		
		Parameters
		----------
		spec : ~CtllDes.core.specs.Specifications
			new Satellite specifications
		"""
		self.spec = spec


	def update_instruments(self,instruments,f=False):
		"""Updates satellite's insruments
		
		Parmaeters
		----------
		instruments : ~CtllDes.core.instrument.Instrument 
			instrument or list containing instruments
		f : boolean
			force parameter, if True, satellite's instruments
			will be replaced. If false they will be appended.
		"""

		itr = copy.deepcopy(instruments)
		if isinstance(instruments, list):
			if not f:
				itr += self.instruments 
			self.instruments = itr
		elif isinstance(instruments, Instrument):
			itr = [itr]
			if not f:
				itr += self.instruments
			self.instruments = itr	
		else:
			pass
	
		
	def update_status(self,status):
		"""Updates satellite's status
		
		Parameters
		----------
		status : string
			SAT_STATUS string
		"""

		self.status = status


	def info(self):
		print(f"id : {self.id}\t[a ecc inc raan argp nu] : "+
			f"[{self.orbit.a:.1f} {self.orbit.ecc:.1f} {self.orbit.inc:.1f}"+
			f" {self.orbit.raan:.1f} {self.orbit.argp:.1f} {self.orbit.nu:.1f}]"+
			f"\t\tstatus : {self.status}")
		 





class Propagator(object):

	def __init__(
		self,
		orbit,
	 	T = 0.01,
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

		self.tofs = None
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
		if method not in propagation.ALL_PROPAGATORS:
			raise Exception("Unknown propagation method") 
		self._method = method


	@property
	def dt(self):
		return self._dt
	
	@dt.setter
	def dt(self,dt):
		if not isinstance(dt,float): raise Exception("dt has to be float")
		if dt <= 0: raise Exception("dt must be positive and !=0")	
		self._dt = dt

	



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

		self.coords = propagation.propagate(self.orbit,
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
		if kwargs != self.kwargs:
			flag = True
			self._kwargs = kwargs
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
	

	def get_ssps_from_rr(self, T, dt, w, rr , method=propagation.cowell, lon_offset=0, **kwargs):
		"""Return subsatellite points"""
		
		flag = False		
		if kwargs != self.kwargs:
			flag = True
			self._kwargs = kwargs
		if T > self.T:
			flag = True 
			self.T = T
		if method != self.method:
			flag = True
			self.method = method

		if flag:
			self._setcoords()

		rs,lats,lons = trigsf.c2s(rr[:,0].value, rr[:,1].value, rr[:,2].value)

		tofs = np.linspace(0,self.T*24*3600*u.s,num=int(self.T*24*3600/dt))
		
		#de-rotation
		lons = np.array([((lon*u.rad)-(t*w)+lon_offset).value%(2*np.pi) for lon,t in zip(lons,tofs)])*u.rad
		lats = lats*u.rad

		return lons,lats


	def get_ssps(self,T,dt,w,method=propagation.cowell,lon_offset=0,**kwargs):
		"""Return subsatellite points"""

		flag = False		
		if kwargs != self.kwargs:
			flag = True
			self._kwargs = kwargs
		if T > self.T:
			flag = True 
			self.T = T
		if method != self.method:
			flag = True
			self.method = method

		if flag:
			self._setcoords()


		#ephemerides = ephem.Ephem(self.coords,self.tofs, Planes.EARTH_EQUATOR)
		rr,vv = self.get_rv(T,dt,method=propagation.cowell, **kwargs)#ephemerides.rv()
		rs,lats,lons = trigsf.c2s(rr[:,0].value, rr[:,1].value, rr[:,2].value)

		tofs = np.linspace(0,self.T*24*3600*u.s,num=int(self.T*24*3600/dt))
		
		#de-rotation
		lons = np.array([((lon*u.rad)-(t*w)+lon_offset).value%(2*np.pi) for lon,t in zip(lons,tofs)])*u.rad
		lats = lats*u.rad


		#this manipulation is needed to fold the interpolation 
		# adder = 2*np.pi * u.rad
		# for i in range(1,len(lons)):
		# 	if np.abs(lons[i].value-lons[i-1].value) > 1:
		# 		sgn=np.sign(lons[i].value-lons[i-1].value)
		# 		lons[i:] += -sgn*adder

		
		#interpolation of manipulated and folded points 
		#flats = interp1d(tofs,lats)
		#flons = interp1d(tofs,lons)
		#tofs = np.linspace(0,T*3600*24*u.s,int(T*3600*24/dt))
		#lats = flats(tofs) * u.rad
		#lons = (lons.value) % (2*np.pi) * u.rad 

		return lons,lats






rho0 = rho0_earth.to(u.kg / u.km ** 3).value
H0 = H0_earth.to(u.km).value

@jit
def J2_and_drag(t0, state, k, J2, R, C_D, A_over_m, H0, rho0):
    return J2_perturbation(t0, state, k, J2, R) + atmospheric_drag_exponential(
        t0, state, k, R, C_D, A_over_m, H0, rho0
    )

