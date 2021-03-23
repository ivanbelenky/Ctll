from poliastro.bodies import Earth
from poliastro.twobody import Orbit, propagation

import uuid

PROPAGATOR_DT = 100

class Sat(object):

	def __init__(self,
		state,	
		status,
		spec,
		epoch=J2000,
		instruments=None,

	):
    	#TODO masa de satélite CD and whatnot.

		"""Constructor.

		Parameters
		----------
		state : ~poliastro.twobody.states.ClassicalState
			State for satellite orbit
		status : string
			SAT_STATUS string
		spec : ~Ctlldes.core.Spec
			specifications
		epoch : ~astropy.time.Time, optional
		    Epoch, default to J2000.
		instruments : list
			List of Instrument objects

		"""
		self._state = state
		self._status = status
		self._spec = spec

		self._epoch = epoch

		self._instruments = instruments


		self._id = uuid.uuid4()
		self._orbit = Orbit(state,epoch)
		self._Propagator = Propagator(self.orbit)

	@property
	def state(self):
		return self._state


	@property
	def status(self):
		return self._status
		

	@property
	def instruments(self):
		return self._instruments

	@property
	def id(self):
		return self._id

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
	def CovInstruments(self):
		covInstr = [instr for instr in self.instruments
		if "_coverage" in set(dir(instr))]
		return self.covInstr
	
	#TODO: property setters for Instruments objects 



	# @classmethod 
	# def from_HelioSync(cls,):

	# 	return cls()

	# @classmethod
	# def GeoSync():

	# 	return cls()

	# @classmethod
	# def GeoHelioSync(cls, ):

	# 	return cls()


	def rv(self,T,dt=1.,method=propagation.cowell,**kwargs):
	"""Propagates orbit specified amount of days.
	

	"""
		rr,vv = self.Propagator.get_rv(self.orbit,T,dt,
			method=method,**kwargs)

		return rr,vv

	def ssps(self,T,dt=1.,method=propagation.cowell,**kwargs):
		"""Get subsatellite points for satellite
		
		"""

		tofs = np.linspace(0,T*3600*24*u.s,int(T*3600*24/dt))
		rr,vv = self.rv(T,dt,method,**kwargs)
		rs,lats,lons = trigsf.c2s(
			rr[0,:].value,
			rr[1,:].value,
			rr[2,:].value
		)

		w = self.attractor.angular_velocity
		lons = np.array([lon-t*w for lon in lons])

		return lats,lons


	def UpdateInstruments(self,instruments):

		self.instruments = instruments

	def UpdateStatus(self,newStatus):
		"""Updates satellite's status
		
		Parameters
		----------
		newStatus : string
			SAT_STATUS string
		"""

		self.status=newStatus


	def Info(self):
		print(f"id:{self.id}\t[a ecc inc raan argp nu ]:"+
			f"[ {self.a} {self.ecc} {self.inc}"+
			f" {self.raan} {self.argp} {self.nu}]"+
			f"\tstatus:{self.status}")
		 





class Propagator(object):

	def __init__(
		self,
		orbit,
	 	T = 1,
	 	method=propagation.cowell,	   
		**kwargs
	):
		
		"""Propagator object builder. Build coordinates"""

		self._orbit = orbit 
		self._T = T 
		self._kwargs = kwargs
		self._method = method
		#TODO: check if this is ok unpacking.

		self._tofs = None
		self._coords = self._setcoords()
		

		#hardsetted to optimize performance on interpolation
		self._propagatorDT = PROPAGATOR_DT 



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
		if not isinstance(T,float): 	
			raise Exception("T has to be float or int")
		if T <= 0: raise Exception("T must be positive and != 0")
		self._T = T
	
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
		 num=int(self.T*24*3600/self._propagatorDT)))
		
		self.coords = propagation.propagate(self.orbit,
			self.tofs,method=self.method,**self.kwargs)


	def get_rv(self,T,dt=1.,*,method=propagation.cowell,**kwargs):

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



		if not isinstance(T,float): 
			raise Exception("T has to be float or int")
		if T <= 0: raise Exception("T must be positive and != 0")


		if not isinstance(dt,float): raise Exception("dt has to be float")
		if dt <= 0: raise Exception("dt must be positive and !=0")	



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

		tofs = TimeDelta(np.linspace(0,T*24*3600*u.s,
		 num=int(T*24*3600/dt)))


		ephemerides = ephem.Ephem(self.coords,self.tofs,
			Planes.EARTH_EQUATOR)
		
		return ephemerides.rv(tofs)
	



