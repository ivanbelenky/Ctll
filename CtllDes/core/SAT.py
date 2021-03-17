from poliastro.bodies import Earth
from poliastro.twobody import Orbit, propagation

import uuid


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
	









	@classmethod 
	def from_HelioSync(cls,):

		return cls()

	@classmethod
	def GeoSync():

		return cls()

	@classmethod
	def GeoHelioSync(cls, ):

		return cls()


	def Propagate(self,T,dt=1.):
	"""Propagates orbit specified amount of days"""


		#TODO: check cached properties, see best way to fit
		#perturbations
		
		tofs = TimeDelta(np.linspace(0, T*24*3600*u.s,
		 num=int(T*24*3600/dt)))

		#TODO find propagation best way

		rr = coords.xyz.T.to(u.km).value
		vv = coords.differentials["s"].d_xyz.T.to(u.km / u.s).value

		return [rr,vv]



	def UpdateStatus(self,newStatus):
		"""Updates satellite's status
		
		Parameters
		----------
		newStatus : string
			SAT_STATUS strings
		"""

		self.status=newStatus


	def Info(self):
		print(f"id:{self.id}\t[a ecc inc raan argp nu ]:"+
			f"[ {self.a} {self.ecc} {self.inc}"+
			f" {self.raan} {self.argp} {self.nu}]"+
			f"\tstatus:{self.status}")
		 

