from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit, propagation, states
from poliastro.twobody.states import ClassicalState as CS 
from poliastro.constants import J2000
from poliastro.frames import Planes
from astropy.time import Time, TimeDelta


import numpy as np
from astropy import time, units as u


from .satellite import Sat, SAT_ST 
from .specs import Specifications

try:
    from functools import cached_property  # type: ignore
except ImportError:
    from cached_property import cached_property  # type: ignore



class Ctll(object):

	def __init__(
		self,
		states,
		statuss=None,
		specs=None,
		instrumentss=None,
		pattern=None,
		epoch=J2000,		
	):

		"""Constructor.

		Parameters
		----------
		cops : list 
			list of poliastro.twobody.states BaseState objects childs
		states : list
			list of ~poliastro.twobody.states.BaseState childs
		statuss : list, optional
			list of SAT_STATUS strings, one for each sat 
		specs : list, optional
			list of satellite spec objects 
		epoch : ~astropy.time.Time, optional
		    Epoch, default to J2000.
		instrumentss : list,optional
			list of lists of instruments, sublist for each sat
		pattern : string 
			Pattern Type for constellation, pattern must be in PAT

		"""
		self._states = states
		self.statuss = statuss if statuss else []
		self.specs = specs if specs else []

		self._epoch = epoch

		self.instrumentss = instrumentss if instrumentss else []

		self._pattern = [(pattern, self.N)] if pattern \
		else [(PAT['NP'], self.N)]


		self._sats = self._set_sats()


	@property
	def states(self):
		return self._states
	
	@states.setter
	def states(self,statess):
		if not isinstance(statess,list):
			raise Exception("States must be a list")  
		elif len(statess) == 0:
			raise Exception("No state provided")
		self._statess = statess

	@property
	def statuss(self):
		return self._statuss

	
	#TODO: implement warning control, check if this lists are valid
	#this is just a patch to initialize without giving instruments and sepcs
	#lists. Maybe a way is to incorporate in the Instruments module and the specs
	#and status one, a generator for N satellites with one state.

	@statuss.setter
	def statuss(self,arg):
		if not isinstance(arg,list):
			raise Exception("Satellite status must be a list")
		
		N_statuss = len(arg)

		if N_statuss < self.N and N_statuss != 0:
			arg.append([SAT_ST["On"] for _ in range(self.N-N_statuss)])
			print("Not enought status, for remainining status ONLINE was used")
			self._statuss = arg
		elif N_statuss > self.N:
			print("Too many status, list has been sliced")
			self._statuss = arg[:self.N] 
		elif N_statuss == 0:
			self._statuss = [SAT_ST["On"] for _ in range(self.N)]
		else:
			self._statuss = arg

	@property
	def specs(self):
		return self._specs
	
	@specs.setter
	def specs(self,specs):
		if not isinstance(specs,list):
			raise Exception("Satellite specifications must be a list")
		
		N_specs = len(specs)

		if N_specs < self.N and N_specs > 0:
			specs.append([Specifications() for _ in range(self.N-N_specs)])
			print("Not enought specifications, for remainining specs.default was used")
			self._specs = specs
		elif N_specs > self.N:
			print("Too many specifications, list has been sliced")
			self._specs = specs[:self.N] 
		elif N_specs == 0:
			self._specs = [Specifications() for _ in range(self.N)]
		else:
			self._specs = specs

	@property
	def instrumentss(self):
		return self._instrumentss
	
	@instrumentss.setter
	def instrumentss(self,instrumentss):
		if not isinstance(instrumentss,list):
			raise Exception("Satellites instruments must be a list")
		
		N_instr = len(instrumentss) 
		
		if N_instr == 0:
			self._instrumentss = [ [] for _ in range(self.N)] 
		elif N_instr < self.N:
			print("Not enoough instruments")
			instrumentss.append([[] for _ in range(self.N-N_instr)]) 
		elif N_instr > self.N:
			print("Too many instruments")
			self._instrumentss = instrumentss[:self.N]

	@property
	def sats(self):
		return self._sats
	
	@property
	def SatsId(self):
		return [sat.id for sat in self.sats]
	
	@cached_property
	def N(self):
		"""Number of sats in the constellation """
		return len(self.states)

	@property
	def pattern(self):
		return self._pattern
	

	@property
	def epoch(self):
		"""Epoch of the orbit. """
		return self._epoch


	@classmethod
	def from_WalkerDelta(
		cls,
		T,
		P,
		F,
		p,
		ecc,
		inc,
		argp,
		raans_offset = 0,
		statuss=None,
		specs=None,
		instrumentss=None,
		epoch=J2000,
		plane=Planes.EARTH_EQUATOR,
		attractor=Earth

	):
		"""Returns Ctll from TPF, N, and classical orbit
		parameters.

		This method creates a WalkerDelta constellation type.
		IT sets the raan and nu classical parameters from  T/P/F
		Note that you must specify classical parameters to use this
		function directly. Is encouraged you use the IO module to
		manage your requests, the IO uses mainly this clsmethod
		whenever the Walker Pattern is selected.

		Parameters
		----------
		T : int
		    # of total satellites in the constellation
		P : int
			# of planes
		F : int 
			Nodal phase angle between planes, must be in [0,T-1]
		a : ~astropy.units.Quantity
		    Semi-major axis.
		ecc : ~astropy.units.Quantity
		    Eccentricity.
		inc : ~astropy.units.Quantity
		    Inclination
		argp : ~astropy.units.Quantity
		    Argument of the pericenter.
		statuss : list
			list of SAT_STATUS strings, one for each sat  
		epoch : ~astropy.time.Time, optional
		    Epoch, default to J2000.
		plane : ~poliastro.frames.Planes
		    Fundamental plane of the frame.
		
		Note
		----	
		raan and nu are omitted since TPF determines them

		"""	
		S = T/P
		
		raans = [raans_offset+int(j/S)*360*u.deg/P for j in range(T)]
		
		nus = [(j%S)*360*u.deg/S+int(j/S)*360*u.deg*F/T 
		for j in range(T)
		]

		states = [CS(attractor,p,ecc,inc,raans[j],argp,nus[j],plane) 
		for j in range(T)]
		

		return cls(
			states,
			statuss,
			specs,
			instrumentss,		
			pattern=PAT['WD']+f' T/P/F = {T}/{P}/{F}',
			epoch=epoch
		)


	def _set_sats(self):
		"""Returns sats attribute: list of sats"""

		return [Sat(st,status,spec,instr) for st,status,spec,instr
		 in zip(self.states,self.statuss,self.specs,self.instrumentss)]



	def rv(self,T,dt=1.,method=propagation.cowell,**kwargs):
		""" Propagates orbit for online satellites,
		T days of flight.

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

		rrvv = [sat.rv(T,dt,method,**kwargs) for sat in self.sats
		if sat.status is SAT_ST["On"]]	
		return rrvv

	def ssps(self,T,dt=1.,method=propagation.cowell,**kwargs):
		""" Get subsatellite points for online Satellites,
		T days of flight.

		Parameters
		----------
		T : float
			Desired time of propagation in days
		dt : float
			Size of intervals in seconds
		
		Returns
		-------
		sspss : list
			list of tuples (lat,long), lat and lon are Quantity 
			objects array, size = floor(T*24*3600/dt)
		
		"""

		sspss = [sat.ssps(T,dt,method,**kwargs) for sat in self.sats
		if sat.status is SAT_ST["On"]]	
		
		return sspss


	def getSat(self,uuid):
		"""Get satellite by uuid
		
		Parameters
		----------
		uuid : UUID
			id of desired satellite

		"""

		for sat in self.sats:
			if sat.id == uuid:
				return sat
		raise Exception("Satellite not found")


	def UpdateStatus(self,newStatuss):
		"""Updates ctll sats status.
		
		Parameters
		----------
		newStatuss : list
			list of CtllDes.SAT.SAT_ST values
		"""
		
		self.statuss = newStatuss 
		for sat,status in zip(self.sats,newStatuss):
			sat.UpdateStatus(status)

	def UpdateSpecs(self,newSpecs):
		"""Updates ctll sats specifications.
		
		Parameters
		----------
		newSpecs : list
			list of CtllDes.core.specs.Specifications 
		"""
		
		self.statuss = newSpecs 
		for sat,spec in zip(self.sats,newSpecs):
			sat.UpdateSpec(spec)

	def UpdateInstruments(self,newInstr):
		"""Updates ctll sats Instruments.
		
		Parameters
		----------
		newInstr : list
			list of CtllDes.core.instrument.Instrument
		"""
		
		self.instrumentss = newInstr 
		for sat,instr in zip(self.sats,newInstr):
			sat.UpdateInstruments(instr)


	def getSatsId(self):
		"""Returns satellites ids list in order."""
		return [sat.id for sat in self.sats]

	def getOnlineSatsId(self):
		"""Returns online satellites ids list in order.
		"""
		return [sat.id for sat in self.sats 
		if sat.status is SAT_ST["On"]]

	def getOfflineSatsId(self):
		"""Returns offline satellites ids list in order.
		"""
		return [sat.id for sat in self.sats 
		if sat.status is SAT_ST["Off"]]


	def Info(self,v=False):
		"""Prints out Constellation info.
		
		Parameters
		----------
		v : boolean
			verbose option, display all SATs info
		
		"""

		if not v:
			for patt in self.pattern:
				print(f"\n{patt[1]} satellites within {patt[0]}")

		else:
			offset = 0 
			for j in range(len(self.pattern)):	
				print(f"\n{self.pattern[j][1]} satellites in "+
					f"{self.pattern[j][0]}\n")
				for i in range(self.pattern[j][1]):
					self.sats[offset+i].Info()
				offset += self.pattern[j][1]
				





PAT = {
	'NP':'no Pattern',
	'WD':'Walker Delta Pattern'
}
