from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit, propagation, states
from poliastro.twobody.states import ClassicalState as CS 
from poliastro.constants import J2000
from poliastro.frames import Planes
from astropy.time import Time, TimeDelta

import pandas as pd
import numpy as np
from astropy import time, units as u
import time
import copy



from .satellite import Sat, SAT_ST 
from .specs import Specifications
from .instrument import Instrument
try:
    from functools import cached_property  # type: ignore
except ImportError:
    from cached_property import cached_property  # type: ignore



class Ctll(object):

	def __init__(
		self,
		states,
		status=None,
		specs=None,
		instruments=None,
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
		status : list, optional
			list of SAT_STATUS strings, one for each sat 
		specs : list, optional
			list of satellite spec objects 
		epoch : ~astropy.time.Time, optional
		    Epoch, default to J2000.
		instruments : list,optional
			list of lists of instruments, sublist for each sat
		pattern : string 
			Pattern Type for constellation, pattern must be in PAT

		"""
		self._states = states
		self.status = status if status else []
		self.specs = specs if specs else []
		self.instruments = instruments if instruments else []
		self._pattern = pattern if pattern \
		else [{"PAT":PAT['NP'], "N":self.N}]

		self._epoch = epoch

		self._sats = self._set_sats()


	@property
	def states(self):
		return self._states
	
	@states.setter
	def states(self,statess):
		#TODO: chheckk that all states attractor is the same one.
		if not isinstance(statess,list):
			raise Exception("States must be a list")  
		elif len(statess) == 0:
			raise Exception("No state provided")
		self._statess = statess

	@property
	def status(self):
		try:
			a = [sat.status for sat in self.sats]
			return a
		except:
			return self._status

	
	#lists. Maybe a way is to incorporate in the Instruments module and the specs
	#and status one, a generator for N satellites with one state.

	@status.setter
	def status(self,arg):
		if not isinstance(arg,list):
			raise Exception("Satellite status must be a list")
		
		N_status = len(arg)

		if N_status < self.N and N_status != 0:
			arg += [SAT_ST["On"] for _ in range(self.N-N_status)]
			print("Not enought status, for remainining status ONLINE was used")
			self._status = arg
		elif N_status > self.N:
			print("Too many status, list has been sliced")
			self._status = arg[:self.N] 
		elif N_status == 0:
			self._status = [SAT_ST["On"] for _ in range(self.N)]
		else:
			self._status = arg

	@property
	def specs(self):	
		try:
			a = [sat.spec for sat in self.sats]
			return a
		except:
			return self._specs
		
	
	@specs.setter
	def specs(self,specs):
		if not isinstance(specs,list):
			raise Exception("Satellite specifications must be a list")
		
		N_specs = len(specs)

		if N_specs < self.N and N_specs > 0:
			specs += [Specifications() for _ in range(self.N-N_specs)]
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
	def instruments(self):
		try:
			a = [sat.instruments for sat in self.sats]
			return a
		except:
			return self._instruments
		
		
	@instruments.setter
	def instruments(self,instruments):
		if not isinstance(instruments,list) and not isinstance(instruments,Instrument):
			raise Exception("Satellites instruments must be a list or Instrument object")
		
		if isinstance(instruments,Instrument):
			instruments = [ instruments for _ in range(self.N)]
			#TODO: check the implications of this hackerino, nasty it is.

		N_instr = len(instruments) 
		
		if N_instr == 0:
			self._instruments = [ [] for _ in range(self.N)] 
		elif N_instr < self.N:
			print(f"Not enough instruments, no instruments for last {self.N-N_instr} satellites")
			instruments += [[] for _ in range(self.N-N_instr)]
			self._instruments = instruments 
		elif N_instr > self.N:
			print("Too many instruments")
			self._instruments = instruments[:self.N]
		else:
			self._instruments = instruments

	@property
	def sats(self):
		return self._sats

	@cached_property
	def N(self):
		"""Number of sats in the constellation """
		return len(self.states)

	@property
	def sats_id(self):
		return [sat.id for sat in self.sats]
	
	@property
	def online_id(self):
		return [sat.id for sat in self.sats 
		if sat.status is SAT_ST["On"]]
	
	@property
	def offline_id(self):
		return [sat.id for sat in self.sats 
		if sat.status is SAT_ST["Off"]]


	@property
	def pattern(self):
		return self._pattern
	

	@property
	def epoch(self):
		"""Epoch of the orbit. """
		return self._epoch


	@classmethod
	def from_sats(cls,sats,pattern=None,epoch=J2000):
		"""Build constellation from satellite list or unique satellite.
		
		sats : ~CtllDes.core.satellite.Sat list
			satellites
		pattern : string
			Satellite pattern
		"""
		if isinstance(sats,Sat):
			return cls([sats.state],
			status=[sats.status],
			specs=[sats.spec],
			instruments=sats.instruments,
			pattern=pattern,
			epoch=epoch)


		states = [sat.state for sat in sats]
		status = [sat.status for sat in sats]
		specs = [sat.spec for sat in sats]
		instruments = [sat.instruments for sat in sats]
		
		return cls(states,
			status=status,
			specs=specs,
			instruments=instruments,
			pattern=pattern,
			epoch=J2000)

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
		raan_offset = 0*u.deg,
		nu_offset = 0*u.deg,
		status=None,
		specs=None,
		instruments=None,
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
		argp : ~astropy.units.Quantity, optional
		    Argument of the pericenter.
		raan_offset : ~astropy.units.Quantity
			right ascension of the ascending node offset
		nu_offset : ~astropy.units.Quantity, optional
			True anomaly offset
		status : list
			list of SAT_STATUS strings, one for each sat  
		epoch : ~astropy.time.Time, optional
		    Epoch, default to J2000.
		plane : ~poliastro.frames.Planes
		    Fundamental plane of the frame.
		
		Note
		----	
		raan and nu are omitted since TPF determines them unless
		the offset is granted. This values force the ctll to pass
		over those points at least one time, at the beggining of 
		the orbit. 

		"""	
		
		S = T/P
		
		raans = [((raan_offset.to(u.deg).value+int(j/S)*360/P)%360)*u.deg for j in range(T)]
		
		nus = [( ((j%S)*360/S+int(j/S)*360*F/T +nu_offset.to(u.deg).value) % 360)*u.deg 
		for j in range(T)
		]

		nu_ = []
		for nu in nus:
			if nu > 180*u.deg:
				nu_.append(nu-360*u.deg)
			else:
				nu_.append(nu)

		nus = nu_

		states = [CS(attractor,p,ecc,inc,raans[j],argp,nus[j],plane) 
		for j in range(T)]
		

		return cls(
			states,
			status,
			specs,
			instruments,		
			pattern=[{"PAT":PAT['WD']+f' T/P/F = {T}/{P}/{F}',"N":T}],
			epoch=epoch
		)


	def _set_sats(self):
		"""Returns sats attribute: list of sats"""

		return [Sat(st,status,spec,instr) for st,status,instr,spec
		 in zip(self.states,self.status,self.instruments,self.specs)]



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

	def ssps(self,T,dt=1.,method=propagation.cowell, lon_offset = 0, **kwargs):
		""" Get subsatellite points for online Satellites,
		T days of flight.

		parameters
		----------
		T : float
			Desired time of propagation in days
		dt : float
			Size of intervals in seconds
		
		Returns
		-------
		sspss : list
			list of tuples (lon,lat), lat and lon are Quantity 
			objects array, size = floor(T*24*3600/dt)
		
		"""

		sspss = [sat.ssps(T,dt,method,lon_offset=lon_offset,**kwargs) for sat in self.sats
		if sat.status is SAT_ST["On"]]	
		
		return sspss

	def __str__(self):
		
		string = ""
		for pat in self.pattern:
			string += f"{pat['N']} satellites within {pat['PAT']}\n"
		
		return string


	def __add__(self, other):
		if isinstance(other, __class__):
			new_pattern = other.pattern + self.pattern
			new_states = other.states + self.states 
			new_status = other.status + self.status
			new_specs = other.specs + self.specs
			new_instruments = other.instruments + self.instruments
			epoch = self.epoch

			return Ctll(new_states, new_status, new_specs, new_instruments,
			new_pattern,epoch)

		else:
			raise TypeError("Invalid types")

	#TODO: check if it works or if is neccesary
	def get_sat(self,uuid):
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

	#TODO: check if it works or if is neccesary
	def update_status(self,newStatuss):
		"""Updates ctll sats status.
		
		Parameters
		----------
		newStatuss : list
			list of CtllDes.SAT.SAT_ST values
		"""
		
		self.status = newStatuss 
		for sat,status in zip(self.sats,newStatuss):
			sat.UpdateStatus(status)
	

	#TODO: check if it works or if is neccesary 
	def update_specs(self,new_specs):
		"""Updates ctll sats specifications.
		
		Parameters
		----------
		new_specs : list
			list of CtllDes.core.specs.Specifications 
		"""
		
		self.status = new_specs 
		for sat,spec in zip(self.sats,new_specs):
			sat.update_spec(spec)


	#TODO: check if it works or if is neccesary
	def update_instruments(self,newInstr,f=False):
		"""Updates ctll sats Instruments.
		
		Parameters
		----------
		newInstr : list
			list of CtllDes.core.instrument.Instrument
		"""
		
		self.instruments = newInstr 
		for sat,instr in zip(self.sats,newInstr):
			sat.update_instruments(instr,f)


	def info(self,v=False):
		"""Prints out Constellation info.
		
		Parameters
		----------
		v : boolean
			verbose option, display all SATs info
		
		"""

		if not v:
			for patt in self.pattern:
				print(f"\n{patt['N']} satellites within {patt['PAT']}")

		else:
			offset = 0 
			for j in range(len(self.pattern)):	
				print(f"\n{self.pattern[j]['N']} satellites in "+
					f"{self.pattern[j]['PAT']}\n")
				for i in range(self.pattern[j]['N']):
					self.sats[offset+i].info()
				offset += self.pattern[j]['N']
				





PAT = {
	'NP':'no Pattern',
	'WD':'Walker Delta Pattern'
}
