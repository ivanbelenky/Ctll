from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit, propagation, states
from states import ClassicalState as CS 

import numpy as np
from astropy import time, units as u

from .SAT import Sat 




class Ctll(object):

	def __init__(
		self,
		states,
		statuss,
		specs,
		epoch=J200,
		instrumentss=None,
		pattern=None,		
	):

		"""Constructor.
	
		Parameters
		----------
		cops : list 
			list of dictionaries with clasiccal orbit parameters
		statuss : list
			list of SAT_STATUS strings, one for each sat 
		specs : list
			list of satellite spec objects 
		epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.
		instrumentss : list,optional
			list of lists of instruments, sublist for each sat
		pattern : string 
			Pattern Type for constellation, pattern must be in PAT
		
		"""
		self._states = states
		self._statuss = statuss
		
		self._epoch = epoch

		self._instrumentss = instrumentss

		self._pattern = [(pattern, self.N)] if pattern
		else [(PAT['NP'], self.N)]


		self._sats = self._set_sats()



	@property
	def statuss(self):
		return self._statuss

	@property
	def sats(self):
		return self._sats
	
	@property
	def SatsId(self):
		return [sat.id for sat in self.sats]
	
	@property
	def N(self):
		"""Number of sats in the constellation """
		return len(self._states)

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
        statuss,
        epoch=J2000,
        plane=Planes.EARTH_EQUATOR,
        attractor=Earth,
        instruments=None
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
   		
   		raan and nu are omitted since the TPF determines them

		"""	
		S = T/P
		
		raans = [int(j/S)*360*u.deg/P for j in range(T)]
		
		nus = [(j%S)*360*u.deg/S+int(j/S)*360*u.deg*F/T 
		for j in range(T)
		]

		states = [CS(attractor,p,ecc,raans[j],argp,nus[j],plane) 
		for j in range(T)]
		

		return cls(
			states,
			statuss,
			instrumentss,
			epoch,
			pattern=PAT['WD']+f' T/P/F = {T}/{P}/{F}',
		)


	def _set_sats(self):
		"""Returns sats attribute: list of sats"""

		return [Sat(st,instr,status) for st,instr,status
		 in zip(states,instrumentss,statuss)]



	def Propagate(self,T,dt=100):
		""" Propagates orbit for all Online satellites.

		Parameters
		----------
		T : int
			Time desired for propagation
		dt : float
			Size of intervals in seconds
		Returns
		-------
		propagation : list
			list of tuples (rr,vv), rr and vv are Quantity 
			objects array, size=T*24*3660/dt

		Propagate returns list of tuples. Each tuple contains the 
		coordinates of the   

		Whenever needed outside the class the classic orbit elements can be
		obtained using built-in poliastro.core.elements.rv2coe.
		The propagation results are the value in km, but this
		not being a Quantity.
	
		TODO: Find best way to 
		add	perturbations.
		"""

		propagation = [sat.Propagate(T,dt) for sat in self.sats
		if sat.status is SAT_STATUS_ONLINE]	
		
		return propagation

	def UpdateStatus(self,newStatuss):
		"""Updates ctll sats status."""

		self.statuss = newStatuss
		for sat,status in zip(self.sats,newStatuss):
			sat.UpdateStatus(status)


	def getOnlineSatsId(self):
		"""Returns online sats ids list in order of satellites
		Online
		"""
		return [sat.id for sat in self.sats 
		if sat.status is SAT_STATUS_ONLINE]



	def Info(self,v=False):
		"""Prints out Constellation info.
		
		Parameters
		----------
		v : boolean
			verbose option to display all SATs info
		
		TODO: nicer way for 
		"""

		if not v:
			for patt in self.pattern:
				print(f"{patt[1]} satellites within {patt[0]}")

		for j in range(len(self.pattern)):	
			print(f"{self.pattern[j][1]} satellites in"+
				f"{self.pattern[j][0]}\n")
			for i in range(patt[1]):
				self.sats[j+i].Info()





SAT_STATUS_ONLINE = 'Status:Online'
SAT_STATUS_OFFLINE = 'Status:Offline'
PAT={
	'NP':'no Pattern',
	'WD':'Walker Delta Pattern'
	}
