# from ..requests import coverage as cob

import numpy as np

from astropy import units as u 

import CtllDes.requests.coverage as cov

import uuid


class Instrument(object):
	def __init__(self):
		self._id = uuid.uuid4()
		
	def coverage(self):
		raise NotImplemented 

	def communications(self):
		raise NotImplemented	

	@property
	def id(self):
		return self._id
	

#distancia focal tama;o de pixel, cantidad de pixels
#P/F 
class Camera(Instrument):
	
	def __init__(self,pixels,resol):
		super().__init__()
		self._pixels = pixels
		self._resol = resol
		self._FOV = (np.pi/10)*u.rad

	@property
	def FOV(self):

		#TODO: calculate FOV from pixels and resol (maybe Ive missed an arguemnt)
		#FOV must be set with units radians, units of astropy
		#may fix this but is nice to be consistent.

		return self._FOV               
	
	def coverage(self):
		return cov.symmetric(self.FOV)
