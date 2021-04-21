
from astropy import units as u


DefaultSpec = {
	"m":1*u.kg,
	"Cd":1*u.one,
	"A":1*u.m*u.m
}

class Specifications(object):
	"""Satellite Specifications for physical information about satellites"""
	def __init__(
		self,
		m=DefaultSpec["m"],
		Cd=DefaultSpec["Cd"],
		A=DefaultSpec["A"],
		#args for self calculation of drag coefficient
	):

		self._m = m
		self._Cd = Cd
		self._A = A
		self._A_over_m = self.A/self.m
	

	@property
	def m(self):
		return self._m
	
	@property
	def Cd(self):
		return self._Cd
	
	@property
	def A(self):
		return self._A

	@property
	def A_over_m(self):
		return self._A_over_m

	
	

	def __str__(self):
		return f" m:{self.m}, Cd:{self.Cd}, A:{self.A}"





