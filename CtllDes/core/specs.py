

DefaultSpec = {
	"m":1,
	"Cd":1,
	"A_over_m":1
}

class Specifications(object):
	"""Satellite Specifications for physical information about satellites"""
	def __init__(
		self,
		m=DefaultSpec["m"],
		Cd=DefaultSpec["Cd"],
		A_over_m=DefaultSpec["A_over_m"],
		#args for self calculation of drag coefficient
	):

		self._m = m
		self._Cd = Cd
		self._A_over_m = A_over_m
	
	@property
	def m(self):
		return self._m
	
	@property
	def Cd(self):
		return self._Cd
	
	@property
	def A_over_m(self):
		return self._A_over_m
	





