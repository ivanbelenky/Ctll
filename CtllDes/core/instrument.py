# from ..requests import coverage as cob

class Instrument():
	def __init__(self):
		pass

	def coverage(self):
		raise NotImplemented 

	def communications(self):
		raise NotImplemented	



class Camera(Instrument):
	
	
	def __init__(self,pixels,resol):
		super().__init__()
		self._pixels = pixels
		self._resol = resol
	
	@cached_property
	def FOV(self):
		#TODO: calculate FOV from pixels and resol (maybe Ive missed an arguemnt)
		return 0               
	
	def coverage(self):
		return symmetric(self.FOV)
