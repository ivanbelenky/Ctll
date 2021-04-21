import numpy as np
import matplotlib.pyplot as plt

from ..core import ctll


distannce_clean(ctll,T):
	lst = []
	sats = ctll.sats
	sat1 = sats[0]
	r1,v1 = sat1.rv(T,dt=100)

	for i in range(1,len(sats)):
		r2,v2 = sat[i].rv(T,dt=100)
		for r in r2:
			lst.append(Distance(r1,r,sat1.id,sat[i].id,sat.attractor))

	return lst

#this is nasty
def distances_nasty(ctll,T):
	lst = []
	k = 0
	sats = ctll.sats
	for i in range(len(sats)-1):
		k += 1
		for j in range(k,len(sats)):
			for r1,v in sats[i].rv(T,dt=100):
				for r2,v in sat[k].rv(T,dt=100):
					lst.append(Distance(r1,r2,sats[i].id,sats[k].id, sats[i].attractor))




class Distance(object):

	def __init__(self,r1,r2,id_1=None,id_2=None,R=None):
		self._distance = np.sqrt(np.sum((r1-r2)**2))
		if id_1 and id_2:
			self._id = (id_1, id_2)
		if R:
			self._R = R

		self.view = self._inview()



		@property
		def distance(self):
			return self._distance
		
		@property
		def view(self):
			return self._view
		
		@view.setter
		def view (self,value):
			self._view = value


	def _inview(self):
		try:
			sphere = Sphere(self.R)
			line = Line(r1,r2)

			if line.intersects(sphere):
				return False
			else:
				return True

		except AttributeError as e:
			return True


