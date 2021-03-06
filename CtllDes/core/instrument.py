# from ..requests import coverage as cob

import numpy as np

from astropy import units as u 

from CtllDes.requests.coverage import symmetric, push_broom, _push_broom, symmetric_with_roll

import uuid


class Instrument(object):
	def __init__(self):
		self._id = uuid.uuid4()
		
	def coverage(self,lons,lats,r,v,target,R):
		"""Coverage functions is associated with the coverage module.
		Any overwrited child method coverage must accept the specified parameters
		and return a list or iterable with 1 or 0, in view or not respectively.

		Parameters
		----------
		lons : ~astropy.units.quantity.Quantity 
			array of longitudes as they come from the ssps method.
		lats : ~astropy.units.quantity.Quantity 
			array of latittudes as they come from the ssps method
		r : ~astropy.units.quantity.Quantity
			satellite's positions
		v : ~astropy.units.quantity.Quantity
			satellite's velocities
		target : ~CtllDes.targets.targets.Target
			desired target of coverage analysis
		R : ~astropy.units.quantity.Quantity 
			attractor mean radius 

		Returns
		-------
		cov : Iterable
			elements from iterable must be 1 or 0 indicating if target is in view 
			or not. 

		"""
		raise NotImplemented 

	def communications(self):
		raise NotImplemented	

	@property
	def id(self):
		return self._id
	


class Camera(Instrument):
	
	def __init__(self, f_l, s_w):
		super().__init__()
		self.f_l = f_l
		self.s_w = s_w
		
		self.FOV = np.arctan(self.s_w/2/self.f_l)*u.rad 
			
	def coverage(self,lons,lats,r,v,target,R):
		return symmetric(self.FOV,lons,lats,r,v,target,R)


class RollCamera(Instrument):
	def __init__(self, f_l, s_w, roll_angle):
		"""Constructor for RollCamera.

		Parameters
		----------
		f_l : ~astropy.units.quantity.Quantity
		    focal length
		s_W : ~astropy.units.quantity.Quantity
		    sensor width
		roll_angle : ~astropy.units.quantity.Quantity
		    maximum rolling angle
		"""

		super().__init__()
		self.f_l = f_l
		self.s_w = s_w 
		self.FOV = np.arctan(self.s_w/2/self.f_l)
		self.roll = roll_angle.to(u.rad)
        
	def coverage(self, lons, lats, r, v, target, R):
		return symmetric_with_roll(self.FOV,
									lons,
									lats,
									r,
									v,
									target,
									R,
									roll_angle = self.roll)




class GodInstrument(Instrument):
    def __init__(self):
        super().__init__()

    def coverage(self, lons, lats, r, v, target, R):
        return [1 for _ in range(len(r))]



class PushBroom(Instrument):
	""" Constructor for PushBroom.

	Parameters
	----------
	f_l : ~astropy.units.quantity.Quantity
	    focal length
	s_W : ~astropy.units.quantity.Quantity
	    sensor width
	"""
	
	def __init__(self, fl, sw):
		super().__init__()
		self.f_l = fl
		self.s_w = sw 
		self.FOV = np.arctan(self.s_w/2/self.f_l)

	def coverage(self, lons, lats, r, v, target, R):
		return push_broom(self.FOV, lons, lats, r, target, R)


class RollPushBroom(Instrument):
	""" Constructor for PushBroom.

	Parameters
	----------
	f_l : ~astropy.units.quantity.Quantity
	    focal length
	s_W : ~astropy.units.quantity.Quantity
	    sensor width
	roll_angle : ~astropy.units.quantity.Quantity
	    maximum rolling angle
	"""
	
	def __init__(self, fl, sw, roll_angle):
		super().__init__()
		self.f_l = fl
		self.s_w = sw
		self.roll_angle = roll_angle 
		self.FOV = np.arctan(self.s_w/2/self.f_l)

	def coverage(self, lons, lats, r, v, target, R):
		return push_broom(self.FOV + self.roll_angle, lons, lats, r, target, R)



