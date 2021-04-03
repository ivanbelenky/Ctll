from ..core import ctll, Instrument, sat
import collections.abc
from collections.abc import Iterable


try:
    from functools import cached_property  # type: ignore
except ImportError:
    from cached_property import cached_property  # type: ignore

TOL = 1**(-10)

def isCovered(ssps,r,target,R,coverage_method):
	"""The CoverageMethod, returns an arbitrary length tuple of
	comparing functions. 

	Parameters
	----------
	

	This functions must be such that their input
	are, the sub satellite points, the position and the target in question. 
	It returns a list the same length of subsatellite points  """

	outputs = np.array(coverage_method(ssps,r,target,R))
	#column-wise multiplication, checks all requirements.
	cov = [ np.prod(outputs[:,i]) for i in outputs.shape[0] ]
	return cov	


#list of built-in coverage methods

def symmetric(FOV):
	def _symmetric(ssps,r,target,R):
		"""TODO docstring"""
		radiis = np.sqrt(np.sum(a**2,axis=1))

		rho = np.arcsin(R/radiis)*u.rad
    	eps = np.arccos((np.sin(FOV))/(np.sin(rho)))*u.rad
    	lam = (np.pi/2)*u.rad - FOV - eps

		#TODO: missing math 

    	boolean = []
    	for _ in angles:
    		if _ <= lam:
    			boolean.append(1)
    		else:
    			boolean.append(0)
    	return boolean

    return _symmetric

def symmetric_disk(FOV_max,FOV_min):
	"""TODO docstring"""
	def _symmetric_disk(ssps,r,target,R)
		radiis = np.sqrt(np.sum(a**2,axis=1))

		rho = np.arcsin(R/radiis)*u.rad
    	eps = np.arccos((np.sin(FOV))/(np.sin(rho)))*u.rad
    	lam_min = (np.pi/2)*u.rad - FOV_min - eps
    	lam_max = (np.pi/2)*u.rad - FOV_max - eps
    	
		#TODO: missing math 

    	boolean = []
    	for _ in angles:
    		if  lam_min <= _ <= lam_max:
    			boolean.append(1)
    		else:
    			boolean.append(0)
    	
    	return boolean

    return _symmetric_disk




def symmetric_with_roll():
	pass



COVERAGE_COMPARING_METHODS = [symmetric, symmetric_disk, symmetric_with_roll]


class Coverages(collections.abc.Set):
		def __init__(self,covs,tag=None):
		self._covs = lst = list()
		self._tag = tag if tag else "No Tag"

		if isinstance(covs,Coverage):
			lst.append(covs)
		elif not isinstance(covs,Iterable):
			raise TypeError("covs must be Coverage, or iterable collection of Coverage objects")
		else: 	
			for cov in covs:
				if not isinstance(cov,Coverage):
					raise TypeError("covs must be a collection of Coverage objects") 
				if cov not in lst:
					lst.append()

		self._targets = {cov.target for cov in self.covs}

	@property
	def covs(self):
		return self._covs
	
	def __iter__(self):
		return iter(self.covs)

	def __contains__(self,value):
		return value in self.covs

	def __len__(self):
		len(self.covs)

	def __str__(self):
		return self.tag


	@cached_property
	def data(self):
		return self.to_data()
		#TODO: check if it is worth doing named tuple or other thing


	@classmethod
	def from_ctll(cls,ctll,targets,T,dt):
		"""Get coverages from constellation.
		
		Parameters
		----------
		ctll : CtllDes.core.ctll.Ctll
			CtllDes constellation object
		targets : CtllDes.targets.targets.Targets
			Desired targets to build coverages
		T : float | int
			Desired time of coverage in days
		dt : float | int, optional
			time of sampling in seconds

		Returns
		-------
		Coverages object containing only the coverage from
		instruments that have _coverage function implemented.
		
		TODO: provide duck typing or not for instruments in ctll?
		"""

		ctll_covs = []
		for sats in ctll.sats:
			try:
				sat_covs = Coverages.from_sat(sat,targets,T,dt)
			except Exception:
				pass
			else:
				ctll_covs += (sat_covs)

		if not len(Coverages): 
			raise Exception("Constellation has no Coverage Instruments")

		return Coverages(ctll_covs,tag=ctll.__str__())


	@classmethod
	def from_sat(sat, targets,T, dt=1.):
		"""Build list of coverage objects from satellite
	
		Parameters
		----------
		sat : ~CtllDes.core.sat
			sat object
		targets : ~CtllDes.targets.targets.Targets 
			Desired targets of coverage
		T : float
			Desired Time of analysis in days.
		dt : float | int, optional
			time of sampling in seconds

		Returns
		-------
		Coverages : list
			List of Coverage objects, one for each target and  
			coverage instrument.

		"""


		cov_instruments = [sat.instr for instr in sat.cov_instruments]
		if not len(cov_instruments):
			raise Exception("No coverage instruments found on" +
				f" satID : {sat.id}")	

		ssps = sat.ssps(T,dt)
		r,v = sat.rv(T,dt)

		sat_coverages = []
		for instr in cov_instruments:
			for target in targets:
				cov = isCovered(ssps,r,target,sat.attractor.R_mean,instr.coverage())
				sat_coverages.append(Coverage(cov,target,T,dt,instr.id))

		return sat_coverages
		

	
	def filter_by_target(self,target):
		"""Obtain coverages filtered by target
		
		Parameters
		----------
		target : CtllDes.targets.targets.Target		
		"""
		lst = [cov for cov in self.covs if cov.target == target]
		if not len(lst):
			print(f"{target} not found in coverages")
		else:
			return Coverages(lst,tag=f"{self.tag} filtered by target = {target}")
		
	def filter_by_instrument(self,instr_id):
		"""Obtain coverages filtered by instrument_id
		
		Parameters
		----------
		instr_id : UUID
			desired instrument id to filter out
		"""
		lst = [cov for cov in self.covs if cov.id == instr_id]
		if not len(lst):
			print(f"Instrument with ID : {instr_id} not found")
		else:
			return Coverages(lst,tag=f"{self.tag} filtered by target = {instr_id}")

	def collapse_instruments(self):
		"""Obtain coverages regardless of instrument.

		Returns
		-------
		collapsed : CtllDes.requests.coverage.Coverages
			Coverages object with Covs merged for target

		"""
		tgts = self.targets
		cov_tgt = [[] for tgt in tgts]
		for cov in self.covs:
			idx = tgts.index(cov.target)
			cov_tgt[idx].append(cov)

		collapsed = []	
		for covs in cov_tgt: 
			new_cov = covs[0]
			for i in range(1,len(covs)):
				new_cov = new_cov + covs[i]
			collapsed.append(new_cov)

		return Coverages(collapsed,tag=f"{self.tag} collapsed")



		

class Coverage(object):
	#TODO: docstring
	
	def __init__(self,
		cov,
		target,
		T,
		dt,
		instr_id=None,
	):
		self._cov = cov
		self._target = target
		self._instr_id = instr_id 
		self._T = T
		self._dt = dt


	def __add__(self,other):
		if self.target != other.target:
			raise ValueError("Targets must be equal")
		elif self.T-other.T<TOL or self._dt-other.dt<TOL:
			raise ValueError("T and dt must be equal")
		elif len(self.cov) != len(other.cov):
			raise ValueError("Tolerance may be ill-defined")
		else:
			new_cov = np.array(self.cov) | np.array(other.cov)
			return Coverage(list(new_cov),self.target,self.T,self.dt)

	@property
	def cov(self):
		return self._cov
	
	@property
	def target(self):
		return self._target
	
	@property
	def T(self):
		return self._T
	
	@property
	def dt(self):
		return self._dt
	
		

	#TODO: implement
	@property
	def accumulated(self):
		"""Accumulated time of view in seconds"""
		return self._accum()

	def _accum(self):
		return self.dt*sum(self.cov)


	@property
	def mean_gap_light(self):
		return self._mean_gap(1)

	@property 
	def mean_gap_dark(self)Ñ
		return self._mean_gap(0)

	def _mean_gap(self,view):
		gap = 0 
		switch = 0 
		c_gap = 0
		
		for c in self.cov:
			if switch == 0 and c == view:
				switch = 1
				c_gap += 1
				gap += 1
			elif switch == 1 and c != view:
				switch = 0
			elif switch == 1 and c == view:
				gap += 1
			elif switch == 0 and c != view:
				continue

		if c_gap == 0:
			return 0
		else:
			return gap*self.dt/c_gap	


	
	@property
	def response_time(self):	
		return self._resp_t()
	
	def _resp_t(self):
		idx = []
		resp_time = 0
		gap = 0

		for c in self.cov:
			if switch == 0 and c == 0:
				switch = 1
				gap += 1
			elif switch == 1 and c != 0:
				switch = 0
				resp_time += (gap+1)*gap/2
				gap = 0
			elif switch == 1 and c == 0:
				gap += 1
			elif switch == 0 and c != 0:
				continue

		return resp_time/len(self.cov)



		

	def plot_lapida(self, **kwargs):
		"""Grafico de lapidas. 

		Funcion escalon en los intervalos donde el objetivo es cubierto.
		En esencia es la misma funcion step de pyplot. En caso que quiera ser 

		Recibe: 
			**kwargs: argumentos propios de la funcion matplotlib.pyplot.step
			(para mas informacion recurrir a documentacion oficial de matplotlib
			https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.pyplot.step.html)
		Devuelve:
			 
		"""
		
		merged = self.merge()
		x = np.array( [self.dt*i/3600 for i in range(len(merged))] ) 
		y = np.array(merged)

		fig = plt.figure()
		plt.step(x,y,**kwargs)
		plt.title(str(self.acumulado()))

		return fig
