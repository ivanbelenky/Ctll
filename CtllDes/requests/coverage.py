from ..core import ctll, Instrument, sat
import collections.abc
from collections.abc import Iterable


try:
    from functools import cached_property  # type: ignore
except ImportError:
    from cached_property import cached_property  # type: ignore




def coverages_from_sat(sat, targets,T, dt=1.):
	"""Build list of coverage objects from satellite
	
	Parameters
	----------
	sat : ~CtllDes.core.sat
		sat object
	targets : ~CtllDes.targets.Targets 
		Desired targets of coverage
	T : float
		Desired Time of analysis in days.
	dt : float, optional
		time interval between interpolations of ssp
	
	Returns
	-------
	Coverages : list
		List of Coverage objects, one for each target, and one for 
		each coverage instrument in the satellite.

	"""


	CovInstruments = [sat.instr for instr in sat.CovInstruments]
	if not len(CovInstruments):
		raise Exception("No coverage instruments found on" +
			f" satID : {sat.id}")	

	ssps = sat.ssps(T,dt)
	r,v = sat.rv(T,dt)

	Coverages = []
	for instr in CovInstrument:
		for target in targets:
			cov = isCovered(ssps,r,target,sat.attractor.R_mean,instr.coverage())
			Coverages.append(Coverage(cov,target,instr.id,dt))

	return Coverages

def coverages_from_ctll(ctll,targets,T,dt=1.):
	"""Get coverages from Constellation Object.
		
	"""
	Coverages = []
	for sats in ctll.sats:
		try:
			covs = coverages_from_sat(sat,targets,T,dt)
		except Exception:
			pass
		else:
			Coverages.append(covs)

	if not len(Coverages): 
		raise Exception("Constellation has no Coverage Instruments")

	return Coverages

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



COVERAGE_COMPARING_METHODS = [symmetric, symmetric_with_roll]


class Coverages(collections.abc.Set):
		def __init__(self,covs,tag=None):
		self._targets = lst = list()
		self._tag = tag if tag else "No Tag"

		if isinstance(covs,Coverage):
			lst.append(covs)
		elif not isinstance(covs,Iterable):
			raise TypeError("covs must be Coverage, or iterable collection of Coverage objects")
		else: 	
			for cov in covs:
				if not isinstance(cov,Coverage):
					raise TypeError("Targets must be Target or iterable collection of Target objects") 
				if cov not in lst:
					lst.append()

	@property
	def covs(self):
		return self._covs
	
	def __iter__(self):
		return iter(self.covs)

	def __contains__(self,value):
		return value in self.covs

	def __len__(self):
		len(self.covs)

	@cached_property
	def data(self):
		return self.to_data()
		#TODO: check if it is worth doing named tuple or other thing

	#TODO: implement
	@classmethod
	def from_ctll
		"""TODO docstring"""
		pass

	#TODO: implement
	@classmethod
	def from_sat
		"""TODO docstring"""
		pass

	#TODO: implement
	def filter_by_target(self,target):
		#TODO: implement	
		pass

	#TODO: implement
	def filter_by_instrument(self,instrument):
		#TODO: implement	
		pass

	#TODO: implement
	def collapse(self):
		#TODO: implement	
		pass


class Coverage(object):
	#TODO: docstring
	
	def __init__(self,
		cov,
		target,
		dt
	):
		#TODO: T should be added
		self._cov = cov
		self._target = target 
		self._dt = dt
		
	#TODO: implement
	@cached_property
	def accumulated(self):
		return self._accum()

	#TODO: implement
	@cached_property
	def mean_gap(self):
		return self._mean_gap()

	#TODO: implement
	@cached_property
	def response_time(self):
		return self._resp_t()
	
	#TODO: implement
 
	#TODO: do this functions all over again, change language	
	def acumulado(self):
		"""Tiempo acumulado de cobertura.

		El tiempo acumulado de cobertura depende de la resolucion temporal
		
		Recibe:

		Devuelve:
			acum: tiempo acumulado de cobertura en el lapso dado."""

		acum = 0
				
		for c in merged:
			acum += c*dt
		return acum


	def gap_medio(self, vista = 1):
		"""
		Calcula gap medio

		Devuelve el valor del gap medio de cobertura o no cobertura respectivamente
		dependiendo del valor de vista. 

		Recibe: 
			vista: si es 1, calculara el valor medio del gap de cobertura. 
			el complemento en caso de ser 0
		Devuelve:
			promedio de duracion del gap, en segundos."""

		merged = self.merge()
		
		gap = 0 
		switch = 0 
		c_gap = 0
		
		for c in merged:
			if switch == 0 and c == vista:
				switch = 1
				c_gap += 1
				gap += 1
			elif switch == 1 and c != vista:
				switch = 0

			elif switch == 1 and c == vista:
				gap += 1
			elif switch == 0 and c != vista:
				continue

		if c_gap == 0:
			return 0
		else:
			return gap*self.dt/c_gap	
		

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
