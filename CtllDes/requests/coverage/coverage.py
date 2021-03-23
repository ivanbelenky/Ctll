from ..core import CTLL, Instrument, SAT


def coverages_from_SAT(SAT, targets,T, dt=1.):
	"""Build list of coverage objects from satellite
	
	Parameters
	----------
	SAT : ~CtllDes.core.SAT
		SAT object
	targets : ~CtllDes.targets 
		Desired targets of coverage
	T : float
		Desired Time of analysis in days.
	dt : float, optional
		time interval between interpolations of ssp
	
	Returns
	-------
	Coverages : list
		List of Coverage objects. 

	"""


	CovInstruments = [SAT.instr for instr in SAT.CovInstruments]
	if not len(CovInstruments):
		raise Exception("No coverage instruments found on" +
			f" SATID : {SAT.id}")	

	ssps = SAT.ssps(T,dt)

	Coverages = []
	for instr in CovInstrument:
		for target in targets:
			cov = isCovered(ssps,target,instr._coverage)
			Coverages.append(Coverage(cov,target,instr.id))

	return Coverages

def coverages_from_CTLL(CTLL,targets,T,dt=1.):
	"""Get coverages from Constellation Object.

	"""
	Coverages = []
	for sats in CTLL.sats:
		try:
			covs = coverages_from_SAT(sat,targets,T,dt)
		except Exception:
			pass
		else:
			Coverages.append(covs)

	if not len(Coverages): 
		raise Exception("Constellation has no Coverage Instruments")

	return Coverages

def isCovered(ssps,target,CoverageMethod):



class Coverage(object):
	
	def __init__(self,
		cov,
		target,
		instrumentID = None
	):

		self._cov = cov
		self._target = target 
		self._instrumentID = instrumentID if instrumentID 
		else uuid.uuid4()


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
