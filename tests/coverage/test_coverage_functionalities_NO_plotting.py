import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

import targets as t
import coverage as c
import CtllDes 
from CtllDes.core import ctll, satellite

p = 7000 * u.km
ecc = 0 * u.one 
inc = 95 * u.deg
raan = 0 * u.rad
argp = 0 * u.rad
nu = 0 * u.rad
plane = Planes.EARTH_EQUATOR 

pi = np.pi

def test_coverage_basic():


	#constellation
	T = 10
	P = 5
	F = 0
	constellation = ctll.Ctll.from_WalkerDelta(T,P,F,p,ecc,inc,argp)
	
	tgts = t.Targets.from_country('Argentina',N=50)
	covs = c.Coverage.from_ctll(constellation,tgts,10)
	covs.plot()

	#get data from cov, that is a list of named tuples, with targets, coverages figures
	#y el instrument id que use
	
	covs_data = covs.data
	acumulated = covs_data.accumulated
	mean_gap = covs_data.mean_gap

