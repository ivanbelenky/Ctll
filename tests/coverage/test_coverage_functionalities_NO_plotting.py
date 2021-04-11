import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

from astropy import units as u


from poliastro.twobody import states as st
from poliastro.bodies import Earth, Mars, Sun
from poliastro.frames import Planes
from poliastro.twobody import propagation,Orbit

import numpy as np

import CtllDes
import CtllDes.targets.targets as t
 
from CtllDes.core import ctll, satellite, instrument
from CtllDes.requests.coverage import Coverages



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
	T = 4
	P = 2
	F = 0

	cam = instrument.Camera(0,0)
	
	constellation = ctll.Ctll.from_WalkerDelta(T,P,F,p,ecc,inc,argp,instrumentss=cam)
	
	tgts = t.Targets.from_country('Argentina',N=4)
	
	covs = Coverages.from_ctll(constellation,tgts,5)
	
	df = covs.to_df()

	print(df)

	#covs.plot()

	#get data from cov, that is a list of named tuples, with targets, coverages figures
	#y el instrument id que use
	
	#covs_data = covs.data
	#acumulated = covs_data.accumulated
	#mean_gap = covs_data.mean_gap

