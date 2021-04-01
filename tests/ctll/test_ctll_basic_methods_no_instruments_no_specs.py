import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

from astropy import units as u

from poliastro.twobody import states
from poliastro.bodies import Earth, Mars, Sun
from poliastro.frames import Planes


import CtllDes 
from CtllDes.core import CTLL  



p = 35000 * u.km
ecc = 0 * u.one
inc = 0 * u.rad
raan = 0 * u.rad
argp = 0 * u.rad
nu = 0 * u.rad
plane = Planes.EARTH_EQUATOR 

state = [states.ClassicalState(Earth,p,ecc,inc,raan,argp,nu,plane)]




def test_one_satellite_ctll_without_specs():
	ctll = CTLL.Ctll(state)
	print("this are the constellation intrumentss",ctll.instrumentss)
	print("this are the constellation specs",ctll.specs)
	print("this are the constellation patterns",ctll.pattern)
	print("this is the constellation status",ctll.statuss)
	ctll.Info(v=True)

def test_