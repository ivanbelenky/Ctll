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




def test_empty_constellation_constructor():
	ctll = CTLL.Ctll(state)

	print("this are the specs:",ctll.specs)

	# ctll.Info(v=True)