import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

from astropy import units as u

from poliastro.twobody import states as st
from poliastro.bodies import Earth, Mars, Sun
from poliastro.frames import Planes


import numpy as np

import CtllDes 
from CtllDes.core import CTLL, SAT


p = 35000 * u.km
ecc = 0 * u.one
inc = 0 * u.rad
raan = 0 * u.rad
argp = 0 * u.rad
nu = 0 * u.rad
plane = Planes.EARTH_EQUATOR 

pi = np.pi




def test_one_satellite_ctll_without_specs_instruments_status():
	
	state = [st.ClassicalState(Earth,p,ecc,inc,raan,argp,nu,plane)]
	ctll = CTLL.Ctll(state)

	print("this are the constellation intrumentss",ctll.instrumentss)
	print("this are the constellation specs",ctll.specs)
	print("this are the constellation patterns",ctll.pattern)
	print("this is the constellation status",ctll.statuss)
	print("Online satellite are",ctll.getOnlineSatsId())
	
	ctll.Info(v=True)

def test_one_satellite_ctll_with_status_without_specs_instruments():
	state = [st.ClassicalState(Earth,p,ecc,inc,raan,argp,nu,plane)]
	ctll = CTLL.Ctll(state)
	
	print("this are the constellation intrumentss",ctll.instrumentss)
	print("this are the constellation specs",ctll.specs)
	print("this are the constellation patterns",ctll.pattern)
	print("this is the constellation status",ctll.statuss)
	print("Online satellite are",ctll.getOnlineSatsId())
	
	ctll.Info(v=True)

def test_arbitrary_satellite_without_specs_isntruments_status():

	states = []
	N = 6
	for i in range(1,N):
		states.append(st.ClassicalState(Earth,p,ecc,inc,raan,argp,nu+2*pi/i*u.rad,plane))
	
	ctll = CTLL.Ctll(states)

	print("this are the constellation intrumentss",ctll.instrumentss)
	print("this are the constellation specs",ctll.specs)
	print("this are the constellation patterns",ctll.pattern)
	print("this is the constellation status",ctll.statuss)
	print("Online satellite are",ctll.getOnlineSatsId())
	
	ctll.Info(v=True)

def test_ctll_class_method_from_Walker_Delta_without_instruments_specs_status():
	T = 10
	P = 5
	F = 0

	ctll = CTLL.Ctll.from_WalkerDelta(T,P,F,p,ecc,inc,argp)

	print("this are the constellation intrumentss",ctll.instrumentss)
	print("this are the constellation specs",ctll.specs)
	print("this are the constellation patterns",ctll.pattern)
	print("this is the constellation status",ctll.statuss)
	print("Online satellite are",ctll.getOnlineSatsId())
	
	ctll.Info(v=True)


