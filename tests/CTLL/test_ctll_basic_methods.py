import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

from astropy import units as u

from poliastro.twobody import states as st
from poliastro.bodies import Earth, Mars, Sun
from poliastro.frames import Planes


import numpy as np

import CtllDes 
from CtllDes.core import ctll, satellite


p = 35000 * u.km
ecc = 0 * u.one
inc = 0 * u.rad
raan = 0 * u.rad
argp = 0 * u.rad
nu = 0 * u.rad
plane = Planes.EARTH_EQUATOR 

pi = np.pi




def test_one_satellite_constellation_without_specs_instruments_status():
	
	state = [st.ClassicalState(Earth,p,ecc,inc,raan,argp,nu,plane)]
	constellation = ctll.Ctll(state)

	print("this are the constellation intrumentss",constellation.instrumentss)
	print("this are the constellation specs",constellation.specs)
	print("this are the constellation patterns",constellation.pattern)
	print("this is the constellation status",constellation.statuss)
	print("Online satellite are",constellation.online_id)
	
	constellation.info(v=True)

def test_one_satellite_constellation_with_status_without_specs_instruments():
	state = [st.ClassicalState(Earth,p,ecc,inc,raan,argp,nu,plane)]
	constellation = ctll.Ctll(state)
	
	print("this are the constellation intrumentss",constellation.instrumentss)
	print("this are the constellation specs",constellation.specs)
	print("this are the constellation patterns",constellation.pattern)
	print("this is the constellation status",constellation.statuss)
	print("Online satellite are",constellation.online_id)
	
	constellation.info(v=True)

def test_arbitrary_satellite_without_specs_isntruments_status():

	states = []
	N = 6
	for i in range(1,N):
		states.append(st.ClassicalState(Earth,p,ecc,inc,raan,argp,nu+2*pi/i*u.rad,plane))
	
	constellation = ctll.Ctll(states)

	print("this are the constellation intrumentss",constellation.instrumentss)
	print("this are the constellation specs",constellation.specs)
	print("this are the constellation patterns",constellation.pattern)
	print("this is the constellation status",constellation.statuss)
	print("Online satellite are",constellation.online_id)
	
	constellation.info(v=True)

def test_constellation_class_method_from_Walker_Delta_without_instruments_specs_status():
	T = 10
	P = 5
	F = 0

	constellation = ctll.Ctll.from_WalkerDelta(T,P,F,p,ecc,inc,argp)

	print("this are the constellation intrumentss",constellation.instrumentss)
	print("this are the constellation specs",constellation.specs)
	print("this are the constellation patterns",constellation.pattern)
	print("this is the constellation status",constellation.statuss)
	print("Online satellite are",constellation.online_id)
	
	constellation.info(v=True)




def test_adding_constellation_and_individual_satellites():

	#constellation1
	states = []
	N = 6
	for i in range(1,N):
		states.append(st.ClassicalState(Earth,p,ecc,inc,raan,argp,nu+2*pi/i*u.rad,plane))	
	constellation1 = ctll.Ctll(states)
	
	#constellation2
	T = 8
	P = 2
	F = 3
	p2 = p - 5000 * u.km
	constellation2 = ctll.Ctll.from_WalkerDelta(T,P,F,p,ecc,inc,argp)
	 
	constellation = constellation1 + constellation2 

	print("this are the constellation intrumentss",constellation.instrumentss)
	#print("this are the constellation specs",constellation.specs)
	print("this are the constellation patterns",constellation.pattern)
	print("this is the constellation status",constellation.statuss)
	print("Online satellite are",constellation.online_id)
	
	constellation.info(v=True)





