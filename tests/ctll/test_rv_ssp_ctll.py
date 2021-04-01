import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

from astropy import units as u

from poliastro.twobody import states as st
from poliastro.bodies import Earth, Mars, Sun
from poliastro.frames import Planes


import numpy as np
import matplotlib.pyplot as plt
import tkinter
import matplotlib
matplotlib.use('TkAgg')

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

def test_rv_ctll():


	#constellation
	T = 10
	P = 10
	F = 0
	p2 = p * u.km
	constellation = ctll.Ctll.from_WalkerDelta(T,P,F,p,ecc,inc,argp)
	constellation.info(v=True)

	rv = constellation.rv(6)

	for rrvv in rv:
		print("\nr:\n",rrvv[0])
		print("\nv:\n",rrvv[1])

	print(type(rv[0]))

	fig= plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.scatter(0,0,0,zdir='z',s=10,c='r')
	for rrvv in rv:
		ax.plot(rrvv[0][::100,0],rrvv[0][::100,1],rrvv[0][::100,2],zdir='z',c='k')
	plt.show()

def test_ssps_ctll():


	#constellation
	T = 1
	P = 1
	F = 0
	constellation = ctll.Ctll.from_WalkerDelta(T,P,F,p,ecc,inc,argp)
	constellation.info(v=True)

	latlon = constellation.ssps(0.041)

	fig = plt.figure()

	print("type of latlon objects", type(latlon[0][1]))
 
	plt.scatter(latlon[0][1],latlon[0][0],s=0.3,c='k')
	plt.xlabel("longitude [rad]")
	plt.ylabel("latitude [rad]")
	plt.show()




