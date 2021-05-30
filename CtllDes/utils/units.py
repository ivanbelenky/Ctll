import numpy as np

from astropy import units as u
from poliastro.bodies import Earth

#Constantes del esquema iterativo
ITERMAX = 1000
TOL = 10**-8

#Constantes matemáticas
pi = np.pi

#Constantes físicas para la tierra
G = 6.6743*(10**(-11)) 
M = Earth.mass.to(u.kg).value
R_e = Earth.R.to(u.m).value
mu = G*M

J_2 = +0.00108263
J_4 = 0 #-1.649*10**-6

#Constantes temporales 
s_m = 86460 # tiempo medio en un día solar, segundos
D_sid = 365.256360 # días siderales por año solar
w_s = 2*pi/(s_m*D_sid) # frecuencia angular deseada para la línea de nodos
#en orbita heliosincr[onica]

we = 7.2921150*10**-5

rad2deg = 180/pi