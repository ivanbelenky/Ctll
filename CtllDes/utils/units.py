import numpy as np



#Constantes del esquema iterativo
ITERMAX = 100
TOL = 10**-5

#Constantes matemáticas
pi = np.pi

#Constantes físicas para la tierra
G = 6.67408*(10**(-11)) 
M = 5.9722*(10**24)
R_e = 6371*1000
mu = G*M

J_2 = +1.08262668*10**-3
J_4 = -1.649*10**-6

#Constantes temporales 
s_m = 86400 # tiempo medio en un día solar, segundos
D_sid = 365.256360 # días siderales por año solar
w_s = 2*pi/(s_m*D_sid) # frecuencia angular deseada para la línea de nodos
#en orbita heliosincr[onica]

we = 7.2921150*10**-5

rad2deg = 180/pi