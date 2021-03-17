import os
import sys
sys.path.insert(0, os.path.abspath('..'))


import pi
from pi.WalkerDelta import constell
from pi import status as st

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

ctll = constell.WalkerDelta(S=1,FOV=20,P=1)

fig=plt.figure()
ax=plt.axes(title=ctll.__str__())


step=0.1
lats = np.arange(0,90,step)

jet= plt.get_cmap('hot')
colors = iter(jet(np.linspace(0,1,12)))

dda=0.008
for i in range(5):
	aux=ctll
	d_cov = np.array([aux.coverage(target=[lat,0],dt=5).acumulado()
	 for lat in lats])
	#d_cov = d_cov/sum(d_cov)/step
	ax.plot(lats,d_cov,color=next(colors),
		label=f'Semiejemayor degradado: {100-(1-dda)**i*100:.2f} %')
	aux.degrade(j=0,da=dda)
	print(aux.sats[0].a)

ax.legend()
ax.set_xlabel("latitud [°]")
ax.set_ylabel("tiempo acumulado de cobertura [s]")
plt.show()