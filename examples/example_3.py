import os
import sys
sys.path.insert(0, os.path.abspath('..'))


import pi
from pi.WalkerDelta import constell
import pi.units as u


import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

ctll = constell.WalkerDelta(S=3,FOV=20,P=3)

fig=plt.figure()
ax=plt.axes(title=ctll.__str__())


step=0.2
lats = np.arange(0,90,step)

d_cov=[]
d_gap_0=[]
d_gap_1=[]
for lat in lats:
	cov = ctll.coverage(target=[lat,-30],dt=1)
	d_cov.append(cov.acumulado())
	d_gap_0.append(cov.gap_medio(vista=0))
	d_gap_1.append(cov.gap_medio(vista=1))
	print(lat)

d_cov = np.array(d_cov)
d_cov =  d_cov/np.sum(d_cov)
d_gap_0 = np.array(d_gap_0)
d_gap_0 = d_gap_0/np.sum(d_gap_0)
d_gap_1 = np.array(d_gap_1)
d_gap_1 = d_gap_1/np.sum(d_gap_1)
ax.plot(lats,d_cov,color='k',label="acumulado")
ax.plot(lats,d_gap_0,color='r',label="gap medio nocob")
ax.plot(lats,d_gap_1,color='g',label="gap medio cob")

ax.legend()
plt.show()