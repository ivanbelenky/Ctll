import os
import sys
sys.path.insert(0, os.path.abspath('..'))


import pi
from pi.WalkerDelta import constell


import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

ctll = constell.WalkerDelta(S=1,FOV=20,P=5)

fig=plt.figure()
ax=plt.axes(title=ctll.__str__())


step=2
lats = np.arange(0,90,step)
longs = np.arange(0,360,40)
jet= plt.get_cmap('jet')
colors = iter(jet(np.linspace(0,1,12)))

for longg in longs:
	d_cov = np.array([ctll.coverage(target=[lat,longg],dt=5).acumulado()
	 for lat in lats])
	d_cov = d_cov/sum(d_cov)/step
	ax.plot(lats,d_cov,color=next(colors),label=longg)

ax.legend()
plt.show()