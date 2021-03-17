import os
import sys
sys.path.insert(0, os.path.abspath('..'))


import pi
from pi.WalkerDelta import constell


import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

ctll = constell.WalkerDelta(S=1,FOV=30,P=1)

latshape = 120
longshape = 200

lats = np.array([180*i/latshape-90 for i in range(latshape+1)])
longs = np.array([360*i/longshape for i in range(longshape+1)])
a = np.zeros(shape = (latshape+1,longshape+1))


for i in range(latshape+1):
	for j in range(longshape+1):
		cobaux = ctll.coverage(target=[lats[i],longs[j]],dt=10)
		print("lat: ",i,"	long:",j)
		a[i][j] = np.log(cobaux.acumulado())
		#a[i][j] = np.log(cobaux.gap_medio(vista=0))

X,Y = np.meshgrid(longs,lats)
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())

ax.coastlines()
ctll.plot_gt()

plt.pcolormesh(X,Y,a,cmap='magma')
plt.title(ctll.__str__())
plt.colorbar()
#plt.imshow(a, cmap='Reds', interpolation = 'bilinear')
plt.show()