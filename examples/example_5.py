import os
import sys
sys.path.insert(0, os.path.abspath('..'))


import pi
from pi.WalkerDelta import constell
from pi import status as st

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


fig=plt.figure()



dda=0.0005
x = np.arange(0,100,1)
y = []
for ex in x:
	ctll = constell.WalkerDelta(S=1,FOV=20,P=1)
	y.append(ctll.coverage(target=[70,-30],dt=1).gap_medio())
	ctll.degrade(j=0,da=dda*ex)
	print("done")

ax=plt.axes(title=ctll.__str__())
ax.legend()
ax.plot(x*0.0005*100,y)
ax.set_xlabel("porcentaje degradadacion [°]")
ax.set_ylabel("tiempo medio de cobertura [s]")
plt.show()