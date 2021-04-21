import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

import matplotlib.pyplot as plt
from shapely.geometry import Point
import CtllDes
from CtllDes.targets.targets import Targets, Target 



def test_targets_builds_from_point():
	point = Point(1,2)
	tgts = Targets(point,tag='test')
	print(point.x)
	tgts.plot()
	try:
		tgts.plot(use_3d=True)
	except:
		pass

	plt.show()

def test_targets_builds_from_Target():
	tgt = Target(1,2)
	tgts = Targets(tgt,tag='test')
	tgts.plot()
	try:
		tgts.plot(use_3d=True)
	except:
		pass
	plt.show()

def test_targets_build_from_country():
	tgts = Targets.from_country('Argentina')
	fig=tgts.plot()
	try:
		tgts.plot(use_3d=True)
	except:
		pass
	plt.show()

def test_targets_build_from_state():
	tgts = Targets.from_state('Buenos Aires')
	tgts.plot()
	try:
		tgts.plot(use_3d=True)
	except:
		pass
	plt.show()
