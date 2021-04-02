import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

import matplotlib.pyplot as plt
from shapely.geometry import Point
from targets.targets import Targets, Target 

def test_targets_builds_from_point():
	point = Point(1,2)
	tgts = Targets(Point,tag='test')
	tgts.plot()
	tgts.plot(use_3d=True)
	plt.show()

def test_targets_builds_from_point():
	tgt = Target(1,2)
	tgts = Targets(tgt,tag='test')
	tgts.plot()
	tgts.plot(use_3d=True)
	plt.show()

def test_targets_build_from_country():
	tgts = Targets.from_country('Argentina')
	tgts.plot()
	tgts.plot(use_3d=True)
	plt.show()

def test_targets_build_from_state():
	tgts = Targets.from_country('Buenos Aires')
	tgts.plot()
	tgts.plot(use_3d=True)
	plt.show()
