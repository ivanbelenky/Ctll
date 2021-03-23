import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

import CtllDes 
from CtllDes.core import CTLL  



def test_empty_constellation_constructor():
	ctll = CTLL.Ctll()
