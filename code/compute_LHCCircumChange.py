#python script to compute LHC Circumference changes
#copyright R. Rietbroek
#Date: 20 June 2016
#License: pending opensource  license to be determined, all rights reserved

from SphereHarm import Pnm,sanity_check
import numpy as np
np.set_printoptions(linewidth=130)


sanity_check()

P10=Pnm(8)

Pvals=P10.at(2/3)
print(Pvals)

Pdotvals=P10.derivTheta(2/3.5)
print(Pdotvals)
