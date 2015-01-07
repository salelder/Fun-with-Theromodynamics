# Python 2.7.6
# Simulates a Carnot engine and tests that efficiency is 1 - TC/TH
from engines import *

TC= 273.0
TH= 400.0

# One mole of a monatomic gas
g= Gas(NA, 3)
g.set_pT(atm, TC)
d1= g.adiabat_T(TH)
d2= g.isotherm_V(25*L)
d3= g.adiabat_T(TC)
d4= g.isotherm_p(atm)

dtot= d1.plus(d2).plus(d3).plus(d4)
Qhot= d2.dQ
W= dtot.dW

print("Simulated work = "+str(W)+" J")
print("Simulated efficiency = "+str(W/Qhot))
print("Theoretical efficiency = "+str(1.0-TC/TH))
