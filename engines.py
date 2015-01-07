# Python 2.7.6
import math

NA= 6.022e23 # Avogadro's number
k= 1.3806488e-23 # in SI base units
atm= 101325 # atmospheric pressure in Pa
L= .001 # 1 L = .001 m^3

class Deltas:
  # A delta object contains {dW, dQ, dS} for work performed by gas,
  # heat absorbed by the gas, and entropy increase of the universe, respectively
  # d represents delta; these are not necessarily infinitesimal changes.
  def __init__(self,dW,dQ,dS):
    self.dW= dW
    self.dQ= dQ
    self.dS= dS
  def plus(self,d):
    return Deltas(self.dW+d.dW, self.dQ+d.dQ, self.dS+d.dS)
class Gas:
  def __init__(self,N,f):
    self.N= float(N)
    self.f= float(f)
    self.gamma= (self.f+2)/self.f
  def getU(self):
    return self.N*k*self.f*self.T/2
  def set_pV(self,p,V):
    self.p= float(p)
    self.V= float(V)
    self.T= self.p*self.V/(self.N*k)
  def set_VT(self,V,T):
    self.V= float(V)
    self.T= float(T)
    self.p= self.N*k*self.T/self.V
  def set_pT(self,p,T):
    self.p= float(p)
    self.T= float(T)
    self.V= self.N*k*self.T/self.p
  def adiabat_T(g,T):
    # Given gas g and target temperature T, returns delta dictionary for an
    # adiabat to the target temperature. Also moves g to the appropriate
    # endpoint in pVT-space.
    Ui= g.getU()
    g.T= T
    dW= Ui - g.getU()
    # Still need to set new p, V
    g.V= (g.p*g.V**g.gamma/(g.N*k*g.T))**(1/(g.gamma-1))
    g.set_VT(g.V, g.T)
    return Deltas(dW,0,0)
  def adiabat_p(g,p):
    # Gas g and target pressure p
    K= g.p*g.V**g.gamma
    Vi= g.V
    g.V= Vi*(g.p/p)**(1/g.gamma)
    dW= (K/(1-g.gamma))*((g.V)**(1-g.gamma)-Vi**(1-g.gamma))
    g.set_pV(p, g.V)
    return Deltas(dW,0,0)
  def adiabat_V(g,V):
    # Gas g and target volume V
    K= g.p*g.V**g.gamma
    dW= (K/(1-g.gamma))*(V**(1-g.gamma)-g.V**(1-g.gamma))
    g.p= K/(V**g.gamma)
    g.set_pV(g.p, V)
    return Deltas(dW,0,0)
  def isotherm_p(g,p):
    # Gas g and target pressure p
    Vi= g.V
    g.set_pT(p, g.T)
    dW= g.N*k*g.T*math.log(g.V/Vi)
    dQ= dW
    dS= abs(dQ/g.T)
    return Deltas(dW,dQ,dS)
  def isotherm_V(g,V):
    # Gas g and target volume V
    dW= g.N*k*g.T*math.log(V/g.V)
    g.set_VT(V, g.T)
    dQ= dW
    dS= abs(dQ/g.T)
    return Deltas(dW,dQ,dS)
  def isochor_p(g,p):
    dW= 0
    Ti= g.T
    Ui= g.getU()
    g.set_pV(p, g.V)
    dQ= g.getU() - Ui
    dS= abs((g.N*k*g.f/2)*log(g.T/Ti))
    return Deltas(dW,dQ,dS)
  def isochor_T(g,T):
    dW= 0
    Ti= g.T
    Ui= g.getU()
    g.set_VT(g.V, T)
    dQ= g.getU() - Ui
    dS= abs((g.N*k*g.f/2)*log(g.T/Ti))
    return Deltas(dW,dQ,dS)
  def isobar_V(g,V):
    Ti= g.T
    Vi= g.V
    g.set_pV(g.p, V)
    dW= g.p*(g.V-Vi)
    dQ= g.N*k*((g.f+2)/2)*(g.T-Ti)
    dS= abs(((g.f+2)/2)*g.N*k*math.log(g.T/Ti))
    return Deltas(dW,dQ,dS)
  def isobar_T(g,T):
    Ti= g.T
    Vi= g.V
    g.set_pT(g.p, T)
    dW= g.p*(g.V-Vi)
    dQ= g.N*k*((g.f+2)/2)*(g.T-Ti)
    dS= abs(((g.f+2)/2)*g.N*k*math.log(g.T/Ti))
    return Deltas(dW,dQ,dS)
    
