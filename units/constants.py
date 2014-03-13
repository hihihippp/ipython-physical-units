# Important physical constants
# TODO: 
# - check if physics_extension is loaded
# - make unload work
# - extend constants with explanation and reference

from math import pi

import numpy as np

# Standard constants 

c0 = Q(299792458., 'm/s')
mu0 = Q(4.e-7 ,'pi*N/A**2').base
eps0 = Q(8.854188e-12 ,'F/m').base
Grav = Q(6.67384e-11 ,'m**3/kg/s**2')
hpl = Q(6.62606957e-34 ,'J*s')
hbar = Q(6.62606957e-34 ,'J*s')/(2*pi)
e0 = Q(1.602176565e-19 ,'C')
me = Q(9.10938291e-31 ,'kg')
mp = Q(1.672621777e-27 ,'kg')
mn = Q(1.674927351e-27 ,'kg')
NA = Q(6.02214129e23 ,'1/mol')
kb = Q(1.3806488e-23 ,'J/K')
g0 = Q(9.80665 ,'m/s**2')


