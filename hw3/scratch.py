"""
HW 3 Scratch
Stellar
2.24.18
"""


import scipy.constants as c
import numpy as np
import astropy.units as u

rho_0 = (4.55*10**(-31) * u.g * (u.cm)**-3).to(u.kg * u.m**-3)
T_0 = 2.7 * u.Kelvin

m_proton = c.proton_mass * u.kg
m_electron = c.electron_mass * u.kg
h = c.Planck * u.Joule * u.second
k_B = c.Boltzmann * u.Joule/u.Kelvin
Chi = (13.6*u.electronvolt).to(u.Joule)



# Note that proton mass --> hydrogen mass (electron negligeable)
consts1 = (8 * m_proton) / (3 * rho_0) * ( 2 * np.pi * m_electron * k_B * T_0 * h**-2) **1.5
consts1


consts2 = Chi / (k_B * T_0)
consts2


(13.6*u.electronvolt).to(u.Joule)
c.Boltzmann * u.Joule * u.Kelvin**-1


# Solve x^-3 = x^-1.5 * 1.05*10^23 * e^(-58452 * x) for x
a = 1./1386.
a = 0.0007215007215

rho_0
n_e = (3. * rho_0 * a**-3)/(8. * m_proton)
n_e
sigma_electron = 6.65 * 10**-29 * u.m**2
s_mfp = 1/(n_e * sigma_electron)

s_mfp







# The End
