"""
Stellar Astrophysics HW6.

Jonas Powell
May 14th
"""

# Import some packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.constants as c
import emcee

# import the differential functions for the integrator
from utils_hw6 import chi2
from T_Tauri import structure_integrator


# Define some general physical constants
m_electron = 9.1e-31                # kg
h = c.h.value                       # J s
G = c.G.value                       # m3 kg-1 s-2
c_vel = 3e8                         # m/s
sig = c.sigma_sb.value              # W / (K4 m2)

mSol = c.M_sun.value
lSol = 3.848e28
kg2earthMass = 1.6744e-25           # kg in units of earth mass
kg2solarMass = 5.0291e-31
m2solarRadius = 1.437e-9

# Define some physical characteristics for our star
luminosity = 1 * lSol           # Watts
mmw = 0.6 * c.u.value           # Mean molecular weight times H mass
gamma = 5./3.                   # For Low Mass white dwarf
m_star = 0.5 * mSol             # kg
t_star = 4000                   # Kelvin
r_star = np.sqrt(luminosity/(4*np.pi * sig * t_star**4))

# Some initial values
m_init = 0.
r_init = 0.00001


"""
Using Ismael's values:
K_ish = 4.1014e8
rho0_ish = 464.097
and P = K * rho0**(5./3.)

we find P0 = 1.141e13 (mks)
"""


# Run an MCMC to find the best values of rho0, K
def lnprob(fitters):
    """
    Get [m, r] predictions for the star of given params.
    """
    
    [K, rho0] = fitters

    if K < 1 or K > 10:
        return -np.inf

    if rho0 < 1 or rho0 > 1e20:
        return -np.inf

    [m, r] = structure_integrator(gamma, K, r_init, m_init, rho0)

    chisq = chi2([m, r], [m_star, r_star])
    ln_p = -0.5 * chisq

    print "New chi2 val: ", chisq
    return ln_p


ndim, nwalkers = 2, 10
# Give initial positions in [K, rho0] space for the walkers
p0 = [np.random.rand(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
sampler.run_mcmc(p0, 100)







# The End
