"""
Find the structure of a T_Tauri star.

Takes values for K, gamma, initial values
Returns mass, radius
"""

# Import some packages
import numpy as np
import pandas as pd
import astropy.constants as c

# Some initial values
m_init = 0.
r_init = 0.00001


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


# Polytropic EOS
def get_new_rho(P, K, gamma):
    rho = (P/K)**(1./gamma)
    return rho


# Mass Conservation: dM/dr = 4pi r**2 rho
def dm_dr(P, K, rho, r, gamma):
    rho = get_new_rho(P, K, gamma)
    dM = 4. * np.pi * rho * r**2.
    return dM


# Hydrostatic Equilibrium: dp/dr = -GM*rho*r**-2
def dp_dr(m, rho, r):
    dP = -G * m * rho * r**(-2.)
    return dP


# Get Temperature at a given radius from pressure and density.
def get_temp(P, K, rho, r):
    T = (P * mmw)/(rho * K)
    return T


# An integrator to integrate those differential functions
def structure_integrator(gamma, K, r0, m0, rho0):
    # Find the starting central pressure
    p0 = K * (rho0 ** gamma)
    t0 = get_temp(p0, K, rho0, r0)

    # Set up a record book for the steps we take
    conditions = []
    initialConditions = {'Radius': r0,
                         'Mass': m0,
                         'Pressure': p0,
                         'Density': rho0,
                         'Temperature': t0
                         }
    conditions.append(initialConditions)

    # Give a step size
    dr = 1.

    # Feed it the initial conditions
    r, m, p, rho = r0, m0, p0, rho0

    while p > 0:
        # Determine the new conditions:
        r_new = r + dr
        m_new = m + (dr * dm_dr(p, K, rho, r, gamma))
        p_new = p + (dr * dp_dr(m, rho, r))
        rho_new = rho + (dr * get_new_rho(p, K, gamma))
        t_new = get_temp(p_new, K, rho_new, r_new)

        # Records them
        new_step = {'Radius': r_new,
                    'Mass': m_new,
                    'Pressure': p_new,
                    'Density': rho_new,
                    'Temperature': t_new
                    }
        conditions.append(new_step)

        # Prepare for iteration
        r = r_new
        m = m_new
        p = p_new
        rho = rho_new
        t = t_new

    conditions_df = pd.DataFrame(conditions)

    # We can just return the last steps from the loop
    return [r, m]
