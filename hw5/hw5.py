"""
Jonas Powell
Stellar HW5
April 22, 2018
"""

# PACKAGES
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as c
import seaborn as sns
import pandas as pd


# CONSTANTS
m_electron = 9.1e-31                # kg
m_h = 1.67e-27                      # kg
h = c.h.value
G = c.G.value
c_vel = 3e8                         # m/s
lm_K = h**2./(5.*m_electron) * (3./(8.*np.pi))**(2./3.) * (1./(2.*m_h))**(5./3.)
hm_K = (h*c_vel / 4.) * (3. / (8.*np.pi)) ** (1./3.) * (1./(2. * m_h))**(4./3.)
lm_gamma = 5./3.                    # For Low Mass white dwarf
hm_gamma = 4./3.

hm_K

kg2earthMass = 1.6744e-25           # kg in units of earth mass
kg2solarMass = 5.0291e-31
m2solarRadius = 1.437e-9

r_init = 0.00001
m_init = 0.
rho_init = 1e8                      # kg m-3


# DIFFERENTIAL FUNCTIONS

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



# An integrator to integrate those differential functions
def Integrator(gamma, K, r0, m0, rho0):
    # Find the starting central pressure
    p0 = K * (rho0 ** gamma)

    # Set up a record book for the steps we take
    conditions = []
    initialConditions = {'Radius': r0,
                         'Mass': m0,
                         'Pressure': p0,
                         'Density': rho0
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

        # Records them
        new_step = {'Radius': r_new,
                    'Mass': m_new,
                    'Pressure': p_new,
                    'Density': rho_new
                    }
        conditions.append(new_step)

        # Prepare for iteration
        r = r_new
        m = m_new
        p = p_new
        rho = rho_new

    conditions_df = pd.DataFrame(conditions)
    return conditions_df


# A quick way to calulate the mean density
def get_mean_density(df):
    mean_density = 3 * df['Mass']/(4 * np.pi * df['Radius']**3)
    return mean_density


# The beginning of some plotters.
def plot_m_vs_r(df):
    ms = df['Mass'] * kg2solarMass
    rs = df['Radius'] * m2solarRadius

    sns.set_style('white')
    plt.plot(ms, rs)
    sns.despine(trim=True)
    plt.show(block=False)


def mass_rad_relation():
    # Initialize the output dataframe
    out = []
    # Set up the list of densities to query
    rhos = 10**(np.arange(9, 16))
    for rho in rhos:
        hm_wd = Integrator(hm_gamma, hm_K, r_init, m_init, rho)
        new_row = {'Rho': rho,
                   'Mass': max(hm_wd['Mass']),
                   'Radius': max(hm_wd['Radius'])
                   }
        out.append(new_row)

    # To be turned into the table for problem 2
    high_mass_white_dwarfs = pd.DataFrame(out)
    return high_mass_white_dwarfs


# CALL STUFF
# Problem 1
lm_white_dwarf = Integrator(lm_gamma, lm_K, r_init, m_init, rho_init)
m_lmwd = max(lm_white_dwarf['Mass']) * kg2solarMass
r_lmwd = max(lm_white_dwarf['Radius']) * m2solarRadius



# Problem 2
hm_white_dwarf = Integrator(hm_gamma, hm_K, r_init, m_init, rho_init)
density_vals = mass_rad_relation()




# The End
