"""
Develop a best-fit model for a T-Tauri star with given parameters.

Using scipy's integrate.odeint and optimize.minimize features, do some stuff.
"""

# Packages
from scipy.integrate import odeint
from scipy.optimize import minimize
import numpy as np
import astropy.constants as c
import matplotlib.pyplot as plt


# Some constants
G = c.G.value                       # m3 kg-1 s-2
sig = c.sigma_sb.value              # W / (K4 m2)
m_sun = c.M_sun.value               # kg
r_sun = c.R_sun.value               # m

mmw = 0.6                           # kg
gamma = 5./3.
l_star = 3.848e28                   # mks whatevers
t_star = 4000                       # Kelvin
m_true = 0.5 * m_sun                # kg
r_true = np.sqrt(l_star/(4*np.pi * sig * t_star**4))

# Some initial values
K_guess = 4.1014e8
rho0_guess = 464.097
m0 = 0
P0_guess = K_guess * rho0_guess**gamma
rs = np.linspace(1, 1e12, int(1e6))

# Set up a list to store visits in K/rho0 space (to be used for plotting)
visited_vals = []

# A list to store the steps through M, P in the integration
steps = []


# The differentials
def differentials(conditions, r, K):
    """Calulate the differential changes in mass and pressure."""
    m, P = conditions

    rho = (P/K)**(1./gamma)
    dM = 4. * np.pi * rho * r**2.
    dP = -G * m * rho * r**(-2.)

    return [dM, dP]


def integrator(vals):
    """Get a chi2 value comparing a model to real values.

    Using given values of K and rho0, make a model star and compare the
    resulting mass and radius to the known correct values.
    """
    K, rho0 = vals

    # Set up a dict to store the visited values
    # visited_vals = {''}

    # Make a model star
    P0_guess = K * rho0**gamma
    sol = odeint(func=differentials, y0=[m0, P0_guess], t=rs, args=(K,))

    # Calculate the chi-squared of this combination
    # Normalize by solar values to get fractional correctness
    ms = np.transpose(sol)[0]
    idx_of_max_m = np.where(ms == max(ms))[0][0]

    m_max = ms[idx_of_max_m]
    r_of_max_m = rs[idx_of_max_m]

    chi2 = ((m_max - m_true)/m_sun)**2 + ((r_of_max_m - r_true)/r_sun)**2
    print K, rho0

    visited_vals.append([K, rho0])

    return chi2


rs = np.linspace(1, 1e12, int(1e6))
r_locs = rs/r_sun


def optimize_fit_scipy(make_plots=False):
    """Search K/rho0 parameter space for best fit with scipy."""
    rs = np.linspace(1, 1e12, int(1e6))

    # Find the best fit:
    guessed_vals = [K_guess, rho0_guess]
    min_vals = minimize(fun=integrator, x0=guessed_vals, method='Nelder-Mead')

    # Interpret the results:
    best_chi2 = min_vals.fun
    best_K, best_rho0 = min_vals.x[0], min_vals.x[1]

    print best_K
    # Make a model with those best fit vals
    P0 = best_K * best_rho0**gamma
    final = odeint(func=differentials, y0=[m0, P0], t=rs, args=(best_K,))

    ms = np.transpose(final)[0]
    idx_of_max_m = np.where(ms == max(ms))[0][0]
    final_m = ms[idx_of_max_m]
    final_r = rs[idx_of_max_m]

    print "Mass of best-fit star (solar masses): ", final_m/m_sun
    print "Radius of best-fit star (solar radii): ", final_r/r_sun
    print "Best chi2 value: ", best_chi2
    print "Final K value: ", best_K
    print "Final rho0 value: ", best_rho0

    pressures = np.transpose(final)[1]

    rhos = []; [rhos.append(best_K * pressures[i] ** (gamma)) for i in range(len(pressures))]

    temps = []; [temps.append(rhos[i] * mmw/(pressures[i] * best_K)) for i in range(len(pressures))]


    if make_plots:
        # Rs already defined

        plt.plot(rs, rhos, color='red', alpha=0.6, marker='.')
        plt.title('Density')
        plt.xlabel('Stellar radius (meters)')
        plt.ylabel('Density (kg/m^3)')
        #plt.xticks(r_locs)
        plt.savefig('rho_r.png')
        plt.close()

        plt.plot(rs, pressures, color='blue', alpha=0.3, marker='.')
        plt.title('Pressure')
        plt.xlabel('Stellar radius (meters)')
        plt.ylabel('Pressure (Pascals)')
        plt.savefig('pressure_r.png')
        plt.close()

        plt.plot(rs, temps, color='green', alpha=0.3, marker='.')
        plt.title('Temperature')
        plt.xlabel('Stellar radius (meters)')
        plt.ylabel('Temperature (Kelvin)')
        plt.savefig('Temperature_r.png')
        plt.close()

    return [pressures, rhos, temps]


def optimize_fit_emcee():
    """Search K/rho0 parameter space for best fit with emcee."""
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

        [m, r] = structure_integrator(gamma, K, r0, m0, rho0)
        chisq = chi2([m, r], [m_star, r_star])
        ln_p = -0.5 * chisq

        print "New chi2 val: ", chisq
        return ln_p


    ndim, nwalkers = 2, 10
    # Give initial positions in [K, rho0] space for the walkers
    p0 = [np.random.rand(ndim) for i in range(nwalkers)]

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
    sampler.run_mcmc(p0, 100)

    return -1


# The End
