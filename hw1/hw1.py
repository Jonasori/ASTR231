"""
HW1 Support Functions
Stellar Physics
Jonas Powell
February 11, 2018
"""


# Problem 3
"""
Photometric data on a star are given in the table below. For your convenience, I also provide the flux density of a zero magnitude star and effective wavelengths of the filters for the photometric system in which the star was observed. The flux densities are in Jansky's (Jy), a favorite unit of infrared and radio astronomers. Be careful, because the Jy is a unit of F$_\nu$, not F$_\lambda$. There is a link on the course Moodle page to a Web site that will help you do the conversions correctly.

Plot the tabulated data on a $F_\lambda$ versus log ($\lambda$) diagram. Make the scale be W/m$^2$/micron, as on the example, with wavelength expressed in microns. On the same diagram plot the blackbody curve (Planck function) with the temperature that you think best fits the data. In other words, you are using the spectral energy distribution (fondly called the SED by most astronomers) to estimate the effective temperature of the star. Be sure to state the effective temperature that you derive for the star by this method.
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const

h = const.h.value
c = const.c.value
k = const.k_B.value
print k

# Arbitrary for now, in K (obvs)
T = 6000

def planck_lam(lam):
    B_nu = (2 * h * c**2 / lam**5) / (np.exp(h*c / (lam * k * T)) - 1)
    return B_nu

planck_lam(10)
# The End

lams = np.arange(200, 1200)

plt.plot(lams, planck_lam(lams))
plt.show()
