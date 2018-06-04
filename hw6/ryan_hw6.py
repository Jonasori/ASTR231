from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

h = 6.62607004e-34 * 1e7 #erg*s
m_e = 9.10938356e-31 * 1e3 #g
m_p = 1.6726219e-27 * 1e3#g
pi = np.pi
G = 6.67408e-11 *(1e2)**3 * (1./1e3) #m^3 kg-1 s-2
K = (h**2. / (5.*m_e))*(3./(8.*pi))**(2./3.) * (1.0/(2.0*m_p))**(5./3.)
rho_0=1e5 #g/cm^3
m_sun = 1.989e30 * 1e3 #g
r_sun = 695508000. *1e2 # cm
r_0 = 1e-6
t = np.linspace(r_0, 1e12, 10000000)
true_r = 1.4525e11  #cm
true_M = .5*m_sun


# There's no way this is right
def dmdr_dpdr(y, t, K):
	M_r, P = y
	rho = (P / K) ** (3./5.)
	# These should be rs not ts
	dpdr= -(G * M_r * rho) / (t**2)
	dmdr = 4.0 * pi * t**2 * rho

	#print dmdr,dpdr
	#print ''
	return [dmdr,dpdr]


def integrate(args):
	K, rho_0 = args
	sol = odeint(dmdr_dpdr,[rho_0*(4./3.)*pi*r_0**3,K*(rho_0**(5./3.))],t,args=(K,))

	max_M_so_far=-float('inf')
	r_at_max_M=-float('inf')

	for i in range(len(sol)):
		s=sol[i]
		if s[0]>max_M_so_far:
			max_M_so_far=s[0]
			r_at_max_M=t[i]
	score = ((max_M_so_far-true_M)/m_sun)**2 + ((r_at_max_M-true_r)/r_sun)**2
	print score,K,rho_0
	return score
	#print 'Mass: '+str(max_M_so_far/m_sun)+' solar masses'
	#print 'Radius: '+str(r_at_max_M/r_sun)+' solar radii'
	#d = max_M_so_far/((4./3.)*pi*r_at_max_M**3)
	#print 'Density: '+str(d)



res = minimize(integrate,[K,rho_0],method = 'Nelder-Mead')

K = res.x[0]
rho_0 = res.x[1]

sol = odeint(dmdr_dpdr,[rho_0*(4./3.)*pi*r_0**3,K*(rho_0**(5./3.))],t,args=(K,))

max_M_so_far=-float('inf')
r_at_max_M=-float('inf')

for i in range(len(sol)):
	s=sol[i]
	if s[0]>max_M_so_far:
		max_M_so_far=s[0]
		r_at_max_M=t[i]

print 'Mass: '+str(max_M_so_far/m_sun)+' solar masses'
print 'Radius: '+str(r_at_max_M/r_sun)+' solar radii'
d = max_M_so_far/((4./3.)*pi*r_at_max_M**3)
print 'Density: '+str(d)



#r = 1.4525e9 #meters
