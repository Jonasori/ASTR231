import numpy as np
import matplotlib.pyplot as plt

###### QUESTION 1 #######

#dm(r)/dr = 4 * np.pi * (r**2) * rho(r)    #Mass of a spherical shell
#dP(r)/dr = (-G * m(r) * rho(r))/(r**2)    #Hydrostatic equilibrium eqn

#Constants and given values:

G = 6.67e-11 #m^3 kg^-1 s^-2
m_sun = 2.e30 #kg
r_sun = 6.9e8 #meters
r_earth = 6.37e6 #meters
h = 6.63e-34 #m**2 kg/s
m_e = 9.1e-31 #kg
m_h = 1.67e-27 #kg
n = 3./2.
gamma = 5./3.

K = ((h**2.)/(5.*m_e))*((3./(8.*np.pi))**(2./3.)) * (1./(2.*m_h))**(5./3.)  #for non-rel, white dwarf case

K

r = 0.000001 #initial r should be close to zero at center of star
m = 0.  #initial mass at center of star should be zero
rho_0 = 1e8 #kg/m**3

#Eqn of state function:
def EOS(gamma,rho):
    return K * rho ** gamma

#eqn for rho:
def rho_func(gamma,P):
    return (P/K) ** (1./gamma)

#hydrostatic eq function:
def hydroeq(m, gamma, P, r):
    rho=rho_func(gamma,P)
    dP = (-G * m * rho)/r**(2.)
    return dP

#mass of spherical shell function:
def mass(r, gamma, P, rho):
    rho=rho_func(gamma,P)
    dm = 4. * np.pi * (r**2.) * rho
    return dm

P = EOS(5./3.,1e8)  #initial pressure should be the max at the center of the star
print "P initial:", P

##arrays to fill with info:
a_r = []  # empty lists that will be filled with the info from the while loop
a_m = []
a_P = []
a_rho = []

#integrating:

def Integrate(r,m,P,rho,n,gamma):
    dr = 10.
    while P>0:          #knows to integrate dm and dP because they are in the eqns:
        r_iplus1 = r + dr
        m_iplus1 = m + (dr * mass(r,gamma, P, rho))
        P_iplus1 = P + (dr * hydroeq(m,gamma, P,r))
        rho_iplus1 = rho + (dr*rho_func(gamma,P))

        a_r.append(r_iplus1)
        a_m.append(m_iplus1)
        a_P.append(P_iplus1)
        a_rho.append(rho_iplus1)

        r = r_iplus1
        m = m_iplus1
        P = P_iplus1
        rho = rho_iplus1
    return np.array(a_r), np.array(a_m), np.array(a_P), np.array(a_rho)
### What is happening: pressure starts at one, as while loop runs the pressure
### will go down in steps of dr util pressure is zero

radius1, mass1, pressure1, rho1 = Integrate(r,m,P,rho_0,n,gamma) #for n = 3/2: gives you  bunch of each rad, mass, pressure as written in the return statement
print "max mass in terms of solar mass:", np.max(mass1)/m_sun

mean_rho = np.max(mass1)/((4./3.) * np.pi * np.max(radius1)**3)
print "average density in kg/m^3:", mean_rho

rad_earthterms= np.max(radius1)/r_earth
rad_sunterms= np.max(radius1)/r_sun
print "radius in meters:", np.max(radius1)
print "radius in terms of earth radii:", rad_earthterms
print "radius in terms of solar radii:", rad_sunterms

print

###### QUESTION 2 #######

#new constants:
c = 3e8 #meters per second
gamma_2 = 4./3.
n_2 = 3.
K_2 = (h*c/4) * (3/(8*np.pi))**(1./3.) * (1./(2.*m_h))**(4./3.)  #ultra relativistic case
rho_2 = 1e12 #kg/m^3
rho_mid = 1e15
rho_3 = 1e18 #kg/m^3

### New funcitons:
#Eqn of state function:
def EOS2(gamma,rho):
    return K_2 * rho ** gamma

#eqn for rho:
def rho_func2(gamma,P):
    return (P/K_2) ** (1./gamma)

#hydrostatic eq function:
def hydroeq2(m, gamma, P, r):
    rho=rho_func2(gamma,P)
    return (-G * m * rho)/r**(2.)

#mass of spherical shell function:
def mass2(r, gamma, P, rho):
    rho=rho_func2(gamma,P)
    return 4. * np.pi * (r**2.) * rho

P_2 = EOS2(4./3.,1e12) #of gamma 2 and rho 2
#P_mid= EOS2(4./3.,1e15) #of gamma 2 and rho 2
#P_3 = EOS2(4./3.,1e18) #of gamma 2 and rho 2

print "P initial, relativistic:", P_2

#integrating:

def Integrate2(r,m,P,rho,gamma):
    dr = 1000.
    while P>0:          #knows to integrate dm and dP because they are in the eqns:
        r_iplus1 = r + dr
        m_iplus1 = m + (dr * mass2(r,gamma, P, rho))
        P_iplus1 = P + (dr * hydroeq2(m,gamma, P,r))
        rho_iplus1 = rho + (dr*rho_func2(gamma,P))

        a_r.append(r_iplus1)
        a_m.append(m_iplus1)
        a_P.append(P_iplus1)
        a_rho.append(rho_iplus1)

        r = r_iplus1
        m = m_iplus1
        P = P_iplus1
        rho = rho_iplus1
    return np.array(a_r), np.array(a_m), np.array(a_P), np.array(a_rho)
### What is happening: pressure starts at one, as while loop runs the pressure
### will go down in steps of dr util pressure is zero

radius2, mass2, pressure2, rho2 = Integrate2(r,m,P_2,rho_2,gamma_2)
#radiusmid, massmid, pressuremid, rhomid = Integrate2(r,m,P_mid,rho_mid,gamma_2)
#radius3, mass3, pressure3, rho3 = Integrate2(r,m,P_3,rho_3,gamma_2)

print "max mass for relativistic/ chandrasekhar limit:", np.max(mass2)

####### ANSWERS: ########
'''
    Question 1:

P initial: 6.76123159365e+19
max mass in terms of solar mass: 0.154434404249
average density in kg/m^3: 16692908.019
radius in meters: 16407830.0
radius in terms of earth radii: 2.57579748823
radius in terms of solar radii: 0.0237794637681

    Question 2:

P initial, relativistic: 3.13828570665e+22
max mass for relativistic/ chandrasekhar limit (kg): 2.87129255231e+30
'''

###########################




#Plotting for prob 2:
'''
#radius - mass:
plt.plot(radius2,mass2)
plt.xlabel('Radius [m]')
plt.ylabel('Mass [kg]')
plt.title('     Mass-radius relation for n=3, gamma=4/3 white dwarf')
plt.savefig('stellarhw5_m_r.png')
plt.close()

#radius density:
plt.plot(radius2,rho2, label='rho = 10^12 kg/m^3 (low)')
#plt.plot(radiusmid,rhomid, label='rho = 10^15 kg/m^3 (mid)')
#plt.plot(radius3,rho3, label='rho = 10^18 kg/m^3 (high)')
plt.xlabel('Radius [m]')
plt.ylabel('Density [kg/m^3]')
plt.title('Denity-radius relation for n=3, gamma=4/3 white dwarf')
plt.legend(loc='best')
plt.savefig('stellarhw5_d_r.png')

plt.close()

#radius pressure:
plt.plot(radius2,pressure2)
plt.xlabel('Radius [m]')
plt.ylabel('Pressure [kg/m^2]')
plt.title('Pressure-radius relation for n=3, gamma=4/3 white dwarf')
plt.savefig('stellarhw5_p_r.png')
plt.close()

'''
