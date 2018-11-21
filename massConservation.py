# Calculating the mass for the different schemes and plot it against time
# to check whether it is conserved or not
import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

def mass(initCond):
    "this function will produce a plot of mass against time for different schemeComparison"
    "and check the conservation of mass for CTCS and clipped CTCS, proving that"
    "CTCS with clipping will not conserve this quantity"
    # Parameters
    xmin = 0
    xmax = 1
    nx = 15
    nt = 10
    c = 0.2
    #parameters of the square wave
    alpha=0
    beta=0.75
    # Derived parameters
    dx = (xmax - xmin)/nx
    t = np.arange(nt+1)

    # spatial and temporal points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)
    # Initial conditions
    phiOld = initCond(x, alpha, beta)
    # calculate the exact mass
    exactMass = np.sum(phiOld)
    # apply the different schemes on the initial condition
    phiFTBS, TV_FTBS, numericalMassFTBS  = FTBS(phiOld.copy(), c, nt)
    phiCTCS, TV_CTCS, numericalMassCTCS  = CTCS(phiOld.copy(), c, nt, False)
    phiCTCS_clipping, TV_CTCS_clipping, numericalMassCTCSClipping  = CTCS(phiOld.copy(), c, nt, True)
    phiBTCS, numericalMassBTCS  = BTCS(phiOld.copy(), c, nt)
    phiLS, TV_SL, numericalMassSL  = SemiLagrangian(phiOld.copy(), c, nt, dx, x)

    # here we plot the first graph where we confront the conservation of mass
    # for all the numerical schemes
    plt.figure(1, figsize = (10, 7))
    plt.clf()
    plt.ion()
    plt.plot(t, numericalMassFTBS, ".-",label='FTBS', color='aqua')
    plt.plot(t, numericalMassCTCS, ".-",label='CTCS', color='red')
    plt.plot(t, numericalMassCTCSClipping, ".-",label='CTCS with clipping', color='magenta')
    plt.plot(t, numericalMassBTCS, ".-",label='BTCS', color='blue')
    plt.plot(t, numericalMassSL, ".-",label='SemiLagrangian', color='navy')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([10,12])
    plt.legend()
    plt.xlabel('$t$')
    plt.ylabel('$Mass$')
    plt.savefig('plots/conservation_of_mass.pdf')

    # here we plot the conservation of mass for CTCS with or without clipping
    plt.figure(2, figsize = (10, 7))
    plt.clf()
    plt.ion()
    plt.plot(t, numericalMassCTCS, ".-",label='CTCS', color='red')
    plt.plot(t, numericalMassCTCSClipping, ".-",label='CTCS with clipping', color='magenta')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([10,12])
    plt.legend()
    plt.xlabel('$t$')
    plt.ylabel('$Mass$')
    plt.savefig('plots/conservation_of_mass_clipping.pdf')

    # print the value of the exact mass top see if it conserved
    print("the exact mass is: ", exactMass)
