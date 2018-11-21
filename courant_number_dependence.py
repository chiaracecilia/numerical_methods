# BTCS is unconditionally stable
# In order to show this we will implement the numerical method with a very big Courant Number,
# and show that this does not affect the stability of the numerical scheme
#We will also show the behaviour of the other numerical schemes

import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

def courant_number_dependence(init_cond, C):
    "Implement the schemes with the initial condition and courant number"
    "given in input. Calling the function will result in two plots for the"
    "implementation of the scheme with different courant number"

    # Parameters
    xmin = 0
    xmax = 1
    nx = 100
    nt = 40
    c = C
    #parameters of the square wave
    alpha=0
    beta=0.75
    # Derived parameters
    dx = (xmax - xmin)/nx

    # spatial and temporal points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)
    t = np.arange(nt+1)
    # Initial conditions
    phiOld = init_cond(x, alpha, beta)
    # Exact solution is the initial condition shifted around the domain
    phiAnalytic = init_cond((x - c*nt*dx)%(xmax - xmin), alpha, beta)

    # Advect the profile using finite difference for all the time steps
    phiFTBS, TV_FTBS, numericalMassFTBS  = FTBS(phiOld.copy(), c, nt)
    phiCTCS, TV_CTCS, numericalMassCTCS = CTCS(phiOld.copy(), c, nt, False)
    phiBTCS, numericalMassBTCS = BTCS(phiOld.copy(), c, nt)
    phiLAGR, TV_SL, numericalMassSL = SemiLagrangian(phiOld.copy(), c, nt, dx, x)

    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)

    # plot the advection with BTCS
    plt.figure(1, figsize=(10,7))
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='--', linewidth=2)
    plt.plot(x, phiBTCS, label='BTCS with courant number %d' %c, color='blue')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.5,1.2])
    plt.legend()
    plt.xlabel('$x$')
    plt.savefig('plots/BTCS_unconditionally_stable.pdf')

    #plot the advection for all the numerical schemes
    plt.figure(2, figsize=(10,7))
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='--', linewidth=2)
    plt.plot(x, phiFTBS, label='FTBS', color='blue')
    plt.plot(x, phiCTCS, label='CTCS', color='aqua')
    plt.plot(x, phiBTCS, label='BTCS', color='green')
    plt.plot(x, phiLAGR, label='Lagrangian', color='navy')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.5,1.2])
    plt.legend()
    plt.xlabel('$x$')
    plt.savefig('plots/dependence_courant_number.pdf')
