#!/usr/bin/python3

# Outer code for setting up the linear advection problem on a uniform
# grid and calling the function to perform the linear advection and plot.

import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

#plot the different implementation of FTBS, CTCS, BTCS and SemiLAGRANGIAN advection schemes on the same graph
#compared to the initial solution.

#print the error norm
#cosinebell initial condition.
def scheme_comparison(init_cond,name):
    "Advect the initial conditions using various advection schemes and"
    "compare results"

    # Parameters
    xmin = 0
    xmax = 1
    nx = 100
    nt = 40
    c = 0.2
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
    phiFTBS, TV_FTBS, numerical_mean_FTBS  = FTBS(phiOld.copy(), c, nt)
    phiCTCS, TV_CTCS, numerical_mean_CTCS = CTCS(phiOld.copy(), c, nt, False)
    phiBTCS, numerical_mean_BTCS = BTCS(phiOld.copy(), c, nt)
    phiLAGR, TV_SL, numerical_mean_SL = SemiLagrangian(phiOld.copy(), c, nt, dx, x)
    phiCTCS_clipping, TV_CTCS_clipping, numerical_mean_CTCS_clipping  = CTCS(phiOld.copy(), c, nt, True)

    # Calculate and print out error norms
    print("FTBS l2 error norm = ", l2ErrorNorm(phiFTBS, phiAnalytic), "CTCS l2 error norm = ", l2ErrorNorm(phiCTCS, phiAnalytic), "BTCS l2 error norm = ", l2ErrorNorm(phiBTCS, phiAnalytic), "Semi-Lagrangian l2 error norm = ", l2ErrorNorm(phiLAGR, phiAnalytic) )
    print("FTBS linf error norm = ", lInfErrorNorm(phiFTBS, phiAnalytic), "CTCS linf error norm = ", lInfErrorNorm(phiCTCS, phiAnalytic), "BTCS linf error norm = ", lInfErrorNorm(phiBTCS, phiAnalytic),"FTBS linf error norm = ", lInfErrorNorm(phiLAGR, phiAnalytic) )

    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)

    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='--', linewidth=2)
    plt.plot(x, phiFTBS, label='FTBS', color='blue')
    plt.plot(x, phiCTCS, label='CTCS', color='red')
    plt.plot(x, phiBTCS, label='BTCS', color='green')
    plt.plot(x, phiLAGR, label='Lagrangian', color='yellow')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])
    plt.legend()
    plt.xlabel('$x$')
    plt.savefig('plots/'+name+'.pdf')

    plt.figure(2)
    plt.clf()
    plt.ion()
    plt.plot(t, TV_FTBS, ".-",label='Total Variation FTBS', color='red')
    plt.plot(t, TV_CTCS, ".-", label='Total Variation CTCS', color='green')
    plt.plot(t, TV_SL, ".-", label='Total Variation SL', color='blue')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([0,20])
    plt.legend()
    plt.xlabel('$x$')
    plt.savefig('plots/'+name+'_total_variation.pdf')

    plt.figure(3)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='--', linewidth=2)
    plt.plot(x, phiCTCS, label='CTCS', color='aqua')
    plt.plot(x, phiCTCS_clipping, label='CTCS with clipping', color='blue')

    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])
    plt.legend()
    plt.xlabel('$x$')
    plt.savefig('plots/'+name+'_clipping.pdf')

scheme_comparison(cosBell, 'cosine_belle')
scheme_comparison(squareWave, 'square_wave')
