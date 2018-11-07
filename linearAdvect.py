#!/usr/bin/python3

# Outer code for setting up the linear advection problem on a uniform
# grid and calling the function to perform the linear advection and plot.

import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

# def main():
#     "Advect the initial conditions using various advection schemes and"
#     "compare results"
#
#     # Parameters
#     xmin = 0
#     xmax = 1
#     nx = 40
#     nt = 40
#     c = 0.2
#
#     # Derived parameters
#     dx = (xmax - xmin)/nx
#
#     # spatial points for plotting and for defining initial conditions
#     x = np.arange(xmin, xmax, dx)
#
#     # Initial conditions
#     phiOld = cosBell(x, 0, 0.75)
#     # Exact solution is the initial condition shifted around the domain
#     phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
#
#     # Advect the profile using finite difference for all the time steps
#     phiCTCS = CTCS(phiOld.copy(), c, nt)
#
#     # Calculate and print out error norms
#     print("CTCS l2 error norm = ", l2ErrorNorm(phiCTCS, phiAnalytic))
#     print("# COMBAK: TCS linf error norm = ", lInfErrorNorm(phiCTCS, phiAnalytic))
#
#     # Plot the solutions
#     font = {'size'   : 20}
#     plt.rc('font', **font)
#     plt.figure(1)
#     plt.clf()
#     plt.ion()
#     plt.plot(x, phiOld, label='Initial', color='black')
#     plt.plot(x, phiAnalytic, label='Analytic', color='black',
#              linestyle='--', linewidth=2)
#     plt.plot(x, phiCTCS, label='CTCS', color='blue')
#     plt.axhline(0, linestyle=':', color='black')
#     plt.ylim([-0.2,1.2])
#     plt.legend(bbox_to_anchor=(1.15 , 1.1))
#     plt.xlabel('$x$')
#     input('press return to save file and continue')
#     plt.savefig('plots/CTCS.pdf')
#
# main()

#!/usr/bin/python3

# Outer code for setting up the linear advection problem on a uniform
# grid and calling the function to perform the linear advection and plot.




# def main():
#     "Advect the initial conditions using various advection schemes and"
#     "compare results"
#
#     # Parameters
#     xmin = 0
#     xmax = 1
#     nx = 40
#     nt = 40
#     c = 0.2
#     #parameters of the square wave
#     alpha=0
#     beta=0.2
#     # Derived parameters
#     dx = (xmax - xmin)/nx
#
#     # spatial points for plotting and for defining initial conditions
#     x = np.arange(xmin, xmax, dx)
#
#     # Initial conditions
#     phiOld = squareWave(x, alpha, beta)
#     # Exact solution is the initial condition shifted around the domain
#     phiAnalytic = squareWave((x - c*nt*dx)%(xmax - xmin), alpha, 0.75)
#
#     # Advect the profile using finite difference for all the time steps
#     phiCTCS = CTCS(phiOld.copy(), c, nt)
#
#     # Calculate and print out error norms
#     print("CTCS l2 error norm = ", l2ErrorNorm(phiCTCS, phiAnalytic))
#     print("# COMBAK: TCS linf error norm = ", lInfErrorNorm(phiCTCS, phiAnalytic))
#
#     # Plot the solutions
#     font = {'size'   : 20}
#     plt.rc('font', **font)
#     plt.figure(1)
#     plt.clf()
#     plt.ion()
#     plt.plot(x, phiOld, label='Initial', color='black')
#     plt.plot(x, phiAnalytic, label='Analytic', color='black',
#              linestyle='--', linewidth=2)
#     plt.plot(x, phiCTCS, label='CTCS', color='blue')
#     plt.axhline(0, linestyle=':', color='black')
#     plt.ylim([-0.2,1.2])
#     plt.legend(bbox_to_anchor=(1.15 , 1.1))
#     plt.xlabel('$x$')
#     input('press return to save file and continue')
#     plt.savefig('plots/CTCS_SQUAREWAVE.pdf')
#
# main()

# def main():
#     "Advect the initial conditions using various advection schemes and"
#     "compare results"
#
#     # Parameters
#     xmin = 0
#     xmax = 1
#     nx = 40
#     nt = 40
#     c = 0.2
#     #parameters of the square wave
#     alpha=0
#     beta=0.2
#     # Derived parameters
#     dx = (xmax - xmin)/nx
#
#     # spatial points for plotting and for defining initial conditions
#     x = np.arange(xmin, xmax, dx)
#
#     # Initial conditions
#     phiOld = squareWave(x, alpha, beta)
#     # Exact solution is the initial condition shifted around the domain
#     phiAnalytic = squareWave((x - c*nt*dx)%(xmax - xmin), alpha, beta)
#
#     # Advect the profile using finite difference for all the time steps
#     phiBTCS = BTCS(phiOld.copy(), c, nt)
#
#     # Calculate and print out error norms
#     print("BTCS l2 error norm = ", l2ErrorNorm(phiBTCS, phiAnalytic))
#     print("# COMBAK: TCS linf error norm = ", lInfErrorNorm(phiBTCS, phiAnalytic))
#
#     # Plot the solutions
#     font = {'size'   : 20}
#     plt.rc('font', **font)
#     plt.figure(1)
#     plt.clf()
#     plt.ion()
#     plt.plot(x, phiOld, label='Initial', color='black')
#     plt.plot(x, phiAnalytic, label='Analytic', color='black',
#              linestyle='--', linewidth=2)
#     plt.plot(x, phiBTCS, label='BTCS', color='blue')
#     plt.axhline(0, linestyle=':', color='black')
#     plt.ylim([-0.2,1.2])
#     plt.legend(bbox_to_anchor=(1.15 , 1.1))
#     plt.xlabel('$x$')
#     input('press return to save file and continue')
#     plt.savefig('plots/BTCS_SQUAREWAVE.pdf')
#
# main()

# def main():
#     "Advect the initial conditions using various advection schemes and"
#     "compare results"
#
#     # Parameters
#     xmin = 0
#     xmax = 1
#     nx = 40
#     nt = 40
#     c = 0.2
#     #parameters of the square wave
#     alpha=0
#     beta=0.2
#     # Derived parameters
#     dx = (xmax - xmin)/nx
#
#     # spatial points for plotting and for defining initial conditions
#     x = np.arange(xmin, xmax, dx)
#
#     # Initial conditions
#     phiOld = cosBell(x, alpha, beta)
#     # Exact solution is the initial condition shifted around the domain
#     phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), alpha, beta)
#
#     # Advect the profile using finite difference for all the time steps
#     phiBTCS = BTCS(phiOld.copy(), c, nt)
#
#     # Calculate and print out error norms
#     print("BTCS l2 error norm = ", l2ErrorNorm(phiBTCS, phiAnalytic))
#     print("# COMBAK: TCS linf error norm = ", lInfErrorNorm(phiBTCS, phiAnalytic))
#
#     # Plot the solutions
#     font = {'size'   : 20}
#     plt.rc('font', **font)
#     plt.figure(1)
#     plt.clf()
#     plt.ion()
#     plt.plot(x, phiOld, label='Initial', color='black')
#     plt.plot(x, phiAnalytic, label='Analytic', color='black',
#              linestyle='--', linewidth=2)
#     plt.plot(x, phiBTCS, label='BTCS', color='blue')
#     plt.axhline(0, linestyle=':', color='black')
#     plt.ylim([-0.2,1.2])
#     plt.legend(bbox_to_anchor=(1.15 , 1.1))
#     plt.xlabel('$x$')
#     input('press return to save file and continue')
#     plt.savefig('plots/BTCS.pdf')
#
# main()

# def main():
#     "Advect the initial conditions using various advection schemes and"
#     "compare results"
#
#     # Parameters
#     xmin = 0
#     xmax = 1
#     nx = 40
#     nt = 40
#     c = 0.2
#     #parameters of the square wave
#     alpha=0
#     beta=0.75
#     # Derived parameters
#     dx = (xmax - xmin)/nx
#
#     # spatial points for plotting and for defining initial conditions
#     x = np.arange(xmin, xmax, dx)
#
#     # Initial conditions
#     phiOld = cosBell(x, alpha, beta)
#     # Exact solution is the initial condition shifted around the domain
#     phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), alpha, beta)
#
#     # Advect the profile using finite difference for all the time steps
#     phiFTBS = FTBS(phiOld.copy(), c, nt)
#
#     # Calculate and print out error norms
#     print("FTBS l2 error norm = ", l2ErrorNorm(phiFTBS, phiAnalytic))
#     print("# COMBAK: TCS linf error norm = ", lInfErrorNorm(phiFTBS, phiAnalytic))
#
#     # Plot the solutions
#     font = {'size'   : 20}
#     plt.rc('font', **font)
#     plt.figure(1)
#     plt.clf()
#     plt.ion()
#     plt.plot(x, phiOld, label='Initial', color='black')
#     plt.plot(x, phiAnalytic, label='Analytic', color='black',
#              linestyle='--', linewidth=2)
#     plt.plot(x, phiFTBS, label='FTBS', color='blue')
#     plt.axhline(0, linestyle=':', color='black')
#     plt.ylim([-0.2,1.2])
#     plt.legend(bbox_to_anchor=(1.15 , 1.1))
#     plt.xlabel('$x$')
#     input('press return to save file and continue')
#     plt.savefig('plots/FTBS.pdf')
#
# main()

def main():
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

    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)

    # Initial conditions
    phiOld = cosBell(x, alpha, beta)
    # Exact solution is the initial condition shifted around the domain
    phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), alpha, beta)

    # Advect the profile using finite difference for all the time steps
    phiFTBS = FTBS(phiOld.copy(), c, nt)
    phiCTCS = CTCS(phiOld.copy(), c, nt)
    phiBTCS = BTCS(phiOld.copy(), c, nt)

    # Calculate and print out error norms
    print("FTBS l2 error norm = ", l2ErrorNorm(phiFTBS, phiAnalytic), "CTCS l2 error norm = ", l2ErrorNorm(phiCTCS, phiAnalytic), "BTCS l2 error norm = ", l2ErrorNorm(phiBTCS, phiAnalytic) )
    print("FTBS linf error norm = ", lInfErrorNorm(phiFTBS, phiAnalytic), "CTCS linf error norm = ", lInfErrorNorm(phiCTCS, phiAnalytic), "BTCS linf error norm = ", lInfErrorNorm(phiBTCS, phiAnalytic), )

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
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])
    plt.legend()
    plt.xlabel('$x$')
    input('press return to save file and continue')
    plt.savefig('plots/comparison_cosBell.pdf')

main()
