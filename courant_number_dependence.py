
import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

# BTCS is unconditionally stable
# In order to show this we will implement the numerical method with a very big Courant Number,
# and show that this does not affect any characteristic of the numerical schemes

#We will also show the behaviour of the other numerical schemes

def scheme_comparison(init_cond, C):
    "Advect the initial conditions using various advection schemes and"
    "compare results"

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
    phiFTBS, TV_FTBS, numerical_mean_FTBS  = FTBS(phiOld.copy(), c, nt)
    phiCTCS, TV_CTCS, numerical_mean_CTCS = CTCS(phiOld.copy(), c, nt, False)
    phiBTCS, numerical_mean_BTCS = BTCS(phiOld.copy(), c, nt)
    phiLAGR, TV_SL, numerical_mean_SL = SemiLagrangian(phiOld.copy(), c, nt, dx, x)

    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)

    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black',
             linestyle='--', linewidth=2)
    plt.plot(x, phiBTCS, label='BTCS with courant number %d' %c, color='blue')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])
    plt.legend()
    plt.xlabel('$x$')
    plt.savefig('plots/BTCS_unconditionally_stable.pdf')


    plt.figure(2)
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
    plt.savefig('plots/dependence_courant_number.pdf')
scheme_comparison(cosBell, 2)
