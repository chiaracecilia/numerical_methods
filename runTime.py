# calculate the running time of the implemented numerical schemeComparison
# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *
import timeit

def runTime():

    # parameters
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
    phiOld = cosBell(x, alpha, beta)
    # Exact solution is the initial condition shifted around the domain
    phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), alpha, beta)

    # Advect the profile using finite difference for all the time steps
    startFTBS = timeit.default_timer()
    phiFTBS, TV_FTBS, numerical_mean_FTBS  = FTBS(phiOld.copy(), c, nt)
    stopFTBS = timeit.default_timer()
    startCTCS = timeit.default_timer()
    phiCTCS, TV_CTCS, numerical_mean_CTCS = CTCS(phiOld.copy(), c, nt, False)
    stopCTCS = timeit.default_timer()
    startBTCS = timeit.default_timer()
    phiBTCS, numerical_mean_BTCS = BTCS(phiOld.copy(), c, nt)
    stopBTCS = timeit.default_timer()
    startLAGR = timeit.default_timer()
    phiLAGR, TV_SL, numerical_mean_SL = SemiLagrangian(phiOld.copy(), c, nt, dx, x)
    stopLAGR = timeit.default_timer()

    print('Time FTBS: ', stopFTBS - startFTBS, 'Time CTCS: ', stopCTCS - startCTCS,\
     'Time BTCS: ', stopBTCS - startBTCS, 'Time LAGR: ', stopLAGR - startLAGR)
