# produces plots to assess the order of accuracy
import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

def accuracy():
    "produce plot to assess the order of accuracy"

    # Parameters
    xmin = 0
    xmax = 1
    nx = 40
    nt = 40
    c = 0.5

    # Derived parameters
    dx = (xmax - xmin)/nx
    # Number of different size count we will consider
    nDx = 10

    # spatial and temporal points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)

    # set up the vectors that will store the error for each method
    errorFTBS = np.zeros(nDx)
    errorBTCS = np.zeros(nDx)
    errorCTCS = np.zeros(nDx)
    errorSL = np.zeros(nDx)
    count = np.zeros(nDx)

    for i in range(nDx):
        nx = (i+1)*10
        nt = nx
        dx = (xmax - xmin)/nx
        x = np.arange(xmin, xmax, dx)
        phiOld = sineWave(x)
        phiAnalytic = sineWave((x - c*nt*dx)%(xmax - xmin))
        # implement the numerical scheme
        phiFTBS, TV_FTBS, numerical_mean_FTBS  = FTBS(phiOld.copy(), c, nt)
        phiCTCS, TV_CTCS, numerical_mean_CTCS = CTCS(phiOld.copy(), c, nt, False)
        phiBTCS, numerical_mean_BTCS = BTCS(phiOld.copy(), c, nt)
        phiLAGR, TV_SL, numerical_mean_SL = SemiLagrangian(phiOld.copy(), c, nt, dx, x)
        # calculate the error
        errorFTBS[i] = l2ErrorNorm(phiFTBS, phiAnalytic)
        errorBTCS[i] = l2ErrorNorm(phiBTCS, phiAnalytic)
        errorCTCS[i] = l2ErrorNorm(phiCTCS, phiAnalytic)
        errorSL[i] = l2ErrorNorm(phiLAGR, phiAnalytic)

        count[i] = dx

    # Plot the solutions
    font = {'size'   : 20}
    plt.rc('font', **font)

    plt.figure(1, figsize=(15,10))
    plt.clf()
    plt.ion()
    plt.loglog(count, errorFTBS, label='FTBS', color='black')
    plt.loglog(count, errorBTCS, label='BTCS', color='red')
    plt.loglog(count, errorCTCS, label='CTCS', color='blue')
    plt.loglog(count, errorSL, label='SemiLagrangian', color='yellow')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-1,10])
    plt.xlim(0.01,0.1)
    plt.legend()
    plt.xlabel('$dx$')
    plt.ylabel('$l_2$ error')
    plt.savefig('plots/accuracy.pdf')

accuracy()
