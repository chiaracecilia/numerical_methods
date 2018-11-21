
import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from advectionSchemes import *
from diagnostics import *

# Parameters
def mean(initCond):

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

    exact_mean = np.sum(phiOld)

    phiFTBS, TV_FTBS, numerical_mean_FTBS  = FTBS(phiOld.copy(), c, nt)
    phiCTCS, TV_CTCS, numerical_mean_CTCS  = CTCS(phiOld.copy(), c, nt, False)
    phiCTCS_clipping, TV_CTCS_clipping, numerical_mean_CTCS_clipping  = CTCS(phiOld.copy(), c, nt, True)
    phiBTCS, numerical_mean_BTCS  = BTCS(phiOld.copy(), c, nt)
    phiLS, TV_SL, numerical_mean_SL  = SemiLagrangian(phiOld.copy(), c, nt, dx, x)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(t, numerical_mean_FTBS, ".-",label='FTBS', color='aqua')
    plt.plot(t, numerical_mean_CTCS, ".-",label='CTCS', color='red')
    plt.plot(t, numerical_mean_CTCS_clipping, ".-",label='CTCS with clipping', color='magenta')
    plt.plot(t, numerical_mean_BTCS, ".-",label='BTCS', color='blue')
    plt.plot(t, numerical_mean_SL, ".-",label='SemiLagrangian', color='navy')

    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([10,12])
    plt.legend()
    plt.xlabel('$x$')
    plt.savefig('plots/conservation_of_mass.pdf')

    plt.figure(2)
    plt.clf()
    plt.ion()
    plt.plot(t, numerical_mean_CTCS, ".-",label='CTCS', color='red')
    plt.plot(t, numerical_mean_CTCS_clipping, ".-",label='CTCS with clipping', color='magenta')

    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([10,12])
    plt.legend()
    plt.xlabel('$x$')
    plt.savefig('plots/conservation_of_mass_clipping.pdf')



    '''phi_CTCS, TV_CTCS  = CTCS(phiOld.copy(), c, nt)
    numerical_mean_CTCS = np.sum(phi_CTCS)

    phi_BTCS  = BTCS(phiOld.copy(), c, nt)
    numerical_mean_BTCS = np.sum(phi_BTCS)

    phiSL, TV_SL  = SemiLagrangian(phiOld.copy(), c, nt, dx, x)
    numerical_mean_SL = np.sum(phiSL)'''

    '''print("exact mean", exact_mean,"numerical mean FTBS", numerical_mean_FTBS, "numerical mean CTCS", numerical_mean_CTCS, "numerical mean BTCS", numerical_mean_BTCS, "numerical mean SL", numerical_mean_SL)'''

mean(squareWave)
