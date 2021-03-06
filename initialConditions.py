# Various different initial conditions for linear advection

import numpy as np

def squareWave(x,alpha,beta):
    "A square wave as a function of position, x, which is 1 between alpha"
    "and beta and zero elsewhere. The initialisation is conservative so"
    "that each phi contains the correct quantity integrated over a region"
    "a distance dx/2 either side of x"

    phi = np.zeros_like(x)

    # The grid spacing (assumed uniform)
    dx = x[1] - x[0]

    # Set phi away from the end points (assume zero at the end points)
    for j in range(1,len(x)-1):
        # edges of the grid box (using west and east notation)
        xw = x[j] - 0.5*dx
        xe = x[j] + 0.5*dx

        #integral quantity of phi
        phi[j] = max((min(beta, xe) - max(alpha, xw))/dx, 0)

    return phi


def cosBell(x, alpha=0, beta=0.5):
    "Function defining a cosine bell as a function of position, x"
    "between alpha and beta with default parameters 0, 0.5"
### The lambda keyword lets you define a function in one line       ###
    width = beta - alpha
    bell = lambda x: 0.5*(1 - np.cos(2*np.pi*(x-alpha)/width))
### chooses bell(x) where condition is true, else chooses zeros     ###
    return np.where((x<beta) & (x>=alpha), bell(x), 0.)


def sineWave(x):
    "A sine wave transformed to take values from 0 to 1"
    return 0.5+0.5*np.sin(2*np.pi*x)
