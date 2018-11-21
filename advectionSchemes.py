# Numerical schemes for simulating linear advection for outer codes
# linearAdvect.py, courant_number_dependence.py, mass_conservation'py, graphs.py

import numpy as np
from scipy.interpolate import lagrange

def FTBS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTBS, Courant number c, "
    "for nt time-steps. It also gives as output the numerical mass at every time"
    "and the total variation of the advected quantity"

    nx = len(phiOld)
    # new array for the total variation at each time step
    TV = np.zeros(nt+1)
    # new array for the mass at each time step
    numericalMassFTBS = np.zeros(nt+1)

    # calculate the total variation of the initial initialConditions
    for j in range(nx):
        TV[0] = TV[0] + abs(phiOld[j] - phiOld[(j-1)%nx])
    # calculate the mass for the initial condition.
    numericalMassFTBS[0] = np.sum(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()

    # FTBS for each time-step
    for it in range(nt):
        # T will be used to calculate the total variation
        T = 0
        # Loop through all space using remainder after division (%)
        # to cope with periodic boundary conditions
        for j in range(nx):
            phi[j] = phiOld[j] - c*\
                     (phiOld[(j)%nx] - phiOld[(j-1)%nx])
            # calculate the partial variation up to the jth space step
            T = T + abs(phi[j] - phi[(j-1)%nx])
            # update arrays for next time-step
        numericalMassFTBS[it+1]= np.sum(phi)
        # store the total variation of the advected quantity
        TV[it+1] = T
            # update arrays for next time-step
        phiOld = phi.copy()

    return phi, TV, numericalMassFTBS

def CTCS(phiOld, c, nt, clippingEnabled):
    "Linear advection of profile in phiOld using CTCS, Courant number c"
    "for nt time-steps, with the possibility of considering clipping"
    "It also gives as output the numerical mass at every time,"

    nx = len(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()
    # new array for the total variation at each time step
    TV = np.zeros(nt+1)
    # new array for the mass at each time step
    numericalMassCTCS = np.zeros(nt+1)

    # in the CTCS to get phi at time t we need to use phi at time t-1
    # in the scheme we will save the previous step in another vector,
    # but to deal with it in the initial condition we chose to apply FTCS
    #scheme and obtain a second initial condition.

    # FTCS for the first time step to define a second initial condition
    for j in range(nx):
        phi[j] = phiOld[j] - c*\
            (phiOld[(j)%nx] - phiOld[(j-1)%nx])
        #total variation at intial time
        TV[0] = TV[0] + abs(phiOld[j] - phiOld[(j-1)%nx])
    # mass at initial time
    numericalMassCTCS[0] = np.sum(phiOld)
    phi2 = phi.copy()

    # CTCS for the other time-steps. We will use phi2 to save the value of phi at time t-1 while we are calculating t+1
    for it in range(nt):
            # Loop through all space using remainder after division (%)
            # to cope with periodic boundary conditions
        # T is used to compute total variation
        T = 0
        for j in range(nx):
            phi[j] = phi2[j] - c*\
                     (phiOld[(j+1)%nx] - phiOld[(j-1)%nx])
            # now check: if the values obtained are out of the bound we choose,
            # we set them to be automatically equal to the value at the boundary
            if clippingEnabled and phi[j] > 1:
                phi[j] = 1
            elif clippingEnabled and phi[j] < 0:
                phi[j] = 0
            T = T + abs(phi[j] - phi[(j-1)%nx])
        # total variation
        TV[it+1] = T
        #update arrays for the next time step
        phi2 = phiOld.copy()
        phiOld = phi.copy()
        # numerical mass is evaluated
        numericalMassCTCS[it+1] = np.sum(phi)

    return phi, TV, numericalMassCTCS

def BTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using BTCS, Courant number c"
    "for nt time-steps"
    "It also gives as output the numerical mass at every time"

    nx = len(phiOld)
    # set up a vector for the numerical mass
    numericalMassBTCS = np.zeros(nt+1)
    # calculate the mass at initial time
    numericalMassBTCS[0] = np.sum(phiOld)

    #create the matrix M
    M = np.zeros((nx,nx), dtype=np.float)
    for j in range (nx):
        M[j][j]=1
        M[(j+1)%(nx)][j]=-c/2
        M[(j-1)%(nx)][j]=c/2

    # new time-step array for phi
    phi = phiOld.copy()

    for j in range(nt):
        phi = np.linalg.solve(M, phi)
        # calculate the numerical mass
        numericalMassBTCS[j+1] = np.sum(phi)

    return phi, numericalMassBTCS

def SemiLagrangian(phiOld, c, nt, dx, x):
    "linear advection of profile in phiOld using SemiLagrangian"
    "Courant number c, for nt time-steps"
    "It also gives as output the numerical mass at every time,"
    "the l2 and linf error norm and the total variation of the advected quantity"

    nx = len(phiOld)
    # set up a vector for the total variation
    TV = np.zeros(nt+1)
    # set up a vector for the numerical mean
    numericalMassSL = np.zeros(nt+1)
    # calculate the mass for the initial condition
    numericalMassSL[0] = np.sum(phiOld)
    # calculate the total variation for th initial condition
    for j in range(nx):
        TV[0] = TV[0] + abs(phiOld[j] - phiOld[(j-1)%nx])

    # new time-step array for phi
    phi = phiOld.copy()

    for it in range(nt):
        # T is used to store the partial variation
        T = 0
        for j in range(nx):
            #k s.t. x_j is in (x_k,x_(k+1))
            k = int(np.floor(j - c))
            #coefficient needed for the interpolation
            beta = j - c - k
            #point that needs to be interpolated
            x_jd = beta*dx + x[k]
            # set up a cubic-Lagrange interpolation of four points
            x_int = np.array([x[k-1],x[k],x[(k+1)%nx],x[(k+2)%nx]])
            # we interpolate using four grid points around the point of interest
            # corrispondent y-point
            phi_int = np.array([phiOld[k-1],phiOld[k],phiOld[(k+1)%nx],phiOld[(k+2)%nx]])
            # calculate the interpolating polynomial
            inter_polinomial = lagrange(x_int, phi_int)
            # interpolated value
            phi[j] = inter_polinomial(x_jd)
            # partial variation
            T = T + abs(phi[j] - phi[(j-1)%nx])
        # total variation
        TV[it+1] = T
        # calculate the mass
        numericalMassSL[it+1] = np.sum(phi)
        # Update for next time step
        phiOld = phi.copy()

    return phi, TV, numericalMassSL
