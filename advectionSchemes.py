# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py

import numpy as np
from scipy.interpolate import lagrange

def FTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTCS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)
    # new time-step array for phi
    phi = phiOld.copy()

    # FTCS for each time-step
    for it in range(nt):
        # Loop through all space using remainder after division (%). (a%d=a if d is bigger than a)
        # to cope with periodic boundary conditions
        for j in range(nx):
            phi[j] = phiOld[j] - 0.5*c*\
                     (phiOld[(j+1)%nx] - phiOld[(j-1)%nx])

        # update arrays for next time-step. We do not care about all the time steps in between.
        phiOld = phi.copy()

    return phi

def FTBS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTBS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)
    TV = np.zeros(nt+1)
    numerical_mean_FTBS = np.zeros(nt+1)
    for j in range(nx):
        TV[0] = TV[0] + abs(phiOld[j] - phiOld[(j-1)%nx])
    numerical_mean_FTBS[0] = np.sum(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()

    # FTBS for each time-step
    for it in range(nt):
        T = 0
        # Loop through all space using remainder after division (%)
        # to cope with periodic boundary conditions
        for j in range(nx):
            phi[j] = phiOld[j] - c*\
                     (phiOld[(j)%nx] - phiOld[(j-1)%nx])
            T = T + abs(phi[j] - phi[(j-1)%nx])
            # update arrays for next time-step
        numerical_mean_FTBS[it+1]= np.sum(phi)
        TV[it+1] = T
            # update arrays for next time-step
        phiOld = phi.copy()

    return phi, TV, numerical_mean_FTBS

def CTCS(phiOld, c, nt, clipping_enabled):
    "Linear advection of profile in phiOld using CTCS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()

    TV = np.zeros(nt+1)
    numerical_mean_CTCS = np.zeros(nt+1)

    # in the CTCS to get phi at time t we need to use phi at time t-1
    # in the scheme we will save the previous step in another vector,
    # but to deal with it in the initial condition we chose to apply FTCS scheme and obtain a second initial condition.
    # FTCS for the first time step to define a second initial condition
    for j in range(nx):
        phi[j] = phiOld[j] - c*\
            (phiOld[(j)%nx] - phiOld[(j-1)%nx])
        TV[0] = TV[0] + abs(phiOld[j] - phiOld[(j-1)%nx])
    numerical_mean_CTCS[0] = np.sum(phiOld)
    phi2 = phi.copy()


    # CTCS for the other time-steps. We will use phi2 to save the value of phi at time t-1 while we are calculating t+1
    for it in range(nt):
            # Loop through all space using remainder after division (%)
            # to cope with periodic boundary conditions
        T = 0
        for j in range(nx):
            phi[j] = phi2[j] - c*\
                     (phiOld[(j+1)%nx] - phiOld[(j-1)%nx])
            if clipping_enabled and phi[j] > 1:
                phi[j] = 1
            elif clipping_enabled and phi[j] < 0:
                phi[j] = 0
            T = T + abs(phi[j] - phi[(j-1)%nx])
        #update arrays for the next time step
        TV[it+1] = T
        phi2 = phiOld.copy()
        phiOld = phi.copy()

        numerical_mean_CTCS[it+1] = np.sum(phi)
    return phi, TV, numerical_mean_CTCS

def BTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using BTCS, Courant number c"
    "for nt time-steps"
    nx = len(phiOld)
    numerical_mean_BTCS = np.zeros(nt+1)
    numerical_mean_BTCS[0] = np.sum(phiOld)
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
        numerical_mean_BTCS[j+1] = np.sum(phi)
    return phi, numerical_mean_BTCS

def SemiLagrangian(phiOld, c, nt, dx, x):
    "linear advection of profile in phiOld using SemiLagrangian"
    "Courant number c, for nt time-steps"
    nx = len(phiOld)
    TV = np.zeros(nt+1)
    numerical_mean_SL = np.zeros(nt+1)
    numerical_mean_SL[0] = np.sum(phiOld)
    for j in range(nx):
        TV[0] = TV[0] + abs(phiOld[j] - phiOld[(j-1)%nx])

    # new time-step array for phi
    phi = phiOld.copy()

    for it in range(nt):
        T = 0
        for j in range(nx):

            k = int(np.floor(j - c)) #k s.t. x_j is in (x_k,x_(k+1))
            beta = j - c - k   #coefficient needed for the interpolation
            x_jd = beta*dx + x[k] #point that needs to be interpolated
            # set up a cubic-Lagrange interpolation of four points
            x_int = np.array([x[k-1],x[k],x[(k+1)%nx],x[(k+2)%nx]]) #we interpolate using four grid points around the point of interest
            phi_int = np.array([phiOld[k-1],phiOld[k],phiOld[(k+1)%nx],phiOld[(k+2)%nx]]) #corrispondent y-point
            inter_polinomial = lagrange(x_int, phi_int) #calculate the interpolating polynomial
            phi[j] = inter_polinomial(x_jd) #interpolated value
            T = T + abs(phi[j] - phi[(j-1)%nx])
        TV[it+1] = T
        numerical_mean_SL[it+1] = np.sum(phi)
        phiOld = phi.copy() # Update for next time step

    return phi, TV, numerical_mean_SL
