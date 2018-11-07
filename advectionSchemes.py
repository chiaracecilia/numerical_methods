# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py

import numpy as np

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

        # update arrays for next time-step. We do not care about all the time steps in between.sdjkl;
        phiOld = phi.copy()

    return phi

def FTBS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTBS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()

    # FTBS for each time-step
    for it in range(nt):
        # Loop through all space using remainder after division (%)
        # to cope with periodic boundary conditions
        for j in range(nx):
            phi[j] = phiOld[j] - c*\
                     (phiOld[(j)%nx] - phiOld[(j-1)%nx])
        # update arrays for next time-step
        phiOld = phi.copy()

    return phi

def CTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using CTCS, Courant number c"
    "for nt time-steps"

    nx = len(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()

    # in the CTCS to get phi at time t we need to use phi at time t-1
    # in the scheme we will save the previous step in another vector, but to deal with it in the initial condition we chose to apply FTCS scheme and obtain a second initial condition.
    # FTCS for the first time step to define a second initial condition
    for j in range(nx):
        phi[j] = phiOld[j] - c*\
            (phiOld[(j)%nx] - phiOld[(j-1)%nx])
    phi2 = phi.copy()


        # CTCS for the other time-steps. We will use phi2 to save the value of phi at time t-1 while we are calculating t+1
    for it in range(nt-1):
            # Loop through all space using remainder after division (%)
            # to cope with periodic boundary conditions
        for j in range(nx):
            phi[j] = phi2[j] - c*\
                     (phiOld[(j+1)%nx] - phiOld[(j-1)%nx])
            # update arrays for next time-step
        phi2 = phiOld.copy()
        phiOld = phi.copy()

    return phi

def BTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using BTCS, Courant number c"
    "for nt time-steps"
    nx = len(phiOld)
    #create the matrix M
    M = np.zeros((nx,nx), dtype=np.float)
    for j in range (nx):
        M[j][j]=1
        M[(j+1)%(nx)][j]=-c/2
        M[(j-1)%(nx)][j]=c/2

    # new time-step array for phi
    phi = phiOld.copy()

    for j in range(nt):
        phi = nplinalg.solve(M, phi)
        phiOld = phi.copy()

    return phi
