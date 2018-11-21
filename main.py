# this code implements all the functions written in the other scripts.
# it is possible to call them with different parameters than the ones specified here

import matplotlib.pyplot as plt
import timeit

start = timeit.default_timer()
# read in all the linear advection schemes, initial conditions and other
# code associated with this applications
from initialConditions import *
from advectionSchemes import *
from diagnostics import *
from graphs import *
from runTime import *
from courant_number_dependence import *
from massConservation import *
from orderAccuracy import *


# run these functions
print('Running main process...')
schemeComparison(cosBell, 'cosineBell')
print('Scheme comparison  with cosine bell initial condition exectuted')
schemeComparison(squareWave, 'squareWave')
print('Scheme comparison  with square wave initial condition exectuted')
courant_number_dependence(cosBell, 3)
print('courant number analysis executed')
mass(squareWave)
print('mass conservation executed')
accuracy()
print('order of accuracy executed')
runTime()
end = timeit.default_timer()

time = end - start
print('...main succesfully executed in: ', time, 'sec')
