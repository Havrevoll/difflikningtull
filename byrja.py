import numpy as np
from annan import annan_funksjon
from f import f
from scipy.integrate import solve_ivp
from particle import Particle

particle = Particle(diameter=0.05, init_position=[-3,2])

solve_ivp(f,(0,10),0, method='BDF', args=(particle))