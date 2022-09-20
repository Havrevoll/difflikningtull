import numpy as np
from annan import annan_funksjon
from f import f
from scipy.integrate import solve_ivp
from particle import Particle

particle = Particle(diameter=0.05, init_position=[-3,2])

resultat = solve_ivp(f,(0,10),[particle.init_position[0],particle.init_position[1],0,0], method='BDF', args=(particle,),vectorized=True,t_eval=np.linspace(0.1,10,100),dense_output=True)

annan_funksjon(1)