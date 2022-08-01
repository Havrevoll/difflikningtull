import numpy as np

from constants import ν, ρ, g

def f(t,x, particle):
    '''inn: x = x,y,u,v
    ut: x = u,v, dudt,dvdt'''
    # v.x = -p.y/(1.+pow(p.x,2.))/(1.+pow(p.y,2.));
    # v.y = p.x/(1.+pow(p.y,2.))/(1.+pow(p.x,2.))-.3;

    number_of_vectors = x.shape[1]
    
    dxdt = x[2:]

    U_f = np.flip(x[:2]*[1,-1])/(1+np.square(x[:2]))/(1+np.square(np.flip(x[:2])))
    
    vel = U_f - dxdt # relativ snøggleik
    # vel_ang = atan2(vel[1], vel[0])
    assert len(vel) == 2
    Re = np.hypot(vel[0],vel[1]) * particle.diameter / ν
    
    # if (Re<1000):
    try:
            # with np.errstate(divide='raise'):
                # cd = 24 / Re * (1+0.15*Re**0.687)
            cd = ( (32 / Re)**(1/1.5) + 1)**1.5
            # Cheng (1997) skildrar Cd for kantete og runde steinar. Dette 
            # er kanskje den viktigaste grunnen til at eg bør gjera dette?
            # Ferguson og Church (2004) gjev nokre liknande bidrag, men viser til Cheng.
    except ZeroDivisionError:
            cd = 2e4
    # else:
    #     cd = 0.44
    
    # print("Re = ", Re," cd= ", cd)
    rho_self_density = ρ / particle.density
    
    drag_component =  3/4 * cd / particle.diameter * rho_self_density * abs(vel)*vel
    gravity_component = (rho_self_density - 1) * g

    # added_mass_component = 0.5 * rho_self_density * dudt_material 

    # pressure_component = rho_self_density * dudt_material
    
    # if np.all(dudt_material == 0.0 ):
    #     addedmass = False

    # lift_component = np.array([[0, -1],[1, 0]]) @ ( 3/4 * 0.5 / particle.diameter * rho_self_density * np.diff(np.square(U_top_bottom), axis=1).reshape(2,number_of_vectors) * norm(drag_component) )
    
    # divisor = 1 + 0.5 * rho_self_density * addedmass
    # divisoren trengst for akselerasjonen av partikkel kjem fram i added 
    # mass på høgre sida, så den må bli flytta over til venstre og delt på resten.
    
    dudt = drag_component + gravity_component 

    # print(f"{t};{x};{drag_component};{gravity_component};{added_mass_component - 0.5 * rho_self_density * dudt};{lift_component};{pressure_component};{dudt}",)

    return np.concatenate((dxdt,dudt))

    return x