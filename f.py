import numpy as np

from constants import ν, ρ, g

def f(t,x, particle):
    '''inn: x = x,y,u,v
    ut: x = u,v, dudt,dvdt'''
    # v.x = -p.y/(1.+pow(p.x,2.))/(1.+pow(p.y,2.));
    # v.y = p.x/(1.+pow(p.y,2.))/(1.+pow(p.x,2.))-.3;
    # https://anvaka.github.io/fieldplay/?cx=0.004050000000000331&cy=0.014800000000000146&w=20.0081&h=20.0081&dt=0.01&fo=0.998&dp=0.001&cm=2&vf=%2F%2F%20p.x%20and%20p.y%20are%20current%20coordinates%0A%2F%2F%20v.x%20and%20v.y%20is%20a%20velocity%20at%20point%20p%0Avec2%20get_velocity%28vec2%20p%29%20%7B%0A%20%20vec2%20v%20%3D%20vec2%280.%2C%200.%29%3B%0A%20%20float%20t%20%3D%205.%3B%20%0A%20%20%2F%2F%20change%20this%20to%20get%20a%20new%20vector%20field%0A%20%20v.x%20%3D%20%20p.y%2F%28.5%20%2B%20t*pow%28p.x%2C2.%29%2Bt*pow%28p.y%2C2.%29%29%3B%0A%20%20v.y%20%3D%20-p.x%2F%28.5%20%2B%20t*pow%28p.y%2C2.%29%2Bt*pow%28p.x%2C2.%29%29%3B%0A%0A%20%20return%20v%3B%0A%7D&code=%2F%2F%20p.x%20and%20p.y%20are%20current%20coordinates%0A%2F%2F%20v.x%20and%20v.y%20is%20a%20velocity%20at%20point%20p%0Avec2%20get_velocity%28vec2%20p%29%20%7B%0A%20%20vec2%20v%20%3D%20vec2%280.%2C%200.%29%3B%0A%20%20float%20t%20%3D%205.%3B%20%0A%20%20%2F%2F%20change%20this%20to%20get%20a%20new%20vector%20field%0A%20%20v.x%20%3D%20%20p.y%2F%28.5%20%2B%20t*pow%28p.x%2C2.%29%2Bt*pow%28p.y%2C2.%29%29%3B%0A%20%20v.y%20%3D%20-p.x%2F%28.5%20%2B%20t*pow%28p.y%2C2.%29%2Bt*pow%28p.x%2C2.%29%29%3B%0A%0A%20%20return%20v%3B%0A%7D&pc=5000
    

    number_of_vectors = x.shape[1]
    # g = np.asarray([[0.],[-g]])
    dxdt = x[2:]

    U_f = np.flip(x[:2]*np.asarray([[1],[-1]]))/(1+np.square(x[:2])+np.square(np.flip(x[:2])))
    
    vel = U_f - dxdt # relativ snøggleik
    # vel_ang = atan2(vel[1], vel[0])
    assert len(vel) == 2
    Re = np.hypot(vel[0],vel[1]) * particle.diameter / ν
    
    # if (Re<1000):
    try:
            # with np.errstate(divide='raise'):
                # cd = 24 / Re * (1+0.15*Re**0.687)
            cd = ( (32 / (Re+0.001))**(1/1.5) + 1)**1.5
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
    gravity_component = (rho_self_density - 1) * np.asarray([[0.],[g]])

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