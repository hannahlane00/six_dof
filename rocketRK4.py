## RocketRK4, Harrison Leece
# t = time, ###Not currently used, important for future expansion###
# v = velocity
# m = mass,
# mFinal = Final mass of rocket, which is important for turning off
#          thrust and mass loss
# Cd = Coefficient of Drag, which is assumed constant, but is really
#      a function of height, but is assumed constant for simplicity
# thrust
# m2 = mass of fuel at 1/2 the step (for runge kutta 4th order k2 and k3 terms)
# h = stepSize
# Diff Eq form:  a = (Thrust/mass) - 32.2ft/s^2 (or gravity) - (.5*Area*rho)*(Cd*v^2)/mass
#                dv/dt = F(t,v); v is a fxn of time, mass is a fxn of t
#The RK4 method constitutes a "time step simulation"


def RK4(time, vel, mass, dry_mass, Cd, thrust, mDot, gravity, A, step_size, density):

####mass2 for k2 and k3, in effect is time + (step_size / 2)
    mass2 = mass - (mDot * step_size / 2)
####mass3 is for k4, in effect is time + step_size
    mass3 = mass - mDot * step_size

    #Thrust is controlled by thrustCurve.py and will turn off when thrust
    #is approx half of initial
    # Try catch trades tiny speed increases for debugging efficiency
    try:
        k1 = step_size * ((thrust / mass) - gravity - (.5 * A * density * Cd / mass * (vel)**2))
        vel2 = vel + (k1 / 2)
        k2 = step_size * ((thrust / mass2) - gravity - (.5 * A * density * Cd / mass2 * (vel2)**2))
        vel3 = vel + (k2 / 2)
        k3 = step_size * ((thrust / mass2) - gravity - (.5 * A * density * Cd / mass2 * (vel3)**2))
        vel4 = vel + (k3)
        k4 = step_size * ((thrust / mass3) - gravity - (.5 * A * density * Cd / mass3 * (vel4)**2))
    except:
        print("\nlikely overflow:  current system properties print")
        print("Thrust: " + str(thrust))
        print("Mass: " + str(mass))
        print("DryMass: " + str(dry_mass))
        print("Velocity: " + str(vel))
        print("Atmospheric Density: " + str(density))
        print("Step Size sanity check: " + str(step_size))
        print("Time sanity check: " + str(time))

    newVelocity = vel + ((k1 + 2 * k2 + 2 * k3 + k4) / 6)

    #print('Velocity: ' + str(newVelocity))
    return newVelocity
