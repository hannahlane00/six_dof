"""
Rocket.py
1/20/20
"""

import math

import numpy as np
from scipy.integrate import odeint

from rocketRK4 import RK4
from SIXDOF import SIXDOF
import matplotlib.pyplot as plt
from atmosphereDensity import seMagic as compDensity
from thrustCurve import newThrustMdot as compThrust
from thrustCurve import pressurantThrust as compThrust2
from waveDrag import compTemp
from waveDrag import compMach
from waveDrag import waveDrag

#compDensity from atmosphereDensity
def seMagic(alt):
    s = -.2457039080 * (alt / 10**4)**1 - .0351850842 * (alt / 10**4)**2 + 0.0021044026 * (alt / 10**4)**3 - 0.0000390562 * (alt / 10**4)**4
    return 23.77 * 10**(-4) * np.exp(s)

#compThrust and compThrust2 from thrustCurve
class thrustCurve():

    def __init__(self,initialThrust, initialMdot, inMass, curMass):
        self.initialThrust = initialThrust
        self.initialMdot = initialMdot
        self.inMass = inMass
        self.curMass = curMass

    def newThrustMdot(self):
        #a linear relationship between tank pressure and mDot is assumed
        #Then a linear elationship between mDot is assumed
        #Mass of helium pressurant is considered neglegible
        #Therefore, a lienar relationship between current propellant mass and thrust
        #and mass flow rate can be established
    
        ratio = self.curMass / self.inMass
        nThrust = self.initialThrust * ratio
        #Remember: initialMdot is computed from thrust/(g*isp) outside of main while loop
        nMdot = self.initialMdot * ratio
    
        #print(initialThrust)
        #print(nThrust)
        #print(nMdot)
    
        #Very important piece of logic; I assume thrust will cut out before tanks are dry
        #Because gasses will not attain sonic velocity at throat
        #Turns off thrust, which will then control RK4 and burnout time logic.
        if nMdot < (.1 * self.initialMdot):
            nThrust = 0
        if nMdot < .0002:
            nMdot = 0
        return nThrust, nMdot

    def pressurantThrust(self):
        thrust = self.initialThrust
        mDot = self.initialMdot
        if self.curMass < .1:
            thrust = 0
            mDot = 0
        return thrust, mDot

#compTemp, compMach, waveDrag from waveDrag
class WaveDrag():

    def __init__(self, ftv, T, ftx, cd, mach):
        self.ftv = ftv
        self.T = T
        self.ftx = ftx
        self.cd = cd
        self.mach = mach

    def compMach(self):
        #convert v from ft/s to m/s
        v = self.ftv * .3048
        #R is kJ/(kg*K)
        """R??, V??, A??"""
        R = .287
        #Specific heat ratio (gamma) is approx 1.4 for whole domain of temperatures
        gamma = 1.4
        #speed of sound = a, meters/s
        a = (gamma * R * 1000 * self.T)**.5
        self.mach = v / a
        return self.mach

    def compTemp(self):
        x = self.ftx * .3048
        """X??,C,T"""
        if(x > 0 and x < 15354.016):
            C = 2.5574 * 10**(-7) * x**2 - 8.8754 * 10**(-3) * x + 18.1061
        elif(x > 15354.016 and x < 49697.695):
            C = 3.6838 * 10**(-8) * x**2 -7.6707 * 10**(-4) * x - 54.7841
        elif(x > 49697.695 and x < 120000):
            C = -2.4347 * 10**(-3) * x + 119.078
        else:
            C = -172
        self.T = C + 273.15
        return self.T
    
    def waveDrag(self):
        """CD??,ADJ??"""
        adj = self.cd - .3
        if(self.mach < .5085):
            nCd = .6827 * self.mach**3 - .4297 * self.mach**2 - .0358 * self.mach + .3 + adj
        elif(self.mach > .5085 and self.mach < 1.3618):
            nCd = -0.7809 * self.mach**4 + 2.324 * self.mach**3 - 2.0189 * self.mach**2 + 0.6793 * self.mach + 0.1837 + adj
        elif(self.mach > 1.3618 and self.mach < 4):
            nCd = -.003495 * self.mach**6 + .07004 * self.mach**5 -.583 * self.mach**4 + 2.564 * self.mach**3 -6.186 * self.mach**2 + 7.466 * self.mach -2.923 + adj
        else:
            nCd = .2184 + adj
        return nCd


#F_wr is wind in r direction
#import drag, gravity, time, mass, thrust, gamma, alpha, and wind in all directions
class SIXDOF:
    def __init__(self, x_1, x_2, x_3, x_4, x_5, x_6, fw, alt, x, step_size, nthrust):
        self.state = [x_1, x_2, x_3, x_4, x_5, x_6]
        self.drag = 0
        self.gravity = 32.174 #ft / s^2
        self.time = 0
        self.thrust = 0
        self.gamma = 88
        self.alphap = 2
        self.stateZ = [0, 0, 0, 0, 0, 0]
        self.fw = fw
        self.alt = alt
        self.x = x
        self.step_size = step_size
        self.nthrust = nthrust
        
    def sixeqdiff(self):
        x_1 = self.state[0]
        x_2 = self.state[1]
        x_3 = self.state[2]
        x_4 = self.state[3]
        x_5 = self.state[4]
        x_6 = self.state[5]

        F_wr = self.fw
        F_wtheta = 0
        F_wphi = 0

        x_dot_1 = x_2
        x_dot_2 = -self.gravity + (-self.drag * sin(self.gamma) + (self.thrust * sin(self.gamma + self.alphap)) + F_wr) / mass + x_1 * x_4**2 * sin(x_5)**2 + x_1 * x_5**2
        x_dot_3 = x_4
        x_dot_4 = ((-self.drag * cos(self.gamma) + self.thrust * cos(self.gamma + self.alphap) + F_wtheta) / m - 2 * x_2 * x_4 * sin(x_5) - 2 * x_1 * x_2 * x_6 * cos(x_5)) * (1 / (x_1 * sin(x_5)))
        x_dot_5 = x_6
        x_dot_6 = ((F_wphi / mass) - 2 * x_2 * x_6 + x_1 * x_4**2 * cos(x_5) * sin(x_5)) * (1 / x_1)
        return  x_dot_1, x_dot_2, x_dot_3, x_dot_4, x_dot_5, x_dot_6
    
    def windForce(self):
        #determines wind velocity value for certain altitude (alt)
        altitude = self.alt
        if altitude <= 1 and altitude >= 0:
            windVel = 5 * (altitude - 1) + 8
        elif altitude <= 2 and altitude > 1:
            windVel = 8
        elif altitude <= 3 and altitude > 2:
            windVel = -6 * altitude + 20
        elif  altitude <= 11 and altitude > 3:
            windVel = ((43 * altitude) - 113.004) / 8
        elif altitude <= 15 and altitude > 11:
            windVel = 45
        elif altitude < 22 and altitude > 15:
            windVel = -6 * altitude + 135
        elif altitude <= 25 and altitude >= 22:
            windVel = 4
        elif altitude <= 30 and altitude > 25:
            windVel = 2 * altitude - 46
        else:
            windVel = 0

        
        #rho = ad.seMagic(alt)
        #rho = 5
        #S = 22 #ft

        return self.alt, windVel


    def gusts(self):
        num_gusts = random.randrange(0, 20)
        alt_range = random.randrange(0, 30000, 1)
        mag_gust = random.randrange(0, 12, 1)
        direction = random.choice((-1, 1))

        #creates an array of values and an empty list for all gust altitudes
        gust_matrix = np.array([0, alt_range, mag_gust, direction])
        gust_altitudes = []
        
        for i in range(num_gusts-1):
            alt_range = random.randrange(0, 30000, 1)
            mag_gust = random.randrange(0, 12, 1)
            direction = random.choice((-1, 1))
            temp_array = np.array([i+1, alt_range, mag_gust, direction])
            gust_matrix = np.vstack((temp_array, gust_matrix))
            gust_altitudes.append(alt_range / 1000)
        gust_matrix = np.flipud(gust_matrix)
        
        return gust_matrix

    def sixdofsolver(self):
        init_state = [x_1, x_2, x_3, x_4, x_5, x_6]
        sixeq = sixeqdiff(self, fw)
        state = odeint(sixeq, init_state)
        return state


def RK4(time, vel, mass, dry_mass, Cd, thrust, mDot, gravity, A, step_size, density, x_3, x_5, x_6):

####mass2 for k2 and k3, in effect is time + (step_size / 2)
    mass2 = mass - (mDot * step_size / 2)
####mass3 is for k4, in effect is time + step_size
    mass3 = mass - mDot * step_size

    Vx = abs(vel) * cos(theta)
    Vy = abs(vel) * sin(b)
    Vz = abs(vel) * sin(theta)
    #Thrust is controlled by thrustCurve.py and will turn off when thrust
    #is approx half of initial
    # Try catch trades tiny speed increases for debugging efficiency

    #Vx
    try:
        k1 = step_size * ((thrust / mass) - gravity - (.5 * A * density * Cd / mass * (Vx)**2))
        velx2 = Vx + (k1 / 2)
        k2 = step_size * ((thrust / mass2) - gravity - (.5 * A * density * Cd / mass2 * (velx2)**2))
        velx3 = Vx + (k2 / 2)
        k3 = step_size * ((thrust / mass2) - gravity - (.5 * A * density * Cd / mass2 * (velx3)**2))
        velx4 = Vx + (k3)
        k4 = step_size * ((thrust / mass3) - gravity - (.5 * A * density * Cd / mass3 * (velx4)**2))
    except:
        print("\nlikely overflow:  current system properties print")
        print("Thrust: " + str(thrust))
        print("Mass: " + str(mass))
        print("DryMass: " + str(dry_mass))
        print("Velocity: " + str(vel))
        print("Atmospheric Density: " + str(density))
        print("Step Size sanity check: " + str(step_size))
        print("Time sanity check: " + str(time))
       
    newVelocity_x = velx + ((k1 + 2 * k2 + 2 * k3 + k4) / 6)

    #Vy
    try:
        k1 = step_size * ((thrust / mass) - gravity - (.5 * A * density * Cd / mass * (Vy)**2))
        vely2 = Vy + (k1 / 2)
        k2 = step_size * ((thrust / mass2) - gravity - (.5 * A * density * Cd / mass2 * (vely2)**2))
        vely3 = Vy + (k2 / 2)
        k3 = step_size * ((thrust / mass2) - gravity - (.5 * A * density * Cd / mass2 * (vely3)**2))
        vely4 = Vy + (k3)
        k4 = step_size * ((thrust / mass3) - gravity - (.5 * A * density * Cd / mass3 * (vely4)**2))
    except:
        print("\nlikely overflow:  current system properties print")
        print("Thrust: " + str(thrust))
        print("Mass: " + str(mass))
        print("DryMass: " + str(dry_mass))
        print("Velocity: " + str(vel))
        print("Atmospheric Density: " + str(density))
        print("Step Size sanity check: " + str(step_size))
        print("Time sanity check: " + str(time))
       
    newVelocity_y = vely + ((k1 + 2 * k2 + 2 * k3 + k4) / 6)

    #Vz
    try:
        k1 = step_size * ((thrust / mass) - gravity - (.5 * A * density * Cd / mass * (Vz)**2))
        velz2 = Vz + (k1 / 2)
        k2 = step_size * ((thrust / mass2) - gravity - (.5 * A * density * Cd / mass2 * (velz2)**2))
        velz3 = Vz + (k2 / 2)
        k3 = step_size * ((thrust / mass2) - gravity - (.5 * A * density * Cd / mass2 * (velz3)**2))
        velz4 = Vz + (k3)
        k4 = step_size * ((thrust / mass3) - gravity - (.5 * A * density * Cd / mass3 * (velz4)**2))
    except:
        print("\nlikely overflow:  current system properties print")
        print("Thrust: " + str(thrust))
        print("Mass: " + str(mass))
        print("DryMass: " + str(dry_mass))
        print("Velocity: " + str(vel))
        print("Atmospheric Density: " + str(density))
        print("Step Size sanity check: " + str(step_size))
        print("Time sanity check: " + str(time))
       
    newVelocity_z = velz + ((k1 + 2 * k2 + 2 * k3 + k4) / 6)
    newVelocity_r = math.sqrt(newVelocity_x**2 + newVelocity_y**2 + newVelocity_z**2)
    #print('Velocity: ' + str(newVelocity))
    return newVelocity_x, newVelocity_y, newVelocity_z, newVelocity_r

#Create a Rocket Class
class Rocket():

    def __init__(self, thrust_input, weight, pr_ratio, isp, diameter, drag_coefficient):
        '''    weight :            total rocket weight
             pr_ratio :          mass of propelants / total wet mass of rockets
             isp :               specific impulse; rating of rocket efficiency
             diameter:           rocket body
             drag_coefficient:   coefficient of drag
        '''
        self.thrust = thrust_input
        self.gravity = 32.174  # ft / s^2
        self.initial_mass = weight / self.gravity
        self.propellant_mass = pr_ratio * self.initial_mass
        self.final_mass = self.initial_mass - self.propellant_mass
        self.area = 3.1416 * (1/4) * diameter**2
        self.m_dot = thrust_input / (isp * self.gravity)
        self.cd = drag_coefficient
        self.pr_ratio = pr_ratio

    def flight_simulation(self):
        step_size = .01
        #initialize rocket data from __init__
        thrust = self.thrust
        gravity = self.gravity  # ft / s^2
        initial_mass = self.initial_mass
        propellant_mass = self.propellant_mass
        final_mass = self.final_mass
        area = self.area
        m_dot = self.m_dot
        cd = self.cd
        pr_ratio = self.pr_ratio
        mass = initial_mass
        #Initializes a bunch of stuff locally
        x_graph = np.array([0])
        v_graph = np.array([0])
        t_graph = np.array([0])
        #Generally, do not change these.  Describes initial conditions and
        #initializes the variables
        t = 0
        vx = 0
        vy = 0
        vz = 0
        x = 0
        rho = 0
        bo_time = 0
        bo_velocity = 0
        current_velocity = 0
        max_q = 0
        max_q_velocity = 0
        max_q_altitude = 0
        max_q_acceleration = 0
        [x_1, x_2, x_3, x_4, x_5, x_6] = [0, 0, 0, 0, 0, 0]
        gamma = 88
        alphap = 0
# initialize force of wind

        #While loop runs until rocket velocity is less than 10 feet/s
        while x >= 0:
            #Calculates atmospheric density at altitude x and passes density to RK4 fxn.
            rho = compDensity(x)
            #same x as displacement x?
            
            #Calculates mDot and thrust as tank pressure drops (if applicable) and  altitude changes
            n_thrust, n_mdot = compThrust2(thrust, m_dot,
                                           (initial_mass * pr_ratio),
                                           propellant_mass)
            #Calculates Cd as a fx of mach number
            K = compTemp(x)
            #again what is this x
            
            mach = compMach(v, K)
            # does this need to be calculated for each new velocity?
            
            n_cd = waveDrag(cd, mach)
            Drag = n_cd * (1/2) * rho * v**2 * area
            # does this need new velocity
            
            print('First print reached')
            sf = SIXDOF(x_1, x_2, x_3, x_4, x_5, x_6)
            [x_1, x_2, x_3, x_4, x_5, x_6] = sf.sixdofsolver(x_1, x_2, x_3, x_4,
                                                             x_5, x_6, Drag,
                                                             n_thrust, x, step_size,
                                                             gravity, gamma, alphap)
            print(x_1, x_2, x_3, x_4, x_5, x_6)
            #RK4 function to solve rocket differential equation.  Returns v(t) approximated using Runge-Kutta 4th order method
            new_velocity = RK4(t, v, mass, final_mass,
                               n_cd, n_thrust, n_mdot,
                               gravity, area, step_size, rho, x_3, x_5, x_6)
            """
            new_velocity = RK4(t, x_2, mass, final_mass,
                                n_cd, n_thrust, n_mdot,
                                gravity, area, step_size, rho)
            """
            #is final mass a constant? Where is it updating

            
            #Computes current displacement using reimen sum
            x = x + ((current_velocity + new_velocity) / 2 * step_size)
            #^^^?
            if n_thrust > 0:
                mass = mass - n_mdot * step_size
                propellant_mass = propellant_mass - n_mdot * step_size
                #Get burnout
                bo_time = t
                bo_velocity = v
                bo_thrust = n_thrust
                bo_mdot = n_mdot
                bo_alt = x
                bo_mass = mass
                # bo because theyre all 0 i assume.. except for mass... is it supposed to be?
                bo_accel = ((thrust / mass)
                            / 32.2
                            - ((area * n_cd * .5 * rho * v**2) / mass)
                            / 32.2 - 1)
                bo_mdot = n_mdot
                temp_max_q = n_cd * rho * 1/2 * v**2
                if max_q < temp_max_q:
                    max_q = temp_max_q
                    max_q_velocity = v
                    max_q_altitude = x
                    max_q_time = t
                    #using the time and step size to find slope of velocity graph
                    if np.size(v_graph) > 3:
                        max_q_acceleration = ((v_graph[int(np.floor(t/step_size))]
                                              - v_graph[int(np.floor(t/step_size)-1)])
                                              / step_size)

            #Thrust is off, but propellant is still being drained by the tank pressue fluid (if applicable, such as blowdown system)
            elif n_thrust == 0 and n_mdot != 0:
                mass = mass - n_mdot * step_size
                propellant_mass = propellant_mass - n_mdot * step_size
            else:
                #Leave this for clarity, though redundant
                mass = mass
            x_graph = np.append(x_graph, x)
            v_graph = np.append(v_graph, v)
            current_velocity = v
            #steps the while loop forward in time
            t = t + step_size
            t_graph = np.append(t_graph, t)


        #returns a tuple containing vectors containing displacement, velocity and time, maximum velocity, burnouttime, and max altitude
        return x_graph, v_graph, t_graph, bo_velocity, bo_time, np.max(x_graph), bo_alt, bo_mass, t_graph[np.argmax(x_graph)], bo_accel, bo_mdot, max_q, max_q_velocity, max_q_altitude, max_q_acceleration, max_q_time
