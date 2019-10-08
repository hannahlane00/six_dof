# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 14:23:41 2019

@author: jwhribal
"""
from math import sin
from atmosphereDensity import seMagic as ad
from scipy.integrate import odeint
import random 
# importing the required module
import matplotlib.pyplot as plt
import numpy as np

drag = 0
gravity = 32.174  # ft / s^2
time = 0
thrust = 0
gamma = 88
alphap = 2
state = [0,0,0,0,0,0]


#F_wr is wind in r direction
#import drag, gravity, time, mass, thrust, gamma, alpha, and wind in all directions
class SIXDOF:
    """def __init__(self, x_1, x_2, x_3, x_4, x_5, x_6):

        self.state = [x_1, x_2, x_3, x_4, x_5, x_6]


    def sixeqdiff(self, fw):
        x_1 = self.state[0]
        x_2 = self.state[1]
        x_3 = self.state[2]
        x_4 = self.state[3]
        x_5 = self.state[4]
        x_6 = self.state[5]

            
        F_wr = Fw
        F_wtheta = 0
        F_wphi = 0

        x_dot_1 = x_2
        x_dot_2 = -gravity + (-Drag * sin(gamma) + (thrust * sin(gamma + alphap)) + F_wr) / (mass) + x_1 * x_4**2 * sin(x_5)**2 + x_1 * x_5**2
        x_dot_3 = x_4
        x_dot_4 = ((-Drag * cos(gamma) + thrust * cos(gamma + alphap) + F_wtheta) / (m) - 2 * x_2 * x_4 * sin(x_5) - 2 * x_1 * x_2 * x_6 * cos(x_5)) * (1 / (x_1 * sin(x_5)))
        x_dot_5 = x_6
        x_dot_6 = ((F_wphi / mass) - 2 * x_2 * x_6 + x_1 * x_4**2 * cos(x_5) * sin(x_5)) * (1 / x_1)
        return  x_dot_1, x_dot_2, x_dot_3, x_dot_4, x_dot_5, x_dot_6"""
    
    def windForce(alt):
        #determines wind velocity value for certain altitude (alt)
        altitude = alt
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

        return alt, windVel


    def gusts():
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

    """def sixdofsolver(self, x_1, x_2, x_3, x_4, x_5, x_6, Drag, n_thrust, x, step_size, gravity, gamma, alphap):
        self.thrust = thrust
        self.altitude = altitude
        self.t = step_size
        self.gravity = gravity
        self.gamma = gamma
        self.alphap = alphap
        self.drag = Drag
        
        init_state = [x_1, x_2, x_3, x_4, x_5, x_6]
        sixeq = sixeqdiff(self, fw)
        state = odeint(sixeq, init_state)
        return state"""

def main():
    if __name__ == "__main__":
        cl = SIXDOF

        #sets arrays and open lists
        matrix = cl.gusts()
        alt_input_array = np.array(np.arange(0, 31))
        array_altitude = []
        array_velocity = []
        igroup = []
        
        #creates lists of wind velocity and corresponding altitude values
        for i in alt_input_array:
            altitude, velocity = cl.windForce(int(i))
            array_altitude.append(altitude)
            array_velocity.append(velocity)

        final_array_altitude = np.array(array_altitude)
        final_array_velocity = np.array(array_velocity)
        
        # imports and plots wind gusts   
        for i in matrix:
            for v in i:
                none, wind = cl.windForce(i[1] / 1000)
                mag = i[2] * i[3]
                final_wind = wind + mag
                initial_x = wind / 50
                final_x = final_wind / 50
                if mag > 0:
                    plt.axhline(y = i[1] / 1000, xmin = initial_x, xmax = final_x, color = 'r', linestyle = '-')
                if mag < 0:
                    plt.axhline(y = i[1] / 1000, xmin = final_x, xmax = initial_x, color = 'r', linestyle = '-')
                    
        #plots arrays
        plt.plot(final_array_velocity, final_array_altitude)
        plt.xlabel('Wind Velocity')
        plt.ylabel('Altitude')
        plt.axis([0, 50, 0, 30])
        plt.show()
    
main()
 


