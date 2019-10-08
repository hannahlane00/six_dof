from math import sin
from atmosphereDensity import seMagic as ad
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import random
import math

def p_vs_a(alt_array):
    p = []
    a = []
    for i in alt_array:
        pressure = press(i)
        p.append(pressure)
        a.append(i)
    p_array = np.array(p)
    a_array = np.array(a)
    x = p_array.sort()
    y = a_array.sort()
    plt.plot(x, y, 'r--')
    return alt_array, p_array

def press(alt):
    pressure = 14.7 * math.exp((-0.366713637542122 * 10**-4) * alt + (-0.000001623765497 * 10**(-4) * alt**2 + (0.000000000007584 * 10**(-4) * alt**3)))
    return pressure
    

if __name__=="__main__":
    alt_values = []
    for i in range(100):
        alt_values.append(random.randrange(0, 150000))
    alt_values_array = np.array(alt_values)
    pvsa(alt_values_array)
    
    plt.xlabel("Pressure")
    plt.ylabel("Altitude")
    plt.axis([0, 15, 0, 150000])
    plt.show()

