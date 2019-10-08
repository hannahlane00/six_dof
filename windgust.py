# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 13:36:47 2019

@author: jwhribal
"""

# Model of wind gusts

import math 
import numpy as np
import random as randrange
import random

num_gusts = random.randrange(0,20)

alt_range = random.randrange(0,30000,1)
mag_gust = random.randrange(0,12,1)
direction = random.choice((-1,1))

gust_matrix = np.array([0, alt_range, mag_gust, direction])


for i in range(num_gusts-1):
    alt_range = random.randrange(0,30000,1)
    mag_gust = random.randrange(0,12,1)
    direction = random.choice((-1,1))
    temp_array = np.array([i+1, alt_range, mag_gust, direction])
    gust_matrix = np.vstack((temp_array, gust_matrix))
    
gust_matrix = np.flipud(gust_matrix)   
print(gust_matrix)

pyth
