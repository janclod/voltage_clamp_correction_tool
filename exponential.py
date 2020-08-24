#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 16:16:28 2020

@author: Claudio Tiecher

This script allows us to generate x-y values for an exponential curve.
"""

import numpy as np

def exponential_rise(tau, amplitude, offset):
    """
    This function generate x-y values of an exponential curve (e.g., current).
    """
    points = [[],[]]
    x_values = np.arange(0.1, 8, 0.5)
    for x in x_values:
        if (x == 0.):
            print("found zero")
        else:
            y = amplitude * np.exp(-tau/x) + offset
            points = np.append(points, [[x],[y]], axis = 1)
    return points

def exponential_decay(tau, amplitude, offset):
    """
    This function generate x-y values of an exponential curve (e.g., current).
    """
    points = [[],[]]
    x_values = np.arange(0.8, 8, 0.025)
    for x in x_values:
        if (x == 0.):
            print("found zero")
        else:
            y = amplitude * np.exp(-tau/x) + offset
        points = np.append(points, [[x],[y]], axis = 1)
    return points

def simulated_data():
    # Rise curve
    values_rise = exponential_rise(1, 1, 0)
    values_rise[0] /= 20
    #subset_values_x = np.array(values_rise[0][2:])
    #subset_values_y = np.array(values_rise[1][2:])
    #values_rise = np.stack((subset_values_x, subset_values_y))
    factor = 1650 / values_rise[1][-1]
    values_rise[1] *= factor
    
    # Decay curve
    values_decay = exponential_decay(2.2, -1.4, 1.4)
    values_decay[0] -= 1.86
    begin_it = 0
    for i in range(0, np.shape(values_decay)[1]):
        if (values_decay[1][i] < 0.87):
            begin_it = i
            break
    subset_values_x = np.array(values_decay[0][begin_it:])
    subset_values_y = np.array(values_decay[1][begin_it:])
    values_decay = np.stack((subset_values_x, subset_values_y))
    values_decay[1] *= factor
    
    # Place all values in one array
    values = np.append(values_rise, values_decay, axis = 1)
    
    return values
