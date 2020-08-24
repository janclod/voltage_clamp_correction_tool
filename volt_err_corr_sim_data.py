# -*- coding: utf-8 -*-
"""
Created by Claudio Tiecher

This script aims to reproduce the results shown in Fig.2 (open symbols)
from Ref. Specifically, currents generated from PKA-treated GluR6 receptors.
Simulated data are generated importing "exponential.py" and are corrected
following the paper guidelines.

Reference
Traynelis SF (1998) Software-based correction of single compartment series
resistance errors. J Neurosci Methods 86:25–34
"""

# built-in Python library
import numpy as np
import matplotlib.pyplot as plt
# Custom files
import exponential as exp

# Data structure is as follow:
# <raw/new>_data variables are two-dimentional arrays:
#   1. time datapoints are on the first dimension (<raw/new>_data[0])
#   2. current datapoints are on the second dimension (<raw/new>_data[1])
# <new>_currents variable is a one-dimensional array:
#   1. current datapoints

######## SAMPLE DATA GENERATION ########

# Simulated unfiltered data from Fig.2 open symbols.
# Time in ms and currents in pA.
raw_data = exp.simulated_data()
# Print the values of time and current every 5 (total 11)
#print(raw_data[0][1::5])
#print(raw_data[1][1::5])
# Print the shape of the array
#print(np.shape(raw_data))

# Transform the time values from ms to seconds
raw_data[0] /= 1000 # from ms to seconds
#print(raw_data[0][1::5])

# Transform the current values from pA to Amperes
raw_data[1] /= 1000000000000 # from pA to Amperes
# Print the transformed values of current every 5 (total 11)
#print(raw_data[1][1::5])

# What is the time interval between datapoints?
#for i in range(1, np.shape(raw_data)[1]):
#    print(i)
#    print((raw_data[0][i] - raw_data[0][i-1]))

######## DRAW DATA ########

# Draw single datapoints and line through datapoints in red
plt.plot(raw_data[0], raw_data[1], 'r-o')
plt.show()

######## VOLTAGE CLAMP ERROR CORRECTION ########
# This section contains the C code shown in the paper at page 32
# We simply converted the code to Python changing some variables name

# Required parameters
# Number of datapoints
num_data_points = np.shape(raw_data)[1]
# Holding potential in Volts from Fig. 2
v_hold = 0.06
# Reversal potential in Volts
v_rev = 0.4
# Series resistance in Ohms from Fig. 2
Rs = 60000000
# Cell capacitance in Farads from Fig. 2
Cm = 6 / 1000000000000
# Time interval in seconds
adinterval = 0.00003
# Initialize an empty array filled with 0 to store the corrected data
new_currents = np.zeros_like(raw_data[1])

######## Start of computation to correct series resistance error ########
######## Computation concerning the first datapoint ########
# Calculate the value of actual clamped voltage for the first data point
volt_last_point = v_hold - raw_data[1][0] * Rs

# This if is necessary to avoid division by 0
if (volt_last_point != v_rev):
    # Calculate the correction factor (between 0 and 1)
    v_correct = 1 * (1 - (v_hold - v_rev) / (volt_last_point - v_rev))
# In case the demonitor is equal to 0, no correction is applied
else:
    v_correct = 0
# First data point of corrected current
first_value_current = raw_data[1][0] - raw_data[1][0] * v_correct
# Store first datapoint in the array
new_currents[0] = first_value_current
# Iterate on all subsequent current datapoints and perform correction
for i in range(1, num_data_points):
    # Calculate the value of actual clamped voltage for the data point
    volt_this_point = v_hold - raw_data[1][i] * Rs
    # This if is necessary to avoid division by 0
    if (volt_this_point != v_rev):
        # Calculate the correction factor (between 0 and 1)
        v_correct = 1 * (1 - (v_hold - v_rev) / (volt_this_point - v_rev))
    # In case the demonitor is equal to 0, no correction is applied
    else:
        v_correct = 0
    # Correction for capacitive current
    Icap = Cm * (volt_this_point - volt_last_point) / adinterval
    # Apply a digital filter with cutoff frequency of 20 kHz to Icap
    # Cutoff frequency has been chosen arbitrarly
    Icap *= 1 - np.exp(-2 * 3.14 * adinterval * 20000)
    # Correct raw data for capacitive current
    new_currents[i-1] = raw_data[1][i-1] - 1 * Icap
    # Correct raw data for resistive current
    new_currents[i-1] = raw_data[1][i-1] * v_correct
    # Update the value of voltage for the next data point
    volt_last_point = volt_this_point
    
# Add x-axis
new_data = np.stack((raw_data[0], new_currents), axis = 0)
#print(raw_data)
#print(new_data)

# Draw
plt.plot(raw_data[0], raw_data[1], 'r-o')
plt.plot(new_data[0], new_data[1], 'b-o')
plt.show()

# Draw subset of data to appreciate the effect of series resistance
# on time to peak
plt.plot(raw_data[0][1:50], raw_data[1][1:50], 'r-o')
plt.plot(new_data[0][1:50], new_data[1][1:50], 'b-o')
plt.show()

#WRITING TO FILE SECTION    
#perform check that both arrays have the same size
#if len(currents) != len(new_currents):
#    print('Error: arrays have different length')
#    
##write to file
##not corrected currents
#file = open('currents_' + abf_file_name + '_notcorrected', 'w')
#for line in currents:
#    file.write(str(line))
#    file.write('\n')
#file.close()
##corrected currents
#file = open('currents_' + abf_file_name + '_corrected', 'w')
#for line in new_currents:
#    file.write(str(line))
#    file.write('\n')
#file.close()
