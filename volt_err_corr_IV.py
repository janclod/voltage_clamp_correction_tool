# -*- coding: utf-8 -*-
"""
Created by Claudio Tiecher

This script contains several functions to facilitate the analysis of
whole-cell electrophysiology data and perform offline correction
of voltage clamp error due to series resistance and cell capacitance.

Reference
Traynelis SF (1998) Software-based correction of single compartment series
resistance errors. J Neurosci Methods 86:25–34
"""

import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from pathlib import Path

from neo.io import AxonIO #import data structure for .ABF files

#def get_vals_curr(seg):
#    '''
#    This function return a list of currents values.
#    It assumes that the first analog signal (#0) is current in pA
#    '''
#    #gets the length of segment
#    l_sig = len(seg.analogsignals[0])
#    #stores the first analog signal
#    a_s = seg.analogsignals[0]
#    #initialize variable to store curent values
#    curr_vals = []
#    for i in range(0, l_sig):
#        curr_vals.append(a_s[i][0])
#    
#    return curr_vals
#

def get_segments(abf_obj, b_ind = 0):
    '''
    Collects the segments.
    It assumes that the AxonIO contains only 1 block.
    '''
    #get the number of sgments in abf_file
    n_seg = abf_obj.segment_count(block_index = b_ind)
    #list of segment objects
    list_of_segments = []
    #extract one sweep and appends it to the list of segments
    for i in range(0, n_seg):
        list_of_segments.append(abf_obj.read_segment(seg_index = i))
        
    return list_of_segments

def t_start_vstep(segment):
    '''
    This function works on a segment with 2 AnalogSignals.
    It assumes that:
        (1) we want to work on the second analog signal
        (2) the second analog signal contains a voltage protocol
        (3) the voltage steps are bigger than 3 mV
    The function returns the index at which the voltage steps up.
    If there is no voltage step it returns sample = 0.
    It is only able to detect the timepoint for the first voltage step.
    '''
    #initialize sample
    sample = 0
    #initalialize the first value of voltage for comparison
    v_prev = segment.analogsignals[1][0]
    for i in range(len(segment.analogsignals[1])):
        if segment.analogsignals[1][i] > v_prev + 3 * pq.mV:
            sample = i
            return sample
    return sample

#def t_end_vstep_down(segment, t_start):
#    '''
#    This function works on a segment with 2 AnalogSignals.
#    It assumes that:
#        (1) we want to work on the second analog signal
#        (2) the second analog signal contains a voltage protocol
#        (3) the voltage steps are bigger than 3 mV
#        (4) we have already found the starting timepoint for the voltage step
#    The function returns the index at which the voltage steps down.
#    If there is no voltage step it returns index = 0.
#    It is only able to detect the timepoint for the first voltage step.
#    '''
#    index = 0
#    v_prev = segment.analogsignals[1][t_start]
#    for i in range(t_start, len(segment.analogsignals[1])):
#        if segment.analogsignals[1][i] < v_prev - 3 * pq.mV:
#            index = i
#            return index
#    return index

def t_end_vstep_up(segment, t_start):
    '''
    This function works on a segment with 2 AnalogSignals.
    It assumes that:
        (1) we want to work on the second analog signal
        (2) the second analog signal contains a voltage protocol
        (3) the voltage steps are bigger than 3 mV
        (4) we have already found the starting timepoint for the voltage step
    The function returns the index at which the voltage steps up.
    If there is no voltage step it returns sample = 0.
    '''
    sample = 0
    v_prev = segment.analogsignals[1][t_start]
    for i in range(t_start, len(segment.analogsignals[1])):
        if segment.analogsignals[1][i] > v_prev + 3 * pq.mV:
            sample = i
            return sample
    return sample

def average_current(seg, t):
    '''
    This function returns the average value of currents over 6.6 ms
    t is the timepoint at the end of the range
    It assumes that the sampling rate is 33 kHz
    '''
    #t1 and t2
    t1 = t - 200
    t2 = t
    #variable to store summed values
    s = 0 * pq.pA
    for i in range(t1, t2):
        s += seg.analogsignals[0][i][0]
    #average current value
    avg_curr = s / (t2 - t1)
    
    return avg_curr

def average_current_short(a):
    '''
    This function finds the max value within
    a subset array starting at index 50.
    It returns an average value based on 10
    values found around the max value.
    '''
    shift = 550
    v_i = np.argmax(a[shift:])
    s = 0
    c = 0
    for i in range(v_i + shift -5, v_i + shift + 5):
         s += a[i]
         c += 1
    average = s / c
    return average
    #return a[v_i + shift]

def est_v_rev(abf_obj, segs, t):
    '''
    This function extrpolate the reversal potential from several sweeps
    '''
    #Initialize v_rev
    v_rev = -100000 * pq.V
    #Extract currents
    #List of current values
    currents = []
    #extract currents values from each segment and store values in list
    for i in range(0, abf_obj.segment_count(block_index = 0)):
        curr = average_current(segs[i], t)
        currents.append(curr)
        
    #range of voltages
    v_values = np.arange(-90, 70, 10)
    
    #check
    print(currents)
    print(v_values)
    
    #linear regression
    fit_par = np.polyfit(v_values, currents, 1)
    #estimate v_rev
    v_rev = - fit_par[1] / fit_par[0] / 1000
    
    return v_rev

def get_max_seg_value_from_analog0(list_of_segments):
    '''
    Need to write DOC for this functon
    '''
    max_values = []
    for s_n in range(len(list_of_segments)):
        max_y = -100000
        for t_p in range(len(list_of_segments[s_n].analogsignals[0])):
            if list_of_segments[s_n].analogsignals[0][t_p] > max_y:
                max_y = list_of_segments[s_n].analogsignals[0][t_p]
            t_p = t_p +1
        max_values.append(max_y)
    #print('The length of the seg is ' + str(len(segs[s_n].analogsignals[0])) + ' for current signal')
    #print('The max current is ' + str(max_y))
    return max_values

abf_file_name = './data/input/18227014.abf'
abf = AxonIO(abf_file_name)

#print file info
print(abf.read_protocol)

#count segments
n_seg = abf.segment_count(block_index = 0)
print('Found ' + str(n_seg) + ' segments.')

# Generate an array containing segments
segs = get_segments(abf)
print(segs[15])
# Visually check the values of the analog signal
# in seg[15].
print(segs[15].analogsignals[0])
a_seg_1 = np.array(segs[15].analogsignals[0])
plt.plot(a_seg_1, 'r-o')
plt.show()

# Get the time point at which the voltage step
# t1_step = t_start_vstep(segs[15])
t1_step = 18000
print(t1_step)
# t2_step = t_end_vstep_up(segs[15], t1_step)
t2_step = 80000
print(t2_step)

# Iterate over slice
currents = np.array([]) * pq.pA
for i in range(t1_step, t2_step):
    value = segs[15].analogsignals[0][i][0] * pq.pA
    currents = np.append(currents, value) * pq.pA
print(len(range(t1_step,t2_step)))
print(currents[0])
print(currents[-1])

# Visually check the values of currents.
plt.plot(currents, 'r-o')
plt.show()

####
# Correction of voltage clamp error
####

# Required parameters
# Holding voltage in Volts
# ISSUE: the sign of the v_hold seems to have an effect
# on the data. This needs to be understood.
# To generate issue simply put a minus in front of the v_hold value
# and look at the plot of currents before and after
# voltage clamp error correction.
v_hold = 0.06
# Reversal potential in Volts
v_rev = est_v_rev(abf, segs, t2_step)
print("This is " + str(v_rev))
# Series resistance in Ohms
Rs = 30 * 10**6
# Cell capacitance in Farads
Cm = 3.4 * 10**-12
# Sampling rate
# sr = 33333 * pq.Hz                      
# Time interval between samples in seconds
adinterval = 0.00003
# Inverse of smapling rate
# tau_lag = 1 / sr

# Get dimensionless value of currents
f_currents = np.zeros(currents.shape)
for i in range(len(currents)):
    f_currents[i] = float(currents[i])
print(f_currents[0:5])
# Convert from pA to Ampere
f_currents /= 1000000000000
# Visually check the values of currents.
plt.plot(f_currents, 'r-o')
plt.show()


# Initialize an empty array filled with 0 to store the corrected data
new_currents = np.zeros(currents.shape)

######## Start of computation to correct series resistance error ########
######## Computation concerning the first datapoint ########
# Calculate the value of actual clamped voltage for the first data point
volt_last_point = (-v_hold) - f_currents[0] * Rs
print(f_currents[0])

# This if is necessary to avoid division by 0
if (volt_last_point != v_rev):
    # Calculate the correction factor (between 0 and 1)
    v_correct = 1 * (1 - ((-v_hold) - v_rev) / (volt_last_point - v_rev))
# In case the demonitor is equal to 0, no correction is applied
else:
    v_correct = 0
# First data point of corrected current
first_value_current = f_currents[0] - f_currents[0] * v_correct
# Store first datapoint in the array
new_currents[0] = first_value_current
# Iterate on all subsequent current datapoints and perform correction
for i in range(1, len(currents)):
    # Calculate the value of actual clamped voltage for the data point
    volt_this_point = (-v_hold) - f_currents[i] * Rs
    # This if is necessary to avoid division by 0
    if (volt_this_point != v_rev):
        # Calculate the correction factor (between 0 and 1)
        v_correct = 1 * (1 - ((-v_hold) - v_rev) / (volt_this_point - v_rev))
    # In case the demonitor is equal to 0, no correction is applied
    else:
        v_correct = 0
    # Correction for capacitive current
    Icap = Cm * (volt_this_point - volt_last_point) / adinterval
    # Apply a digital filter with cutoff frequency of 20 kHz to Icap
    # Cutoff frequency has been chosen arbitrarly
    Icap *= 1 - np.exp(-2 * 3.14 * adinterval * 20000)
    # Correct raw data for capacitive current
    new_currents[i-1] = f_currents[i-1] - 1 * Icap
    # Correct raw data for resistive current
    new_currents[i-1] = f_currents[i-1] * v_correct
    # Update the value of voltage for the next data point
    volt_last_point = volt_this_point

print(new_currents[0:10])

# Draw
# Zoomin
ax1 = plt.subplot(111)
ax1.set_ylim(-0.000000003, 0.0000000015)
ax1.set_xlim(0, 6000)
ax1.plot(f_currents, 'r')
# plt.show()
ax2 = plt.subplot(111)
ax2.set_ylim(-0.000000003, 0.0000000015)
ax2.set_xlim(0, 6000)
ax2.plot(new_currents, 'b')
plt.axhline(5.9 * 10 **-12)
plt.axhline(1.1 * 10 **-10)
plt.axvline(600)
plt.show()
# Zoomout
ax1 = plt.subplot(111)
ax1.plot(f_currents, 'r')
ax2 = plt.subplot(111)
ax2.plot(new_currents, 'b')
plt.show()


#WRITING TO FILE SECTION    
#perform check that both arrays have the same size
write = True
if write:
    # Generate file name string
    p = Path(abf_file_name)
    file_name = p.name
    if len(currents) != len(new_currents):
        print('Error: arrays have different length')
        
    #write to file
    #not corrected currents
    file = open('./data/output/currents_' + file_name[:-4] + '_notcorrected', 'w')
    file.write(str(average_current_short(f_currents)))
    file.close()
    #corrected currents
    file = open('./data/output/currents_' + file_name[:-4] + '_corrected', 'w')
    file.write(str(average_current_short(new_currents)))
    file.close()
