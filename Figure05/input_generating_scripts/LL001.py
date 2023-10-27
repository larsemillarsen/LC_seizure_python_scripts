# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 10:19:00 2022

@author: llarsen
"""

import numpy as np
import sys
sys.path.insert(1, r'E:\Manuscript_analysis_files')

from vies.spike.extract_frequency_bins import extract_frequency_bins
import h5py
from scipy.interpolate import griddata


## DATA FOR PINCH
MUA_filename = r'E:\Manuscript_analysis_files\data\mua_data\LL001.mua.hdf5'
probe_filename = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\small_poly3.prb'
messagefile = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\LL001.events'
LC_mask_file = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\LC_masks\LC_masks\LL001_masktest.npy'

n_pinches = 10 # number of pinch trials
trial_length = 25 # length of pinch trials in seconds
srate = 30000 # samples per second
bin_size = 1 # bin size in seconds

## STEP0: Load pinch timestamps, fill in path to time stamp file from open ephys below
f=open(messagefile, 'r')
events = f.read()
f.close()

events=events.split('pinch')
events=events[3:15]
events = [int(x) for x in events]
events = np.array(events)/srate
adjustments = np.array([-0.4, -0.45, -0.5, -0.45, -0.75, -0.45, -0.3, -0.55, -0.45, -0.5, -0.4, -0.25])
events = events - adjustments

n_pinches = len(events) # number of pinch trials

## STEP1: Load MUA data
with h5py.File(MUA_filename, "r") as f:
    n_channels = len(list(f['spiketimes'].keys()))
    labels = np.zeros(n_channels)
    ids = list(f['spiketimes'].keys())
    
    for i in range(n_channels):
        labels[i] = int(ids[i].split('elec_')[1])
        #print(labels[i])
    
    labels=np.array(labels)        
    labels=labels.astype(int)
    
    n_bins = int(trial_length / bin_size)
    
    mua_matrix = np.zeros((n_bins, len(labels), n_pinches))

    counter = 0
    for i in labels:
        index = int(i)
        e_id = ids[counter] 
        #print(index, e_id)
        spiketimes = (f['spiketimes'][e_id][:]) / srate # recalculates spiketimes in seconds
        for k in range(n_pinches):
            start_bin = events[k] - 20
            stop_bin = events[k] + 5
            time, freq = extract_frequency_bins(spiketimes, bin_size, start_bin, stop_bin)
            if len(freq) == 1199:
                freq = np.append(freq, 0)
            mua_matrix[:,index, k] = freq
        counter += 1

original_muamatrix = mua_matrix
mua_matrix = np.mean(mua_matrix, axis=2)
n_baseline_bins = 19
base_mean_freq = np.mean(mua_matrix[0:n_baseline_bins,:], axis=0)
base_std_freq = np.std(mua_matrix[0:n_baseline_bins,:], axis=0)

z_scored_matrix = np.zeros((np.shape(mua_matrix)[0], np.shape(mua_matrix)[1]))
for i in range(np.shape(mua_matrix)[0]):
    z_scored_matrix[i,:] = (mua_matrix[i,:] - base_mean_freq) / base_std_freq

mua_matrix = z_scored_matrix

## STEP2: Load probe geometry
probe = {}
with open(probe_filename, 'r') as f:
            probetext = f.read()
            exec(probetext, probe)
            del probe['__builtins__']
            
geometry = probe['channel_groups'][1]['geometry']

x = np.zeros(len(geometry))
y = np.zeros(len(geometry))

for i in range(len(geometry)):
    x[i] = int(geometry[i][0])
    y[i] = int(geometry[i][1])
    
dead_channels = []

x = np.delete(x, dead_channels)
y = np.abs(np.delete(y, dead_channels))

pinch = 19
z1 = np.mean(mua_matrix[(pinch-3):(pinch-1),:],axis=0)
z2 = mua_matrix[pinch,:]
z3 = mua_matrix[(pinch+1),:]

xi = np.arange(-25,26,5)
yi = np.arange(0, 276,5)
xi,yi = np.meshgrid(xi,yi)

zi1 = griddata((x,y),z1,(xi,yi),method='cubic')
zi2 = griddata((x,y),z2,(xi,yi),method='cubic')
zi3 = griddata((x,y),z3,(xi,yi),method='cubic')


threshold = 3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z2)):
    if z2[i]>threshold and z3[i]<-threshold:
        counter = counter + 1

significance_matrix = np.zeros((counter, 4))

counter_two = 0
for i in range(len(z2)):
    if z2[i]>threshold and z3[i]<-threshold:
        significance_matrix[counter_two,0] = x[i]
        significance_matrix[counter_two,1] = y[i]
        significance_matrix[counter_two,2] = z2[i]        
        significance_matrix[counter_two,3] = i
        
        counter_two = counter_two + 1 

np.save(LC_mask_file, significance_matrix)
