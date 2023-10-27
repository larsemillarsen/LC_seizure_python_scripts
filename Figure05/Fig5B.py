# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 14:26:40 2022

@author: llarsen
"""

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys
sys.path.insert(1, r'E:\Manuscript_analysis_files')

from vies.spike.extract_frequency_bins import extract_frequency_bins


### General Params

trial_length = 120 # length of pinch trials in seconds
srate = 30000 # samples per second

bin_size = 5 # bin size in seconds
n_bins = int(360 / bin_size)
n_baseline_bins = int(300/bin_size)

Nt = int(trial_length/bin_size) # number of frames for animation
small_probe_filename = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\small_poly3.prb'
big_probe_filename = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\big_poly3.prb'

savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\output'


### LOAD small probe geometry
probe = {}
with open(small_probe_filename, 'r') as f:
            probetext = f.read()
            exec(probetext, probe)
            del probe['__builtins__']
            
geometry = probe['channel_groups'][1]['geometry']

x_small = np.zeros(len(geometry))
y_small = np.zeros(len(geometry))

for i in range(len(geometry)):
    x_small[i] = int(geometry[i][0])
    y_small[i] = int(geometry[i][1])
    
dead_channels = []

x_small = np.delete(x_small, dead_channels)
y_small = np.abs(np.delete(y_small, dead_channels))


xi_small = np.arange(-25,26,5)
yi_small = np.arange(0, 276,5)
xi_small,yi_small = np.meshgrid(xi_small,yi_small)

### LOAD big probe geometry
probe = {}
with open(big_probe_filename, 'r') as f:
            probetext = f.read()
            exec(probetext, probe)
            del probe['__builtins__']
            
geometry = probe['channel_groups'][1]['geometry']

x_big = np.zeros(len(geometry))
y_big = np.zeros(len(geometry))

for i in range(len(geometry)):
    x_big[i] = int(geometry[i][0])
    y_big[i] = int(geometry[i][1])
    
dead_channels = []

x_big = np.delete(x_big, dead_channels)
y_big = np.abs(np.delete(y_big, dead_channels))


xi_big = np.arange(-50,51,10)
yi_big = np.arange(0, 551,10)
xi_big,yi_big = np.meshgrid(xi_big,yi_big)


### LL001
LC_mask_file = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\LC_masks\LL001_masktest.npy'
MUA_filename = r'E:\Manuscript_analysis_files\data\mua_data\LL001.mua.hdf5'
sz1_start = 3782.72


with h5py.File(MUA_filename, "r") as f:
    n_channels = len(list(f['spiketimes'].keys()))
    labels = np.zeros(n_channels)
    ids = list(f['spiketimes'].keys())

    for i in range(n_channels):
        labels[i] = int(ids[i].split('elec_')[1])
    
    labels=np.array(labels)        
    labels=labels.astype(int)
  
    LL001 = np.zeros((n_bins, len(labels)))

    counter = 0
    for i in labels:
        index = int(i)
        e_id = ids[counter] 
        spiketimes = (f['spiketimes'][e_id][:]) / srate # recalculates spiketimes in seconds
        start_bin = sz1_start-300
        stop_bin = start_bin+360
        time, freq = extract_frequency_bins(spiketimes, bin_size, start_bin, stop_bin)
        LL001[:,index] = freq
        
        counter += 1

base_mean_freq = np.mean(LL001[0:n_baseline_bins,:,], axis=0)
base_std_freq = np.std(LL001[0:n_baseline_bins,:,], axis=0)
z_scored_matrix = np.zeros((np.shape(LL001)[0], np.shape(LL001)[1]))
for i in range(np.shape(LL001)[0]):
    z_scored_matrix[i,:] = (LL001[i,:] - base_mean_freq) / base_std_freq

z_LL001 = z_scored_matrix[(n_baseline_bins+2), :]
zi_LL001 = griddata((x_small,y_small),z_LL001,(xi_small,yi_small),method='cubic')



lc_mask_LL001 = np.load(LC_mask_file)

threshold = 3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL001)):
    if z_LL001[i]>threshold:
        counter = counter + 1

significance_matrix_z_LL001 = np.zeros((counter, 4))
counter_two = 0
for i in range(len(z_LL001)):
    if z_LL001[i]>threshold:
        significance_matrix_z_LL001[counter_two,0] = x_small[i]
        significance_matrix_z_LL001[counter_two,1] = y_small[i]
        significance_matrix_z_LL001[counter_two,2] = z_LL001[i]        
        significance_matrix_z_LL001[counter_two,3] = i
        
        counter_two = counter_two + 1 
        
threshold = -3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL001)):
    if z_LL001[i]<threshold:
        counter = counter + 1

significance_matrix_two_z_LL001 = np.zeros((counter, 4))

counter_two = 0
for i in range(len(z_LL001)):
    if z_LL001[i]<threshold:
        significance_matrix_two_z_LL001[counter_two,0] = x_small[i]
        significance_matrix_two_z_LL001[counter_two,1] = y_small[i]
        significance_matrix_two_z_LL001[counter_two,2] = z_LL001[i]        
        significance_matrix_two_z_LL001[counter_two,3] = i  
        
        counter_two = counter_two + 1

### sz1
counter = 0
for i in range(np.shape(significance_matrix_z_LL001)[0]):
    if significance_matrix_z_LL001[i,3] in lc_mask_LL001[:,3]:
        counter = counter + 1

significance_matrix_LC_z_LL001 = np.zeros((counter, 4))

c = 0
for i in range(np.shape(significance_matrix_z_LL001)[0]):
    if significance_matrix_z_LL001[i,3] in lc_mask_LL001[:,3]:
        significance_matrix_LC_z_LL001[c,:] = significance_matrix_z_LL001[i,:]
        c = c + 1
        
counter_two = 0
for i in range(np.shape(significance_matrix_two_z_LL001)[0]):
    if significance_matrix_two_z_LL001[i,3] in lc_mask_LL001[:,3]:
        counter_two = counter_two + 1

significance_matrix_two_LC_z_LL001 = np.zeros((counter_two, 4))

c = 0
for i in range(np.shape(significance_matrix_two_z_LL001)[0]):
    if significance_matrix_two_z_LL001[i,3] in lc_mask_LL001[:,3]:
        significance_matrix_two_LC_z_LL001[c,:] = significance_matrix_two_z_LL001[i,:]
        c = c + 1

### LL009
LC_mask_file = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\LC_masks\LL009.npy'
MUA_filename = r'E:\Manuscript_analysis_files\data\mua_data\LL009.mua.hdf5'
sz1_start = 4408.35


with h5py.File(MUA_filename, "r") as f:
    n_channels = len(list(f['spiketimes'].keys()))
    labels = np.zeros(n_channels)
    ids = list(f['spiketimes'].keys())

    for i in range(n_channels):
        labels[i] = int(ids[i].split('elec_')[1])
    
    labels=np.array(labels)        
    labels=labels.astype(int)
  
    LL009 = np.zeros((n_bins, len(labels)))

    counter = 0
    for i in labels:
        index = int(i)
        e_id = ids[counter] 
        spiketimes = (f['spiketimes'][e_id][:]) / srate # recalculates spiketimes in seconds
        start_bin = sz1_start-300
        stop_bin = start_bin+360
        time, freq = extract_frequency_bins(spiketimes, bin_size, start_bin, stop_bin)
        LL009[:,index] = freq
        
        counter += 1

base_mean_freq = np.mean(LL009[0:n_baseline_bins,:,], axis=0)
base_std_freq = np.std(LL009[0:n_baseline_bins,:,], axis=0)
z_scored_matrix = np.zeros((np.shape(LL009)[0], np.shape(LL009)[1]))
for i in range(np.shape(LL009)[0]):
    z_scored_matrix[i,:] = (LL009[i,:] - base_mean_freq) / base_std_freq

z_LL009 = z_scored_matrix[(n_baseline_bins+2), :]
zi_LL009 = griddata((x_small,y_small),z_LL009,(xi_small,yi_small),method='cubic')



lc_mask_LL009 = np.load(LC_mask_file)

threshold = 3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL009)):
    if z_LL009[i]>threshold:
        counter = counter + 1

significance_matrix_z_LL009 = np.zeros((counter, 4))
counter_two = 0
for i in range(len(z_LL009)):
    if z_LL009[i]>threshold:
        significance_matrix_z_LL009[counter_two,0] = x_small[i]
        significance_matrix_z_LL009[counter_two,1] = y_small[i]
        significance_matrix_z_LL009[counter_two,2] = z_LL009[i]        
        significance_matrix_z_LL009[counter_two,3] = i
        
        counter_two = counter_two + 1 
        
threshold = -3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL009)):
    if z_LL009[i]<threshold:
        counter = counter + 1

significance_matrix_two_z_LL009 = np.zeros((counter, 4))

counter_two = 0
for i in range(len(z_LL009)):
    if z_LL009[i]<threshold:
        significance_matrix_two_z_LL009[counter_two,0] = x_small[i]
        significance_matrix_two_z_LL009[counter_two,1] = y_small[i]
        significance_matrix_two_z_LL009[counter_two,2] = z_LL009[i]        
        significance_matrix_two_z_LL009[counter_two,3] = i  
        
        counter_two = counter_two + 1

### sz1
counter = 0
for i in range(np.shape(significance_matrix_z_LL009)[0]):
    if significance_matrix_z_LL009[i,3] in lc_mask_LL009[:,3]:
        counter = counter + 1

significance_matrix_LC_z_LL009 = np.zeros((counter, 4))

c = 0
for i in range(np.shape(significance_matrix_z_LL009)[0]):
    if significance_matrix_z_LL009[i,3] in lc_mask_LL009[:,3]:
        significance_matrix_LC_z_LL009[c,:] = significance_matrix_z_LL009[i,:]
        c = c + 1
        
counter_two = 0
for i in range(np.shape(significance_matrix_two_z_LL009)[0]):
    if significance_matrix_two_z_LL009[i,3] in lc_mask_LL009[:,3]:
        counter_two = counter_two + 1

significance_matrix_two_LC_z_LL009 = np.zeros((counter_two, 4))

c = 0
for i in range(np.shape(significance_matrix_two_z_LL009)[0]):
    if significance_matrix_two_z_LL009[i,3] in lc_mask_LL009[:,3]:
        significance_matrix_two_LC_z_LL009[c,:] = significance_matrix_two_z_LL009[i,:]
        c = c + 1
        
## LL017
LC_mask_file = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\LC_masks\LL017.npy'
MUA_filename = r'E:\Manuscript_analysis_files\data\mua_data\LL017.mua.hdf5'
sz1_start = 3621.83


with h5py.File(MUA_filename, "r") as f:
    n_channels = len(list(f['spiketimes'].keys()))
    labels = np.zeros(n_channels)
    ids = list(f['spiketimes'].keys())

    for i in range(n_channels):
        labels[i] = int(ids[i].split('elec_')[1])
    
    labels=np.array(labels)        
    labels=labels.astype(int)
  
    LL017 = np.zeros((n_bins, len(labels)))

    counter = 0
    for i in labels:
        index = int(i)
        e_id = ids[counter] 
        spiketimes = (f['spiketimes'][e_id][:]) / srate # recalculates spiketimes in seconds
        start_bin = sz1_start-300
        stop_bin = start_bin+360
        time, freq = extract_frequency_bins(spiketimes, bin_size, start_bin, stop_bin)
        LL017[:,index] = freq
        
        counter += 1

base_mean_freq = np.mean(LL017[0:n_baseline_bins,:,], axis=0)
base_std_freq = np.std(LL017[0:n_baseline_bins,:,], axis=0)
z_scored_matrix = np.zeros((np.shape(LL017)[0], np.shape(LL017)[1]))
for i in range(np.shape(LL017)[0]):
    z_scored_matrix[i,:] = (LL017[i,:] - base_mean_freq) / base_std_freq

z_LL017 = z_scored_matrix[(n_baseline_bins+2), :]
zi_LL017 = griddata((x_small,y_small),z_LL017,(xi_small,yi_small),method='cubic')



lc_mask_LL017 = np.load(LC_mask_file)

threshold = 3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL017)):
    if z_LL017[i]>threshold:
        counter = counter + 1

significance_matrix_z_LL017 = np.zeros((counter, 4))
counter_two = 0
for i in range(len(z_LL017)):
    if z_LL017[i]>threshold:
        significance_matrix_z_LL017[counter_two,0] = x_small[i]
        significance_matrix_z_LL017[counter_two,1] = y_small[i]
        significance_matrix_z_LL017[counter_two,2] = z_LL017[i]        
        significance_matrix_z_LL017[counter_two,3] = i
        
        counter_two = counter_two + 1 
        
threshold = -3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL017)):
    if z_LL017[i]<threshold:
        counter = counter + 1

significance_matrix_two_z_LL017 = np.zeros((counter, 4))

counter_two = 0
for i in range(len(z_LL017)):
    if z_LL017[i]<threshold:
        significance_matrix_two_z_LL017[counter_two,0] = x_small[i]
        significance_matrix_two_z_LL017[counter_two,1] = y_small[i]
        significance_matrix_two_z_LL017[counter_two,2] = z_LL017[i]        
        significance_matrix_two_z_LL017[counter_two,3] = i  
        
        counter_two = counter_two + 1

### sz1
counter = 0
for i in range(np.shape(significance_matrix_z_LL017)[0]):
    if significance_matrix_z_LL017[i,3] in lc_mask_LL017[:,3]:
        counter = counter + 1

significance_matrix_LC_z_LL017 = np.zeros((counter, 4))

c = 0
for i in range(np.shape(significance_matrix_z_LL017)[0]):
    if significance_matrix_z_LL017[i,3] in lc_mask_LL017[:,3]:
        significance_matrix_LC_z_LL017[c,:] = significance_matrix_z_LL017[i,:]
        c = c + 1
        
counter_two = 0
for i in range(np.shape(significance_matrix_two_z_LL017)[0]):
    if significance_matrix_two_z_LL017[i,3] in lc_mask_LL017[:,3]:
        counter_two = counter_two + 1

significance_matrix_two_LC_z_LL017 = np.zeros((counter_two, 4))

c = 0
for i in range(np.shape(significance_matrix_two_z_LL017)[0]):
    if significance_matrix_two_z_LL017[i,3] in lc_mask_LL017[:,3]:
        significance_matrix_two_LC_z_LL017[c,:] = significance_matrix_two_z_LL017[i,:]
        c = c + 1
        


## LL023
LC_mask_file = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\LC_masks\LL023.npy'
MUA_filename = r'E:\Manuscript_analysis_files\data\mua_data\LL023.mua.hdf5'
sz1_start = 2753.56


with h5py.File(MUA_filename, "r") as f:
    n_channels = len(list(f['spiketimes'].keys()))
    labels = np.zeros(n_channels)
    ids = list(f['spiketimes'].keys())

    for i in range(n_channels):
        labels[i] = int(ids[i].split('elec_')[1])
    
    labels=np.array(labels)        
    labels=labels.astype(int)
  
    LL023 = np.zeros((n_bins, len(labels)))

    counter = 0
    for i in labels:
        index = int(i)
        e_id = ids[counter] 
        spiketimes = (f['spiketimes'][e_id][:]) / srate # recalculates spiketimes in seconds
        start_bin = sz1_start-300
        stop_bin = start_bin+360
        time, freq = extract_frequency_bins(spiketimes, bin_size, start_bin, stop_bin)
        LL023[:,index] = freq
        
        counter += 1

base_mean_freq = np.mean(LL023[0:n_baseline_bins,:,], axis=0)
base_std_freq = np.std(LL023[0:n_baseline_bins,:,], axis=0)
z_scored_matrix = np.zeros((np.shape(LL023)[0], np.shape(LL023)[1]))
for i in range(np.shape(LL023)[0]):
    z_scored_matrix[i,:] = (LL023[i,:] - base_mean_freq) / base_std_freq

z_LL023 = z_scored_matrix[(n_baseline_bins+2), :]
zi_LL023 = griddata((x_small,y_small),z_LL023,(xi_small,yi_small),method='cubic')



lc_mask_LL023 = np.load(LC_mask_file)

threshold = 3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL023)):
    if z_LL023[i]>threshold:
        counter = counter + 1

significance_matrix_z_LL023 = np.zeros((counter, 4))
counter_two = 0
for i in range(len(z_LL023)):
    if z_LL023[i]>threshold:
        significance_matrix_z_LL023[counter_two,0] = x_small[i]
        significance_matrix_z_LL023[counter_two,1] = y_small[i]
        significance_matrix_z_LL023[counter_two,2] = z_LL023[i]        
        significance_matrix_z_LL023[counter_two,3] = i
        
        counter_two = counter_two + 1 
        
threshold = -3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL023)):
    if z_LL023[i]<threshold:
        counter = counter + 1

significance_matrix_two_z_LL023 = np.zeros((counter, 4))

counter_two = 0
for i in range(len(z_LL023)):
    if z_LL023[i]<threshold:
        significance_matrix_two_z_LL023[counter_two,0] = x_small[i]
        significance_matrix_two_z_LL023[counter_two,1] = y_small[i]
        significance_matrix_two_z_LL023[counter_two,2] = z_LL023[i]        
        significance_matrix_two_z_LL023[counter_two,3] = i  
        
        counter_two = counter_two + 1

### sz1
counter = 0
for i in range(np.shape(significance_matrix_z_LL023)[0]):
    if significance_matrix_z_LL023[i,3] in lc_mask_LL023[:,3]:
        counter = counter + 1

significance_matrix_LC_z_LL023 = np.zeros((counter, 4))

c = 0
for i in range(np.shape(significance_matrix_z_LL023)[0]):
    if significance_matrix_z_LL023[i,3] in lc_mask_LL023[:,3]:
        significance_matrix_LC_z_LL023[c,:] = significance_matrix_z_LL023[i,:]
        c = c + 1
        
counter_two = 0
for i in range(np.shape(significance_matrix_two_z_LL023)[0]):
    if significance_matrix_two_z_LL023[i,3] in lc_mask_LL023[:,3]:
        counter_two = counter_two + 1

significance_matrix_two_LC_z_LL023 = np.zeros((counter_two, 4))

c = 0
for i in range(np.shape(significance_matrix_two_z_LL023)[0]):
    if significance_matrix_two_z_LL023[i,3] in lc_mask_LL023[:,3]:
        significance_matrix_two_LC_z_LL023[c,:] = significance_matrix_two_z_LL023[i,:]
        c = c + 1
        



## LL007
LC_mask_file = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\LC_masks\LL007_masktest.npy'
MUA_filename = r'E:\Manuscript_analysis_files\data\mua_data\LL007.mua.hdf5'
sz1_start = 1193.42


with h5py.File(MUA_filename, "r") as f:
    n_channels = len(list(f['spiketimes'].keys()))
    labels = np.zeros(n_channels)
    ids = list(f['spiketimes'].keys())

    for i in range(n_channels):
        labels[i] = int(ids[i].split('elec_')[1])
    
    labels=np.array(labels)        
    labels=labels.astype(int)
  
    LL007 = np.zeros((n_bins, len(labels)))

    counter = 0
    for i in labels:
        index = int(i)
        e_id = ids[counter] 
        spiketimes = (f['spiketimes'][e_id][:]) / srate # recalculates spiketimes in seconds
        start_bin = sz1_start-300
        stop_bin = start_bin+360
        time, freq = extract_frequency_bins(spiketimes, bin_size, start_bin, stop_bin)
        LL007[:,index] = freq
        
        counter += 1

base_mean_freq = np.mean(LL007[0:n_baseline_bins,:,], axis=0)
base_std_freq = np.std(LL007[0:n_baseline_bins,:,], axis=0)
z_scored_matrix = np.zeros((np.shape(LL007)[0], np.shape(LL007)[1]))
for i in range(np.shape(LL007)[0]):
    z_scored_matrix[i,:] = (LL007[i,:] - base_mean_freq) / base_std_freq

z_LL007 = z_scored_matrix[(n_baseline_bins+2), :]
zi_LL007 = griddata((x_big,y_big),z_LL007,(xi_big,yi_big),method='cubic')



lc_mask_LL007 = np.load(LC_mask_file)

threshold = 3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL007)):
    if z_LL007[i]>threshold:
        counter = counter + 1

significance_matrix_z_LL007 = np.zeros((counter, 4))
counter_two = 0
for i in range(len(z_LL007)):
    if z_LL007[i]>threshold:
        significance_matrix_z_LL007[counter_two,0] = x_big[i]
        significance_matrix_z_LL007[counter_two,1] = y_big[i]
        significance_matrix_z_LL007[counter_two,2] = z_LL007[i]        
        significance_matrix_z_LL007[counter_two,3] = i
        
        counter_two = counter_two + 1 
        
threshold = -3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL007)):
    if z_LL007[i]<threshold:
        counter = counter + 1

significance_matrix_two_z_LL007 = np.zeros((counter, 4))

counter_two = 0
for i in range(len(z_LL007)):
    if z_LL007[i]<threshold:
        significance_matrix_two_z_LL007[counter_two,0] = x_big[i]
        significance_matrix_two_z_LL007[counter_two,1] = y_big[i]
        significance_matrix_two_z_LL007[counter_two,2] = z_LL007[i]        
        significance_matrix_two_z_LL007[counter_two,3] = i  
        
        counter_two = counter_two + 1

### sz1
counter = 0
for i in range(np.shape(significance_matrix_z_LL007)[0]):
    if significance_matrix_z_LL007[i,3] in lc_mask_LL007[:,3]:
        counter = counter + 1

significance_matrix_LC_z_LL007 = np.zeros((counter, 4))

c = 0
for i in range(np.shape(significance_matrix_z_LL007)[0]):
    if significance_matrix_z_LL007[i,3] in lc_mask_LL007[:,3]:
        significance_matrix_LC_z_LL007[c,:] = significance_matrix_z_LL007[i,:]
        c = c + 1
        
counter_two = 0
for i in range(np.shape(significance_matrix_two_z_LL007)[0]):
    if significance_matrix_two_z_LL007[i,3] in lc_mask_LL007[:,3]:
        counter_two = counter_two + 1

significance_matrix_two_LC_z_LL007 = np.zeros((counter_two, 4))

c = 0
for i in range(np.shape(significance_matrix_two_z_LL007)[0]):
    if significance_matrix_two_z_LL007[i,3] in lc_mask_LL007[:,3]:
        significance_matrix_two_LC_z_LL007[c,:] = significance_matrix_two_z_LL007[i,:]
        c = c + 1






### Plotting

n_segments = 4
color = np.arange(-5.0,5.1,0.1)


fig, ax = plt.subplots(nrows=1, ncols=int(n_segments), figsize=(int(n_segments*1.7+2),4.8), constrained_layout=False)
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
#plt.xticks(fontsize=14)

plot0 = ax[0].contourf(xi_small,yi_small,zi_LL001,color, extend='both', cmap='jet')    # first image on screen
ax[0].plot(x_small,y_small,'k.')
ax[0].plot(lc_mask_LL001[:,0],lc_mask_LL001[:,1],'yo', markersize=8)
ax[0].plot(significance_matrix_z_LL001[:,0],significance_matrix_z_LL001[:,1],'b*', markersize=12)
ax[0].plot(significance_matrix_two_z_LL001[:,0],significance_matrix_two_z_LL001[:,1],'r*', markersize=12)
ax[0].plot(significance_matrix_LC_z_LL001[:,0],significance_matrix_LC_z_LL001[:,1],'b*', markersize=12, markeredgecolor='yellow')
ax[0].plot(significance_matrix_two_LC_z_LL001[:,0],significance_matrix_two_LC_z_LL001[:,1],'r*', markersize=12, markeredgecolor='yellow')
ax[0].set_ylim([-5, 280])
ax[0].invert_yaxis()
ax[0].set_xticklabels([])
ax[0].set_yticks([-0, 50, 100, 150, 200, 250]) 
ax[0].set_yticklabels([-0, 50, 100, 150, 200, 250], fontsize=13, fontweight='bold') 
#ax[0].set_title('Burst', fontsize=16, fontweight='bold')
ax[0].set_xticks([-25,0,25])
ax[0].set_xticklabels([-25, 0, 25], fontsize=13, fontweight='bold')
#ax[0].set_ylabel('Dorsoventral Axis (\u03BCm)',fontsize=16, fontweight='bold')
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_visible(False)
#ax[0].set_xlabel('Mediolateral Axis (\u03BCm)',fontsize=16, fontweight='bold')

plot1 = ax[1].contourf(xi_small,yi_small,zi_LL009,color, extend='both', cmap='jet')    # first image on screen
ax[1].plot(x_small,y_small,'k.')
ax[1].plot(lc_mask_LL009[:,0],lc_mask_LL009[:,1],'yo', markersize=8)
ax[1].plot(significance_matrix_z_LL009[:,0],significance_matrix_z_LL009[:,1],'b*', markersize=12)
ax[1].plot(significance_matrix_two_z_LL009[:,0],significance_matrix_two_z_LL009[:,1],'r*', markersize=12)
ax[1].plot(significance_matrix_LC_z_LL009[:,0],significance_matrix_LC_z_LL009[:,1],'b*', markersize=12, markeredgecolor='yellow')
ax[1].plot(significance_matrix_two_LC_z_LL009[:,0],significance_matrix_two_LC_z_LL009[:,1],'r*', markersize=12, markeredgecolor='yellow')
ax[1].set_ylim([-5, 280])
ax[1].invert_yaxis()
ax[1].set_yticks([])
ax[1].set_yticklabels([])
#ax[1].set_title('Inhibition', fontsize=16, fontweight='bold')
ax[1].set_xticks([])
ax[1].set_xticklabels([])  
#ax[1].set_ylabel('Dorsoventral Axi_smalls (\u03BCm)',fontsize=16, fontweight='bold')
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_visible(False)
#ax[1].set_xlabel('Mediolateral Axi_smalls (\u03BCm)',fontsize=16, fontweight='bold')

plot2 = ax[2].contourf(xi_small,yi_small,zi_LL017,color, extend='both', cmap='jet')    # first image on screen
ax[2].plot(x_small,y_small,'k.')
ax[2].plot(lc_mask_LL017[:,0],lc_mask_LL017[:,1],'yo', markersize=8)
ax[2].plot(significance_matrix_z_LL017[:,0],significance_matrix_z_LL017[:,1],'b*', markersize=12)
ax[2].plot(significance_matrix_two_z_LL017[:,0],significance_matrix_two_z_LL017[:,1],'r*', markersize=12)
ax[2].plot(significance_matrix_LC_z_LL017[:,0],significance_matrix_LC_z_LL017[:,1],'b*', markersize=12, markeredgecolor='yellow')
ax[2].plot(significance_matrix_two_LC_z_LL017[:,0],significance_matrix_two_LC_z_LL017[:,1],'r*', markersize=12, markeredgecolor='yellow')
ax[2].set_ylim([-5, 280])
ax[2].invert_yaxis()
ax[2].set_yticks([])
ax[2].set_yticklabels([])
#ax[2].set_title('Seizure 1', fontsize=16, fontweight='bold')
ax[2].set_xticks([])
ax[2].set_xticklabels([]) 
#ax[2].set_ylabel('Dorsoventral Axi_smalls (\u03BCm)',fontsize=16, fontweight='bold')
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[2].spines['left'].set_visible(False)
#ax[2].set_xlabel('Mediolateral Axi_smalls (\u03BCm)',fontsize=16, fontweight='bold')


plot3 = ax[3].contourf(xi_small,yi_small,zi_LL023,color, extend='both', cmap='jet')    # first image on screen
ax[3].plot(x_small,y_small,'k.')
ax[3].plot(lc_mask_LL023[:,0],lc_mask_LL023[:,1],'yo', markersize=8)
ax[3].plot(significance_matrix_z_LL023[:,0],significance_matrix_z_LL023[:,1],'b*', markersize=12)
ax[3].plot(significance_matrix_two_z_LL023[:,0],significance_matrix_two_z_LL023[:,1],'r*', markersize=12)
ax[3].plot(significance_matrix_LC_z_LL023[:,0],significance_matrix_LC_z_LL023[:,1],'b*', markersize=12, markeredgecolor='yellow')
ax[3].plot(significance_matrix_two_LC_z_LL023[:,0],significance_matrix_two_LC_z_LL023[:,1],'r*', markersize=12, markeredgecolor='yellow')
ax[3].set_ylim([-5, 280])
ax[3].invert_yaxis()
ax[3].set_yticks([])
ax[3].set_yticklabels([])
#ax[3].set_title('Seizure 2', fontsize=16, fontweight='bold')
ax[3].set_xticks([])
ax[3].set_xticklabels([])  
ax[3].spines['top'].set_visible(False)
ax[3].spines['right'].set_visible(False)
ax[3].spines['bottom'].set_visible(False)
ax[3].spines['left'].set_visible(False)



#fig.text(0.5, 0.04, 'Mediolateral Axis (\u03BCm)',fontsize=16, fontweight='bold' , ha='center')
fig.subplots_adjust(right = 0.80)

cbar_ax = fig.add_axes([0.85, 0.15, 0.035, 0.7])
cbar = fig.colorbar(plot0, cax=cbar_ax)
#cbar = fig.colorbar(plot0, ax=ax)
#cbar = fig.colorbar(plot1)
cbar.set_label('Spike Frequency (Z-score)', fontsize=14)
cbar.ax.tick_params(labelsize=12)
#plt.tight_layout()

savedir = savepath + r'\Fig5B_2.png'
plt.savefig(savedir, dpi=600)
plt.close()



### Plotting

n_segments = 1
color = np.arange(-5.0,5.1,0.1)


fig, ax = plt.subplots(nrows=1, ncols=int(n_segments), figsize=(int(n_segments*2.9),5), constrained_layout=True)
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
#plt.xticks(fontsize=14)

plot4 = ax.contourf(xi_big,yi_big,zi_LL007,color, extend='both', cmap='jet')    # first image on screen
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.plot(x_big,y_big,'k.')
ax.plot(lc_mask_LL007[:,0],lc_mask_LL007[:,1],'yo', markersize=8)
ax.plot(significance_matrix_z_LL007[:,0],significance_matrix_z_LL007[:,1],'b*', markersize=12)
ax.plot(significance_matrix_two_z_LL007[:,0],significance_matrix_two_z_LL007[:,1],'r*', markersize=12)
ax.plot(significance_matrix_LC_z_LL007[:,0],significance_matrix_LC_z_LL007[:,1],'b*', markersize=12, markeredgecolor='yellow')
ax.plot(significance_matrix_two_LC_z_LL007[:,0],significance_matrix_two_LC_z_LL007[:,1],'r*', markersize=12, markeredgecolor='yellow')
ax.set_ylim([-10, 560])
ax.set_xlim([-55, 55])
ax.invert_yaxis()
ax.set_yticks([-0, 100, 200, 300, 400, 500]) 
ax.set_yticklabels([-0, 100, 200, 300, 400, 500], fontsize=16, fontweight='bold') 
ax.set_xticks([-50,0,50])
ax.set_xticklabels([-50, 0, 50], fontsize=16, fontweight='bold') 

#ax[4].set_title('Seizure 3', fontsize=16, fontweight='bold')
#ax[4].set_xticks([-50,0,50])
#ax[4].set_xticklabels([-50, 0, 50], fontsize=14, fontweight='bold') 


#fig.text(0.5, 0.04, 'Mediolateral Axis (\u03BCm)',fontsize=16, fontweight='bold' , ha='center')
#fig.subplots_adjust(right = 0.80)

#cbar_ax = fig.add_axes([0.85, 0.15, 0.035, 0.7])
#cbar = fig.colorbar(plot0, cax=cbar_ax)
#cbar = fig.colorbar(plot0, ax=ax)
#cbar = fig.colorbar(plot1)
#cbar.ax.tick_params(labelsize=12)
#cbar.set_label('Spike Frequency (Z-score)', fontsize=16)
#plt.tight_layout()

savedir = savepath + r'\Fig5B_1.png'
plt.savefig(savedir, dpi=600)


