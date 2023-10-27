# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 12:15:00 2022

@author: llarsen
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.signal as signal
import h5py
from scipy.interpolate import griddata

sys.path.insert(1, r'E:\Manuscript_analysis_files')

from vies.spike.extract_frequency_bins import extract_frequency_bins

### LL001

## DATA FOR PINCH
MUA_filename = r'E:\Manuscript_analysis_files\data\mua_data\LL001.mua.hdf5'
probe_filename = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\small_poly3.prb'
savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\output'
messagefile = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\pinch_events\LL001.events'
LC_mask_file = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure05\input\LC_masks\LL001_masktest.npy'


#### LOAD SEIZURE DATA

trial_length = 120 # length of pinch trials in seconds
srate = 30000 # samples per second

bin_size = 5 # bin size in seconds
n_bins = int(360 / bin_size)
n_baseline_bins = int(300/bin_size)

srate_eeg=30000

sz1_start = 3782.72
sz2_start = 4570.89
sz3_start = 5191.40


Nt = int(trial_length/bin_size) # number of frames for animation

## STEP1: Load MUA data
with h5py.File(MUA_filename, "r") as f:
    n_channels = len(list(f['spiketimes'].keys()))
    labels = np.zeros(n_channels)
    ids = list(f['spiketimes'].keys())

    for i in range(n_channels):
        labels[i] = int(ids[i].split('elec_')[1])
    
    labels=np.array(labels)        
    labels=labels.astype(int)
  
    LL001_sz1 = np.zeros((n_bins, len(labels)))
    LL001_sz2 = np.zeros((n_bins, len(labels)))
    LL001_sz3 = np.zeros((n_bins, len(labels)))

    counter = 0
    for i in labels:
        index = int(i)
        e_id = ids[counter] 
        spiketimes = (f['spiketimes'][e_id][:]) / srate # recalculates spiketimes in seconds
        start_bin = sz1_start-300
        stop_bin = start_bin+360
        time, freq = extract_frequency_bins(spiketimes, bin_size, start_bin, stop_bin)
        LL001_sz1[:,index] = freq
        
        start_bin = sz2_start-300
        stop_bin = start_bin+360
        time, freq = extract_frequency_bins(spiketimes, bin_size, start_bin, stop_bin)
        LL001_sz2[:,index] = freq

        start_bin = sz3_start-300
        stop_bin = start_bin+360
        time, freq = extract_frequency_bins(spiketimes, bin_size, start_bin, stop_bin)
        LL001_sz3[:,index] = freq
        
        counter += 1

base_mean_freq = np.mean(LL001_sz1[0:n_baseline_bins,:,], axis=0)
base_std_freq = np.std(LL001_sz1[0:n_baseline_bins,:,], axis=0)
z_scored_matrix = np.zeros((np.shape(LL001_sz1)[0], np.shape(LL001_sz1)[1]))
for i in range(np.shape(LL001_sz1)[0]):
    z_scored_matrix[i,:] = (LL001_sz1[i,:] - base_mean_freq) / base_std_freq

z_LL001_sz1 = z_scored_matrix[(n_baseline_bins+2), :]

base_mean_freq = np.mean(LL001_sz2[0:n_baseline_bins,:,], axis=0)
base_std_freq = np.std(LL001_sz2[0:n_baseline_bins,:,], axis=0)
z_scored_matrix = np.zeros((np.shape(LL001_sz2)[0], np.shape(LL001_sz2)[1]))
for i in range(np.shape(LL001_sz2)[0]):
    z_scored_matrix[i,:] = (LL001_sz2[i,:] - base_mean_freq) / base_std_freq

z_LL001_sz2 = z_scored_matrix[(n_baseline_bins+2), :]

base_mean_freq = np.mean(LL001_sz3[0:n_baseline_bins,:,], axis=0)
base_std_freq = np.std(LL001_sz3[0:n_baseline_bins,:,], axis=0)
z_scored_matrix = np.zeros((np.shape(LL001_sz3)[0], np.shape(LL001_sz3)[1]))
for i in range(np.shape(LL001_sz3)[0]):
    z_scored_matrix[i,:] = (LL001_sz3[i,:] - base_mean_freq) / base_std_freq

z_LL001_sz3 = z_scored_matrix[(n_baseline_bins+2), :]


## STEP0: Load pinch timestamps, fill in path to time stamp file from open ephys below
f=open(messagefile, 'r')
events = f.read()
f.close()

n_pinches = 10 # number of pinch trials
trial_length = 25 # length of pinch trials in seconds
srate = 30000 # samples per second
bin_size = 1 # bin size in seconds

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

pinch = 19
z_burst = mua_matrix[pinch,:]
z_inhibition = mua_matrix[(pinch+1),:]


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


xi = np.arange(-25,26,5)
yi = np.arange(0, 276,5)
xi,yi = np.meshgrid(xi,yi)





threshold = 3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL001_sz1)):
    if z_LL001_sz1[i]>threshold:
        counter = counter + 1

significance_matrix_z_LL001_sz1 = np.zeros((counter, 4))
counter_two = 0
for i in range(len(z_LL001_sz1)):
    if z_LL001_sz1[i]>threshold:
        significance_matrix_z_LL001_sz1[counter_two,0] = x[i]
        significance_matrix_z_LL001_sz1[counter_two,1] = y[i]
        significance_matrix_z_LL001_sz1[counter_two,2] = z_LL001_sz1[i]        
        significance_matrix_z_LL001_sz1[counter_two,3] = i
        
        counter_two = counter_two + 1 
        
threshold = -3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL001_sz1)):
    if z_LL001_sz1[i]<threshold:
        counter = counter + 1

significance_matrix_two_z_LL001_sz1 = np.zeros((counter, 4))

counter_two = 0
for i in range(len(z_LL001_sz1)):
    if z_LL001_sz1[i]<threshold:
        significance_matrix_two_z_LL001_sz1[counter_two,0] = x[i]
        significance_matrix_two_z_LL001_sz1[counter_two,1] = y[i]
        significance_matrix_two_z_LL001_sz1[counter_two,2] = z_LL001_sz1[i]        
        significance_matrix_two_z_LL001_sz1[counter_two,3] = i  
        
        counter_two = counter_two + 1
        

threshold = 3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL001_sz2)):
    if z_LL001_sz2[i]>threshold:
        counter = counter + 1

significance_matrix_z_LL001_sz2 = np.zeros((counter, 4))
counter_two = 0
for i in range(len(z_LL001_sz2)):
    if z_LL001_sz2[i]>threshold:
        significance_matrix_z_LL001_sz2[counter_two,0] = x[i]
        significance_matrix_z_LL001_sz2[counter_two,1] = y[i]
        significance_matrix_z_LL001_sz2[counter_two,2] = z_LL001_sz2[i]        
        significance_matrix_z_LL001_sz2[counter_two,3] = i
        
        counter_two = counter_two + 1 
        
threshold = -3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL001_sz2)):
    if z_LL001_sz2[i]<threshold:
        counter = counter + 1

significance_matrix_two_z_LL001_sz2 = np.zeros((counter, 4))

counter_two = 0
for i in range(len(z_LL001_sz2)):
    if z_LL001_sz2[i]<threshold:
        significance_matrix_two_z_LL001_sz2[counter_two,0] = x[i]
        significance_matrix_two_z_LL001_sz2[counter_two,1] = y[i]
        significance_matrix_two_z_LL001_sz2[counter_two,2] = z_LL001_sz2[i]        
        significance_matrix_two_z_LL001_sz2[counter_two,3] = i  
        
        counter_two = counter_two + 1 
        

threshold = 3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL001_sz3)):
    if z_LL001_sz3[i]>threshold:
        counter = counter + 1

significance_matrix_z_LL001_sz3 = np.zeros((counter, 4))
counter_two = 0
for i in range(len(z_LL001_sz3)):
    if z_LL001_sz3[i]>threshold:
        significance_matrix_z_LL001_sz3[counter_two,0] = x[i]
        significance_matrix_z_LL001_sz3[counter_two,1] = y[i]
        significance_matrix_z_LL001_sz3[counter_two,2] = z_LL001_sz3[i]        
        significance_matrix_z_LL001_sz3[counter_two,3] = i
        
        counter_two = counter_two + 1 
        
threshold = -3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_LL001_sz3)):
    if z_LL001_sz3[i]<threshold:
        counter = counter + 1

significance_matrix_two_z_LL001_sz3 = np.zeros((counter, 4))

counter_two = 0
for i in range(len(z_LL001_sz3)):
    if z_LL001_sz3[i]<threshold:
        significance_matrix_two_z_LL001_sz3[counter_two,0] = x[i]
        significance_matrix_two_z_LL001_sz3[counter_two,1] = y[i]
        significance_matrix_two_z_LL001_sz3[counter_two,2] = z_LL001_sz3[i]        
        significance_matrix_two_z_LL001_sz3[counter_two,3] = i  
        
        counter_two = counter_two + 1 



threshold = 3.16 # correction for multiple testing, i.e. 32 channels, 0.05/32 = 0.0015625, which equates to a threshold z score of 3.16 for two tailed statistics
counter = 0
for i in range(len(z_burst)):
    if z_burst[i]>threshold:
        counter = counter + 1

significance_matrix_burst = np.zeros((counter, 4))

counter_two = 0
for i in range(len(z_burst)):
    if z_burst[i]>threshold:
        significance_matrix_burst[counter_two,0] = x[i]
        significance_matrix_burst[counter_two,1] = y[i]
        significance_matrix_burst[counter_two,2] = z_burst[i]        
        significance_matrix_burst[counter_two,3] = i
        
        counter_two = counter_two + 1 
        
counter = 0
for i in range(len(z_inhibition)):
    if  z_inhibition[i]<-threshold:
        counter = counter + 1

significance_matrix_inhibition = np.zeros((counter, 4))

counter_two = 0
for i in range(len(z_burst)):
    if  z_inhibition[i]<-threshold:
        significance_matrix_inhibition[counter_two,0] = x[i]
        significance_matrix_inhibition[counter_two,1] = y[i]
        significance_matrix_inhibition[counter_two,2] = z_inhibition[i]        
        significance_matrix_inhibition[counter_two,3] = i
        
        counter_two = counter_two + 1 

lc_mask = np.load(LC_mask_file)

### sz1
counter = 0
for i in range(np.shape(significance_matrix_z_LL001_sz1)[0]):
    if significance_matrix_z_LL001_sz1[i,3] in lc_mask[:,3]:
        counter = counter + 1

significance_matrix_LC_z_LL001_sz1 = np.zeros((counter, 4))

c = 0
for i in range(np.shape(significance_matrix_z_LL001_sz1)[0]):
    if significance_matrix_z_LL001_sz1[i,3] in lc_mask[:,3]:
        significance_matrix_LC_z_LL001_sz1[c,:] = significance_matrix_z_LL001_sz1[i,:]
        c = c + 1
        
counter_two = 0
for i in range(np.shape(significance_matrix_two_z_LL001_sz1)[0]):
    if significance_matrix_two_z_LL001_sz1[i,3] in lc_mask[:,3]:
        counter_two = counter_two + 1

significance_matrix_two_LC_z_LL001_sz1 = np.zeros((counter_two, 4))

c = 0
for i in range(np.shape(significance_matrix_two_z_LL001_sz1)[0]):
    if significance_matrix_two_z_LL001_sz1[i,3] in lc_mask[:,3]:
        significance_matrix_two_LC_z_LL001_sz1[c,:] = significance_matrix_two_z_LL001_sz1[i,:]
        c = c + 1


### sz2
counter = 0
for i in range(np.shape(significance_matrix_z_LL001_sz2)[0]):
    if significance_matrix_z_LL001_sz2[i,3] in lc_mask[:,3]:
        counter = counter + 1

significance_matrix_LC_z_LL001_sz2 = np.zeros((counter, 4))

c = 0
for i in range(np.shape(significance_matrix_z_LL001_sz2)[0]):
    if significance_matrix_z_LL001_sz2[i,3] in lc_mask[:,3]:
        significance_matrix_LC_z_LL001_sz2[c,:] = significance_matrix_z_LL001_sz2[i,:]
        c = c + 1
        
counter_two = 0
for i in range(np.shape(significance_matrix_two_z_LL001_sz2)[0]):
    if significance_matrix_two_z_LL001_sz2[i,3] in lc_mask[:,3]:
        counter_two = counter_two + 1

significance_matrix_two_LC_z_LL001_sz2 = np.zeros((counter_two, 4))

c = 0
for i in range(np.shape(significance_matrix_two_z_LL001_sz2)[0]):
    if significance_matrix_two_z_LL001_sz2[i,3] in lc_mask[:,3]:
        significance_matrix_two_LC_z_LL001_sz2[c,:] = significance_matrix_two_z_LL001_sz2[i,:]
        c = c + 1
        

### sz3
counter = 0
for i in range(np.shape(significance_matrix_z_LL001_sz3)[0]):
    if significance_matrix_z_LL001_sz3[i,3] in lc_mask[:,3]:
        counter = counter + 1

significance_matrix_LC_z_LL001_sz3 = np.zeros((counter, 4))

c = 0
for i in range(np.shape(significance_matrix_z_LL001_sz3)[0]):
    if significance_matrix_z_LL001_sz3[i,3] in lc_mask[:,3]:
        significance_matrix_LC_z_LL001_sz3[c,:] = significance_matrix_z_LL001_sz3[i,:]
        c = c + 1
        
counter_two = 0
for i in range(np.shape(significance_matrix_two_z_LL001_sz3)[0]):
    if significance_matrix_two_z_LL001_sz3[i,3] in lc_mask[:,3]:
        counter_two = counter_two + 1

significance_matrix_two_LC_z_LL001_sz3 = np.zeros((counter_two, 4))

c = 0
for i in range(np.shape(significance_matrix_two_z_LL001_sz3)[0]):
    if significance_matrix_two_z_LL001_sz3[i,3] in lc_mask[:,3]:
        significance_matrix_two_LC_z_LL001_sz3[c,:] = significance_matrix_two_z_LL001_sz3[i,:]
        c = c + 1


ziseizure1 = griddata((x,y),z_LL001_sz1,(xi,yi),method='cubic')
ziseizure2 = griddata((x,y),z_LL001_sz2,(xi,yi),method='cubic')
ziseizure3 = griddata((x,y),z_LL001_sz3,(xi,yi),method='cubic')

zi_burst = griddata((x,y),z_burst,(xi,yi),method='cubic')
zi_inhibition = griddata((x,y),z_inhibition,(xi,yi),method='cubic')

color = np.arange(-5.0,5.1,0.1)

n_segments = 5

fig, ax = plt.subplots(nrows=1, ncols=int(n_segments), figsize=(int(n_segments*1.7+2),5.2))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
#plt.xticks(fontsize=14)

plot0 = ax[0].contourf(xi,yi,zi_burst,color, extend='both', cmap='jet')    # first image on screen
ax[0].plot(x,y,'k.')
ax[0].plot(significance_matrix_burst[:,0],significance_matrix_burst[:,1],'y*', markersize=12)
ax[0].set_ylim([-5, 280])
ax[0].invert_yaxis()
ax[0].set_xticklabels([])
ax[0].set_yticks([-0, 50, 100, 150, 200, 250]) 
ax[0].set_yticklabels([-0, 50, 100, 150, 200, 250], fontsize=12, fontweight='bold') 
ax[0].set_title('Burst', fontsize=14, fontweight='bold')
ax[0].set_xticks([])
ax[0].set_xticklabels([], fontsize=14, fontweight='bold') 
#ax[0].set_ylabel('Dorsoventral Axis (\u03BCm)',fontsize=14, fontweight='bold')
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_visible(False)
#ax[0].set_xlabel('Mediolateral Axis (\u03BCm)',fontsize=16, fontweight='bold')

plot0 = ax[1].contourf(xi,yi,zi_inhibition,color, extend='both', cmap='jet')    # first image on screen
ax[1].plot(x,y,'k.')
ax[1].plot(significance_matrix_inhibition[:,0],significance_matrix_inhibition[:,1],'y*', markersize=12)
ax[1].set_ylim([-5, 280])
ax[1].invert_yaxis()
ax[1].set_xticklabels([])
ax[1].set_yticks([])
ax[1].set_yticklabels([])
ax[1].set_title('Inhibition', fontsize=14, fontweight='bold')
ax[1].set_xticks([])
ax[1].set_xticklabels([], fontsize=14, fontweight='bold') 
#ax[1].set_ylabel('Dorsoventral Axis (\u03BCm)',fontsize=16, fontweight='bold')
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_visible(False)
#ax[1].set_xlabel('Mediolateral Axis (\u03BCm)',fontsize=16, fontweight='bold')

plot0 = ax[2].contourf(xi,yi,ziseizure1,color, extend='both', cmap='jet')    # first image on screen
ax[2].plot(x,y,'k.')
ax[2].plot(lc_mask[:,0],lc_mask[:,1],'yo', markersize=8)
ax[2].plot(significance_matrix_LC_z_LL001_sz1[:,0],significance_matrix_LC_z_LL001_sz1[:,1],'b*', markersize=12, markeredgecolor='yellow')
ax[2].plot(significance_matrix_two_LC_z_LL001_sz1[:,0],significance_matrix_two_LC_z_LL001_sz1[:,1],'r*', markersize=12, markeredgecolor='yellow')
ax[2].set_ylim([-5, 280])
ax[2].invert_yaxis()
ax[2].set_yticks([])
ax[2].set_yticklabels([])
ax[2].set_title('Seizure 1', fontsize=14, fontweight='bold')
ax[2].set_xticks([-25,0,25])
ax[2].set_xticklabels([-25, 0, 25], fontsize=12, fontweight='bold') 
#ax[2].set_ylabel('Dorsoventral Axis (\u03BCm)',fontsize=16, fontweight='bold')
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[2].spines['left'].set_visible(False)
#ax[2].set_xlabel('Mediolateral Axis (\u03BCm)',fontsize=16, fontweight='bold')


plot1 = ax[3].contourf(xi,yi,ziseizure2,color, extend='both', cmap='jet')    # first image on screen
ax[3].plot(x,y,'k.')
ax[3].plot(lc_mask[:,0],lc_mask[:,1],'yo', markersize=8)
ax[3].plot(significance_matrix_LC_z_LL001_sz2[:,0],significance_matrix_LC_z_LL001_sz2[:,1],'b*', markersize=12, markeredgecolor='yellow')
ax[3].plot(significance_matrix_two_LC_z_LL001_sz2[:,0],significance_matrix_two_LC_z_LL001_sz2[:,1],'r*', markersize=12, markeredgecolor='yellow')
ax[3].set_ylim([-5, 280])
ax[3].invert_yaxis()
ax[3].set_yticks([])
ax[3].set_yticklabels([])
ax[3].set_title('Seizure 2', fontsize=14, fontweight='bold')
ax[3].set_xticks([])
ax[3].set_xticklabels([]) 
ax[3].spines['top'].set_visible(False)
ax[3].spines['right'].set_visible(False)
ax[3].spines['bottom'].set_visible(False)
ax[3].spines['left'].set_visible(False)

plot2 = ax[4].contourf(xi,yi,ziseizure3,color, extend='both', cmap='jet')    # first image on screen
ax[4].spines['top'].set_visible(False)
ax[4].spines['right'].set_visible(False)
ax[4].spines['bottom'].set_visible(False)
ax[4].spines['left'].set_visible(False)
ax[4].plot(x,y,'k.')
ax[4].plot(lc_mask[:,0],lc_mask[:,1],'yo', markersize=8)
ax[4].plot(significance_matrix_LC_z_LL001_sz3[:,0],significance_matrix_LC_z_LL001_sz3[:,1],'b*', markersize=12, markeredgecolor='yellow')
ax[4].plot(significance_matrix_two_LC_z_LL001_sz3[:,0],significance_matrix_two_LC_z_LL001_sz3[:,1],'r*', markersize=12, markeredgecolor='yellow')
ax[4].set_ylim([-5, 280])
ax[4].invert_yaxis()
ax[4].set_yticks([])
ax[4].set_yticklabels([])
ax[4].set_title('Seizure 3', fontsize=14, fontweight='bold')
ax[4].set_xticks([])
ax[4].set_xticklabels([]) 

#fig.text(0.5, 0.04, 'Mediolateral Axis (\u03BCm)',fontsize=14, fontweight='bold' , ha='center')
fig.subplots_adjust(right=0.84)

cbar_ax = fig.add_axes([0.87, 0.13, 0.025, 0.7])
cbar = fig.colorbar(plot2, cax=cbar_ax)
#cbar = fig.colorbar(plot1)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Spike Frequency (Z-score)', fontsize=14)
#plt.tight_layout()

savedir = savepath + r'\Fig5A.png'
plt.savefig(savedir, dpi=600)

