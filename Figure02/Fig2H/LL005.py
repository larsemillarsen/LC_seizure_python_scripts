# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 12:42:02 2021

@author: llarsen
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, r'E:\Manuscript_analysis_files')

from vies.parse.phy import phy_data
import scipy
import pandas as pd

eventfile = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\input\pinch_events\LL005.events'
path_phy_data = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'
lightfile = 'NONE'
savedir = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\Fig2H\output_data\LL005.npy'
#if not os.path.exists(savedir):
#    os.makedirs(savedir)
       
p_threshold = 0.05

srate = 30000
data = phy_data(path_phy_data, srate)
times = data.times
good_clusters = [16, 17, 20, 26, 36, 38, 40, 47, 53, 61, 63, 71, 74, 77, 79, 82, 84, 95, 100, 102, 112, 120, 124, 128, 136, 146, 149]

sz1_start = 770.52
sz2_start = 1664.58
sz3_start = 2309.58 

light_trials = 0 # IF = 0, no light trials. IF = 1, light trials were performed

###############################
seg1_start = sz1_start - 300
seg1_stop = sz1_start + 60
seg2_start = sz2_start - 300
seg2_stop = sz2_start + 60
seg3_start = sz3_start - 300
seg3_stop = sz3_start + 60

seizure_times = [seg1_start, seg1_stop, seg2_start, seg2_stop, seg3_start, seg3_stop]

### evaluate neuron firing characteristics
bin_size = 5
spike_count_sz1 = np.zeros((int(360/bin_size), len(good_clusters)))
spike_count_sz2 = np.zeros((int(360/bin_size), len(good_clusters)))
spike_count_sz3 = np.zeros((int(360/bin_size), len(good_clusters)))
counter = 0
for t in good_clusters:
    spike_freq, spike_count_sz1[:,counter], template_id = data.extract_frequency_bins(bin_size, seg1_start, seg1_stop, templates = t)
    spike_freq, spike_count_sz2[:,counter], template_id = data.extract_frequency_bins(bin_size, seg2_start, seg2_stop, templates = t)
    spike_freq, spike_count_sz3[:,counter], template_id = data.extract_frequency_bins(bin_size, seg3_start, seg3_stop, templates = t)
    counter = counter + 1

min_firing = 0
firing_criteria = np.zeros((len(good_clusters)))
for i in range(len(firing_criteria)):
    if np.min(spike_count_sz1[0:59,i]) > min_firing and np.min(spike_count_sz3[0:59,i]) > min_firing and np.min(spike_count_sz3[0:59,i]) > min_firing:
        firing_criteria[i] = 1
        

### pinch response evaluation
f=open(eventfile, 'r')
events = f.read()
f.close()

events=events.split('pinch')
events=events[1:9]
events = [int(x) for x in events]
events = np.array(events)/srate
adjustments = np.array([-0.5, -0.4, -0.4, -0.4, -0.35, -0.3, -0.3, -0.55])
events = events - adjustments

templates = data.templates

spiketimes = times/srate
response_matrix = np.zeros((len(good_clusters), 4))

counter = 0
for m in good_clusters:
    spike_list = [None]
    cluster_id = str(m)
    index = [i for i, x in enumerate(list(templates)) if x == m] # gets index of all matching values in templates
    template_spikes = spiketimes[index]

    for x in range(len(events)):
        pinch = events[x]
        start = pinch - 20
        stop = pinch + 5
        
        index=[i for i in range(len(template_spikes)) if template_spikes[i] > start and template_spikes[i] < stop]
        spike_times_window = template_spikes[index] - start

        if x == 0:
            spike_list[0]=spike_times_window
        else:
            spike_list.append(spike_times_window)
        
    bins = 25
    window = 25
    binsize = window/bins
    
    bin_time = np.arange(0, window, binsize) + binsize/2
    bin_count = np.zeros(bins)
    
    for i in range(len(spike_list)):
        for k in range(len(bin_time)):
            for j in range(len(spike_list[i])):
                if spike_list[i][j] > bin_time[k] - binsize/2 and spike_list[i][j] <= bin_time[k] + binsize/2:
                    bin_count[k]=bin_count[k] + 1


    mean = np.mean(bin_count[0:18])
    std = np.std(bin_count[0:18])
    z_bins = (bin_count - mean) / std
    
    response_matrix[counter, 0] =  z_bins[19]
    response_matrix[counter, 1] =  z_bins[20]
    response_matrix[counter, 2] =  z_bins[19] - z_bins[20]
    response_matrix[counter, 3] = m 
    
    counter = counter + 1


threshold = 3
LC_neurons = 'Based on pinch response, templates '
no_neurons = 0
for i in range(len(good_clusters)):
    if response_matrix[i,0] > threshold and response_matrix[i,1] < -threshold and i<len(good_clusters)-1:
        LC_neurons = LC_neurons + str(int(response_matrix[i,3])) + ', '
        no_neurons = no_neurons + 1
    elif response_matrix[i,0] > threshold and response_matrix[i,1] < -threshold and i==len(good_clusters)-1:
        LC_neurons = LC_neurons + str(int(response_matrix[i,3]))
        no_neurons = no_neurons + 1

pinch_response = np.zeros((len(good_clusters)))
for i in range(len(good_clusters)):
    if response_matrix[i,0] > threshold and response_matrix[i,1] < -threshold:
        pinch_response[i]=1
        

### Load pinch timestamps
light_response = np.zeros((len(good_clusters)))
if light_trials == 1:
    events = np.load(lightfile)
    events=events[0:19, 0]
    
    response_matrix = np.zeros((len(good_clusters), 3))
    counter = 0
    for m in good_clusters:
        spike_list = [None]
        cluster_id = str(m)
        index = [i for i, x in enumerate(list(templates)) if x == m] # gets index of all matching values in templates
        template_spikes = spiketimes[index]
        
        for x in range(len(events)):
            pinch = events[x]
            start = pinch - 14
            stop = pinch + 14 + 7
            
            index=[i for i in range(len(template_spikes)) if template_spikes[i] > start and template_spikes[i] < stop]
            spike_times_window = template_spikes[index] - start
    
            if x == 0:
                spike_list[0]=spike_times_window
            else:
                spike_list.append(spike_times_window)    
            
        bins = int(35/7)
        window = 35
        binsize = window/bins
        
        bin_time = np.arange(0, window, binsize) + binsize/2
        bin_count = np.zeros(bins)
        
        for i in range(len(spike_list)):
            for k in range(len(bin_time)):
                for j in range(len(spike_list[i])):
                    if spike_list[i][j] > bin_time[k] - binsize/2 and spike_list[i][j] <= bin_time[k] + binsize/2:
                        bin_count[k]=bin_count[k] + 1 
        
        mean = np.mean(bin_count[[0, 1, 3, 4]])
        std = np.std(bin_count[[0, 1, 3, 4]])
        z_bins = (bin_count - mean) / std
        
        response_matrix[counter, 0] =  z_bins[2]
        response_matrix[counter, 1] = m 
        
        counter = counter + 1
    
    threshold = 3
    LC_neurons = 'Based on light response, templates '
    no_neurons = 0
    for i in range(len(good_clusters)):
        if response_matrix[i,0] < -threshold and i<len(good_clusters)-1:
            LC_neurons = LC_neurons + str(int(response_matrix[i,1])) + ', '
            no_neurons = no_neurons + 1
        elif response_matrix[i,0] < -threshold  and i==len(good_clusters)-1:
            LC_neurons = LC_neurons + str(int(response_matrix[i,1]))
            no_neurons = no_neurons + 1
    
    for i in range(len(good_clusters)):
        if response_matrix[i,0] < -threshold:
            light_response[i]=1    

response_matrix = np.zeros((len(good_clusters), 4))
response_matrix[:,0] = good_clusters
response_matrix[:,1] = firing_criteria
response_matrix[:,2] = pinch_response
response_matrix[:,3] = light_response

tagged_by = []

LC_neurons = 'Based on firing stability, a positive pinch or light response, templates '
no_neurons=0
lc_neuron_list = []
for i in range(len(good_clusters)):
    if response_matrix[i,2] == 1 or response_matrix[i,3] == 1:
        if response_matrix[i,1] == 1:
            LC_neurons = LC_neurons + str(int(response_matrix[i,0])) + ', '
            lc_neuron_list.append(int(response_matrix[i,0]))
            no_neurons=no_neurons+1
            if response_matrix[i,2] == 1 and response_matrix[i,3] == 1:
                tagged_by.append('Pinch and Light')
            elif response_matrix[i,2] == 1 or response_matrix[i,3] == 0:
                tagged_by.append('Pinch')
            elif response_matrix[i,2] == 0 or response_matrix[i,3] == 1:
                tagged_by.append('Light')
                
LC_neurons = LC_neurons + ' are estimated to be stabily recorded LC neurons. A total of ' + str(no_neurons) + ' LC neurons were found.'
print(LC_neurons)

#####

counter = 0
z_matrix = np.zeros((len(lc_neuron_list), 3))
effect = np.zeros((len(lc_neuron_list), 3))
change = np.zeros((3,3))
for t in lc_neuron_list:
    spike_freq, spike_count_sz1, template_id = data.extract_frequency_bins(bin_size, seg1_start, seg1_stop, templates = t)
    spike_freq, spike_count_sz2, template_id = data.extract_frequency_bins(bin_size, seg2_start, seg2_stop, templates = t)
    spike_freq, spike_count_sz3, template_id = data.extract_frequency_bins(bin_size, seg3_start, seg3_stop, templates = t)
    
    spike_count = np.column_stack((spike_count_sz1, spike_count_sz2, spike_count_sz3))
 
    basemean_spike_count_sz1 = np.mean(spike_count_sz1[0:59])
    basemean_std_sz1 = np.std(spike_count_sz1[0:59])
    spike_count_sz1_norm = spike_count_sz1 / basemean_spike_count_sz1 
    spike_count_sz1_z = (spike_count_sz1 - basemean_spike_count_sz1) / basemean_std_sz1 
    spike_count_sz1_norm = spike_count_sz1_norm[62]
    z_seizure1 = spike_count_sz1_z[62]

    basemean_spike_count_sz2 = np.mean(spike_count_sz2[0:59])
    basemean_std_sz2 = np.std(spike_count_sz2[0:59])
    spike_count_sz2_norm = spike_count_sz2 / basemean_spike_count_sz2 
    spike_count_sz2_z = (spike_count_sz2 - basemean_spike_count_sz2) / basemean_std_sz2 
    spike_count_sz2_norm = spike_count_sz2_norm[62]
    z_seizure2 = spike_count_sz2_z[62]
    
    basemean_spike_count_sz3 = np.mean(spike_count_sz3[0:59])
    basemean_std_sz3 = np.std(spike_count_sz3[0:59])
    spike_count_sz3_norm = spike_count_sz3 / basemean_spike_count_sz3 
    spike_count_sz3_z = (spike_count_sz3 - basemean_spike_count_sz3) / basemean_std_sz3 
    spike_count_sz3_norm = spike_count_sz3_norm[62]
    z_seizure3 = spike_count_sz3_z[62]
    
    z_matrix[counter,0]=z_seizure1
    z_matrix[counter,1]=z_seizure2
    z_matrix[counter,2]=z_seizure3

    effect[counter,0]=spike_count_sz1_norm * 100
    effect[counter,1]=spike_count_sz2_norm * 100
    effect[counter,2]=spike_count_sz3_norm * 100
    
    counter = counter + 1
    
p_values = scipy.stats.norm.sf(abs(z_matrix))*2

for y in range(np.shape(p_values)[1]):
    counter_inhibited = 0
    counter_excited = 0
    counter_nochange = 0
    for x in range(np.shape(p_values)[0]):
        if z_matrix[x,y] < 0 and p_values[x,y] < p_threshold:
            counter_inhibited = counter_inhibited + 1
        elif z_matrix[x,y] > 0 and p_values[x,y] < p_threshold:
            counter_excited = counter_excited + 1
        else:
            counter_nochange = counter_nochange + 1

    change[0,y] = counter_inhibited
    change[1,y] = counter_excited
    change[2,y] = counter_nochange        

   
headers = ['template ID', 'effect sz1', 'effect sz2', 'effect sz3', 'z sz1', 'z sz2','z sz3', 'p sz1', 'p sz2', 'p sz3', 'tagged by',]
a_excel_matrix = np.column_stack((lc_neuron_list, effect, z_matrix, p_values, tagged_by))

np.save(savedir, a_excel_matrix)

a_excel_matrix = pd.DataFrame(a_excel_matrix, columns=headers)

