# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 16:21:25 2021

@author: llarsen
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, r'E:\OneDrive - UGent\python_functions')

from vies.parse.phy import phy_data

savedir = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure01\output' # directory to store figures
if not os.path.exists(savedir):
    os.makedirs(savedir)

lightfile = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure01\input\light017.npy'
path_phy_data = r'E:\Manuscript_analysis_files\data\phy_data\LL017-final.GUI'

### Load spike time stamps
threshold = 3 # significance threshold
srate = 30000
data = phy_data(path_phy_data, srate)
times = data.times
times = data.times
templates = data.templates
selected_neuron = [217] # the phy template used for plotting
spiketimes = times/srate

### Load pinch timestamps
events = np.load(lightfile)
events=events[6:21, 0]
events = np.delete(events, 12)

### loads and arranges spike time stamps per trial
response_matrix = np.zeros((len(selected_neuron), 3))
counter = 0
for m in selected_neuron:
    spike_list = [None]
    cluster_id = str(m)
    index = [i for i, x in enumerate(list(templates)) if x == m] # gets index of all matching values in templates
    template_spikes = spiketimes[index]
    
    for x in range(len(events)):
        pinch = events[x]
        start = pinch - 21
        stop = pinch + 14 + 8
        
        index=[i for i in range(len(template_spikes)) if template_spikes[i] > start and template_spikes[i] < stop]
        spike_times_window = template_spikes[index] - start

        if x == 0:
            spike_list[0]=spike_times_window
        else:
            spike_list.append(spike_times_window)
    
        
    bins = int(42/3)
    window = 42
    binsize = window/bins
    
    bin_time = np.arange(0, window, binsize) + binsize/2
    bin_count = np.zeros(bins)
    
    for i in range(len(spike_list)):
        for k in range(len(bin_time)):
            for j in range(len(spike_list[i])):
                if spike_list[i][j] > bin_time[k] - binsize/2 and spike_list[i][j] <= bin_time[k] + binsize/2:
                    bin_count[k]=bin_count[k] + 1
    
    
    mean = np.mean(bin_count[[0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13]])
    std = np.std(bin_count[[0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13]])
    z_bins = (bin_count - mean) / std
    

    
    
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(6,4))
ax.eventplot(spike_list, color='k', linelengths = 0.7)
ax.axvspan(21, 24, color='b', alpha=0.5, lw=0)
#ax.vlines(14, -0.5, 18.5, colors='r', linewidth=2, linestyles='dashed')
#ax.vlines(14+7, -0.5, 18.5, colors='r', linewidth=2, linestyles='dashed')
#ax.hlines(14.3, 21, 24, colors='b', linewidth=15, linestyles='solid')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(False)
ax.set_xlim(0,42)
ax.set_ylim(-0.6,13.6)
ax.set_yticks([])
#plt.suptitle('Response to Blue Light', fontsize=16, fontweight='bold')
ax.set_xticks(np.arange(0, 42.1, step=3))
ax.set_xticklabels(np.arange(-21, 21.1, step=3).astype(int))
ax.set_xlabel('Time (s)', fontsize=14)
ax.set_ylabel('Trials', fontsize=14)
save_filename = '\Light_response_example_a.png'
savepath = savedir + save_filename
plt.savefig(os.path.abspath(savepath), dpi=300)
plt.close()



fig, ax = plt.subplots(1, figsize=(6,4))
ax.bar(bin_time, z_bins, width=2.8, color='k')
ax.hlines(-3, 0, 42, colors='r', linewidth=2, linestyles='dashed')
ax.hlines(3, 0, 42, colors='r', linewidth=2, linestyles='dashed')
ax.axvspan(21, 24, color='b', alpha=0.5, lw=0)
#ax.vlines(3, -0.2, 9.2, colors='r', linewidth=2, linestyles='dashed')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.set_xlim(0,42)
ax.set_ylim(-15,5)
#plt.yticks([])
ax.set_xticks(np.arange(0, 42.1, step=3))
ax.set_xticklabels(np.arange(-21, 21.1, step=3).astype(int))
ax.set_xlabel('Time (s)', fontsize=14)
ax.set_ylabel('Z-scored Spike Count', fontsize=14)
#plt.rcParams["figure.figsize"] = (7,5)
save_filename = '\Light_response_example.png'
savepath = savedir + save_filename
plt.tight_layout()
#fig.subplots_adjust(top=0.92)
save_filename = '\Light_response_example_b.png'
savepath = savedir + save_filename
plt.savefig(os.path.abspath(savepath), dpi=300)
plt.close()