# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 16:21:25 2021

@author: llarsen
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, r'E:\Manuscript_analysis_files')

from vies.parse.phy import phy_data

savedir = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure01\output' # save directory
if not os.path.exists(savedir):
    os.makedirs(savedir)

eventfile = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure01\input\messages.events'
path_phy_data = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'

threshold = 2

srate = 30000

data = phy_data(path_phy_data, srate)

times = data.times

### Load pinch timestamps
f=open(eventfile, 'r')
events = f.read()
f.close()

events=events.split('pinch')
events=events[1:9]
events = [int(x) for x in events]
events = np.array(events)/srate
adjustments = np.array([-0.5, -0.4, -0.4, -0.4, -0.35, -0.3, -0.3, -0.55]) # alignment of trials
events = events - adjustments

###
times = data.times
templates = data.templates

good_clusters = [38]
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
        start = pinch - 21
        stop = pinch + 5
        
        index=[i for i in range(len(template_spikes)) if template_spikes[i] > start and template_spikes[i] < stop]
        spike_times_window = template_spikes[index] - start - 20
        if x == 0:
            spike_list[0]=spike_times_window
        else:
            spike_list.append(spike_times_window)  
        
    bins = 26
    window = 26
    binsize = window/bins
    
    bin_time = np.arange(0, window, binsize) + binsize/2 - 20
    bin_count = np.zeros(bins)
    
    for i in range(len(spike_list)):
        for k in range(len(bin_time)):
            for j in range(len(spike_list[i])):
                if spike_list[i][j] > bin_time[k] - binsize/2 and spike_list[i][j] <= bin_time[k] + binsize/2:
                    bin_count[k]=bin_count[k] + 1


    mean = np.mean(bin_count[0:18])
    std = np.std(bin_count[0:18])
    z_bins = (bin_count - mean) / std
    




plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(6,4))

ax.eventplot(spike_list, color='k', linelengths = 0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(False)
ax.set_xlim(-20,5)
ax.set_yticks([])
ax.set_xticks(np.arange(-20, 5.1, step=5))
ax.set_xlabel('Time (s)', fontsize=14)
ax.set_ylabel('Trials', fontsize=14)
plt.tight_layout()
save_filename = '\Pinch_response_example_a.png'
savepath = savedir + save_filename
plt.savefig(os.path.abspath(savepath), dpi=300)


fig, ax = plt.subplots(1, figsize=(6,4))
ax.bar(bin_time, z_bins, color='k')
ax.hlines(3, -20, 5, colors='r', linewidth=2, linestyles='dashed')
ax.hlines(-3, -20, 5, colors='r', linewidth=2, linestyles='dashed')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.set_xticks(np.arange(-20, 5.1, step=5))
ax.set_xlim(-20,5)
ax.set_ylim(-10,16)
ax.set_xlabel('Time (s)', fontsize=14)
ax.set_ylabel('Z-scored Spike Count', fontsize=14)
plt.tight_layout()
save_filename = '\Pinch_response_example_b.png'
savepath = savedir + save_filename
plt.savefig(os.path.abspath(savepath), dpi=300)