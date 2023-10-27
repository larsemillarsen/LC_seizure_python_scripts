# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 12:34:08 2022

@author: User
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.signal import find_peaks

sys.path.insert(1, r'E:\Manuscript_analysis_files')

from vies.parse.neuron import load_neuronfile
from vies.lfp.filter import butter_bandpass_filter
from vies.parse.phy import phy_data


### LL007
eeg_path = r'E:\Data\LCSeizureData\LCSeizureNIDAQ\16_LL007\Seizure\Seizure0.mat'
eeg_srate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'
LL007_spike_templates = [81]
seizure_start = 1193.42


filter_low = 200
filter_high = 4500
unit_srate = 30000
eeg_srate = eeg_srate
templates = LL007_spike_templates

time_eeg, hip_eeg = load_neuronfile(eeg_path, srate=eeg_srate, channel=1, gain=200, inputrange=20)
filtered_lfp = butter_bandpass_filter(hip_eeg, filter_low, filter_high, eeg_srate, order=1)

sd_signal = np.std(filtered_lfp)
peaks = find_peaks(-filtered_lfp, height=5*sd_signal, distance=0.01*eeg_srate)
peaks_seconds = peaks[0] / eeg_srate
peaks_height = peaks[1]['peak_heights']

peaks_seconds = peaks_seconds - eeg_stimstart  
path_phy_data = phy_path

data = phy_data(path_phy_data, unit_srate)

spike_times = np.array(data.extract_spike_trains(LL007_spike_templates)) - seizure_start

spike_times = spike_times[spike_times>-120]
spike_times = spike_times[spike_times<120]


plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(8,2))
ax.plot(time_eeg-120, filtered_lfp*2 + 5, 'b', linewidth = 0.6, alpha = 1)
ax.scatter(peaks_seconds, -peaks_height*2+5, marker='x', color='r')
ax.eventplot(spike_times, color='g', linelengths = 4, linewidth=1, lineoffsets = - 5)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticks([])
ax.set_yticklabels([]) 
ax.set_xticks([])
ax.set_xticklabels([])

ax.set_xlim((10, 15))
savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure07\output\Fig7E_2.png'
plt.savefig(os.path.abspath(savepath), dpi=600) 












