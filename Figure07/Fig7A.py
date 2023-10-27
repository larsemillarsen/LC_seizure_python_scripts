# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 11:35:53 2022

@author: llarsen
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, r'E:\OneDrive - UGent\python_functions')

from vies.parse.neuron import load_neuronfile
from vies.lfp.filter import butter_bandpass_filter
from vies.parse.phy import phy_data
from vies.spike.extract_frequency_bins import extract_frequency_bins
from vies.spike_lfp.lfp_phase_spike_coherence import spike_lfp_radians
from astropy.stats.circstats import rayleightest


### LL023
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL023\4_Seizure300_2201218A0004.mat'
eeg_srate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL023-final.GUI'

seizure_start = 2061.56
time_eeg, hip_eeg = load_neuronfile(eeg_path, srate=eeg_srate, channel=1, gain=200, inputrange=20)

filter_low = 8
filter_high = 16

start = int((eeg_stimstart + 10) * eeg_srate)
stop = int((eeg_stimstart + 20) * eeg_srate)
start = int((eeg_stimstart - 10) * eeg_srate)
stop = int((eeg_stimstart + 50) * eeg_srate)


### FULL SEIZURE
fig, ax = plt.subplots(1, figsize=(8,2))
ax.plot(time_eeg[start:stop], hip_eeg[start:stop], 'k', linewidth = 0.15)
ax.vlines(60, np.min(hip_eeg[start:stop]), np.max(hip_eeg[start:stop]), colors='r', linewidth=2, linestyles='dashed', zorder=3)
ax.vlines(70, np.min(hip_eeg[start:stop]), np.max(hip_eeg[start:stop]), colors='r', linewidth=2, linestyles='dashed', zorder=3)
ax.axvspan(77, 80, color='b', alpha=0.6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticks([])
ax.set_yticklabels([]) 
ax.set_xticks([])
ax.set_xticklabels([]) 

plt.tight_layout()
savepath2 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure07\output\Fig7A_LFP.png'
#plt.savefig(os.path.abspath(savepath2), dpi=300)



#### FILTERED EXAMPLES
filtered_lfp = butter_bandpass_filter(hip_eeg, filter_low, filter_high, eeg_srate, order=1)


### LOAD spike data
path_phy_data = phy_path
unit_srate = 30000
data = phy_data(path_phy_data, unit_srate)

spikes_coupled = data.extract_spike_trains(48) - seizure_start
spikes_uncoupled = data.extract_spike_trains(33) - seizure_start

period_on = 17
period_off = 20
start = int((eeg_stimstart + period_on) * eeg_srate)
stop = int((eeg_stimstart + period_off) * eeg_srate)

test = extract_frequency_bins(spikes_coupled, 1, -20, 60)[1]

spikes_coupled = spikes_coupled[spikes_coupled>period_on]
spikes_coupled = spikes_coupled[spikes_coupled<period_off]

spikes_uncoupled = spikes_uncoupled[spikes_uncoupled>period_on]
spikes_uncoupled = spikes_uncoupled[spikes_uncoupled<period_off]


fig, ax = plt.subplots(1, figsize=(8,4))
ax.plot(time_eeg[start:stop], hip_eeg[start:stop], 'k', linewidth = 1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticks([])
ax.set_yticklabels([]) 
ax.set_xticks([])
ax.set_xticklabels([]) 
ax.plot(time_eeg[start:stop], filtered_lfp[start:stop] * 4, 'b', linewidth = 2, linestyle='--', alpha = 0.5)

savepath2 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure07\output\Fig7A_LFP_zoom.png'
plt.savefig(os.path.abspath(savepath2), dpi=600)




fig, ax = plt.subplots(2, figsize=(8,4))
ax[0].plot(time_eeg[start:stop]-60, filtered_lfp[start:stop] + 1, 'b', linewidth = 1, linestyle='--', alpha = 0.5)
ax[0].eventplot(spikes_coupled, color='g', linelengths = 0.7)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_visible(False)
#ax[0].set_xlim(period_on,period_off)
ax[0].set_yticks([])
ax[0].set_yticklabels([]) 
ax[0].set_xticks([])
ax[0].set_xticklabels([])


ax[1].plot(time_eeg[start:stop]-60, filtered_lfp[start:stop] + 1, 'b', linewidth = 1, linestyle='--', alpha = 0.5)
ax[1].eventplot(spikes_uncoupled, color='r', linelengths = 0.7)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_visible(False)
#ax[1].set_xlim(period_on,period_off)
ax[1].set_yticks([])
ax[1].set_yticklabels([]) 
ax[1].set_xticks([])
ax[1].set_xticklabels([])


plt.tight_layout()
savepath2 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure07\output\Fig7A_spikes.png'
plt.savefig(os.path.abspath(savepath2), dpi=600)





spikes_coupled = data.extract_spike_trains(48) - seizure_start - 10
spikes_uncoupled = data.extract_spike_trains(33) - seizure_start - 10

period_on = 0
period_off = 10
start = int((eeg_stimstart + period_on) * eeg_srate)
stop = int((eeg_stimstart + period_off) * eeg_srate)

test = extract_frequency_bins(spikes_coupled, 1, -20, 60)[1]

spikes_coupled = spikes_coupled[spikes_coupled>period_on]
spikes_coupled = spikes_coupled[spikes_coupled<period_off]

spikes_uncoupled = spikes_uncoupled[spikes_uncoupled>period_on]
spikes_uncoupled = spikes_uncoupled[spikes_uncoupled<period_off]

start = int((eeg_stimstart + 10) * eeg_srate)
stop = int((eeg_stimstart + 20) * eeg_srate)

spike_radians_coupled = spike_lfp_radians(filtered_lfp[start:stop], eeg_srate, 0, spikes_coupled)
spike_radians_uncoupled = spike_lfp_radians(filtered_lfp[start:stop], eeg_srate, 0, spikes_uncoupled)

spike_radians_coupled_angles = (spike_radians_coupled + np.pi) % (2*np.pi) - np.pi
spike_radians_uncoupled_angles = (spike_radians_uncoupled + np.pi) % (2*np.pi) - np.pi

count_coupled, bin = np.histogram(spike_radians_coupled_angles, bins=16)
count_uncoupled, bin = np.histogram(spike_radians_uncoupled_angles, bins=16)

# By default plot density (frequency potentially misleading)
density = None
if density is None or density is True:
    # Area to assign each bin
    area_coupled = count_coupled / count_coupled.size
    area_uncoupled = count_uncoupled / count_uncoupled.size
    # Calculate corresponding bin radius
    radius_coupled = (area_coupled / np.pi)**.5
    radius_uncoupled = (area_uncoupled / np.pi)**.5
else:
    radius_coupled = count_coupled
    radius_uncoupled = count_uncoupled

widths = np.diff(bin)
#fig, ax = plt.subplots(2, figsize=(8,4), subplot_kw=dict(projection='polar'))

#rose_plot(ax[0], spike_radians_coupled)
#rose_plot(ax[1], spike_radians_uncoupled)

#plt.tight_layout()


plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(2, figsize=(4,6), subplot_kw=dict(projection='polar'))

ax[0].bar(bin[:-1], radius_coupled, zorder=1, align='edge', width=widths,
       edgecolor='C0', linewidth=0, color = 'g')
ax[0].set_theta_offset(0)

ax[0].set_yticks([])


ax[1].bar(bin[:-1], radius_uncoupled, zorder=1, align='edge', width=widths,
       edgecolor='C0', linewidth=0, color='r')
ax[1].set_theta_offset(0)

ax[1].set_yticks([])

plt.tight_layout()

savepath2 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure07\output\Fig7B.png'
plt.savefig(os.path.abspath(savepath2), dpi=600)

p_coupled = rayleightest(spike_radians_coupled)
p_uncoupled = rayleightest(spike_radians_uncoupled)
