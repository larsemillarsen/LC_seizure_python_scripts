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
from vies.parse.phy import phy_data
from vies.lfp.filter import butter_bandpass_filter
from vies.spike.extract_frequency_bins import extract_frequency_bins
from vies.spike_lfp.lfp_phase_spike_coherence import spike_lfp_radians
from astropy.stats.circstats import rayleightest

from vies.general.roseplot import rose_plot
from scipy.signal import hilbert

### LL023
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL007\Seizure0.mat'
eeg_srate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'

seizure_start = 1193.42


time_eeg, hip_eeg = load_neuronfile(eeg_path, srate=eeg_srate, channel=1, gain=200, inputrange=20)

filter_low = 200
filter_high = 4500

start = int((eeg_stimstart + 10) * eeg_srate)
stop = int((eeg_stimstart + 30) * eeg_srate)


start = int((eeg_stimstart - 10) * eeg_srate)
stop = int((eeg_stimstart + 30) * eeg_srate)


### FULL SEIZURE
fig, ax = plt.subplots(1, figsize=(6,3))
ax.plot(time_eeg[start:stop], hip_eeg[start:stop], 'k', linewidth = 0.25)
ax.vlines(120, np.min(hip_eeg[start:stop]), np.max(hip_eeg[start:stop]), colors='r', linewidth=2, linestyles='dashed', zorder=3)
ax.vlines(130, np.min(hip_eeg[start:stop]), np.max(hip_eeg[start:stop]), colors='r', linewidth=2, linestyles='dashed', zorder=3)
ax.axvspan(70+60, 75+60, color='b', alpha=0.6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticks([])
ax.set_yticklabels([]) 
ax.set_xticks([110, 120, 130, 140, 150])
ax.set_xticklabels([-10, 0, 10, 20, 30]) 
ax.set_xlabel('Time (s)',fontsize=14, fontweight='bold')

plt.tight_layout()
savepath2 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure07\output\Fig7E.png'
plt.savefig(os.path.abspath(savepath2), dpi=600)


