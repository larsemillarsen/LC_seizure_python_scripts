# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 20:00:53 2021

@author: Lars
"""


import sys
sys.path.insert(1, r'C:\Users\llarsen\OneDrive - ugentbe\python_functions')
import numpy as np
from astropy.stats.circstats import rayleightest, circmean

from vies.spike_lfp.lfp_phase_spike_coherence import spike_lfp_radians
from vies.lfp.filter import butter_bandpass_filter
from vies.parse.neuron import load_neuronfile
from vies.parse.phy import phy_data


def spike_lfpphase_coupling(lfp_path, lfp_srate, lfp_stimstart, phy_path, unit_srate, templates, seizure_start):
    eeg_srate = lfp_srate
    eeg_stimstart = lfp_stimstart
    
    eeg_path = lfp_path
    time_eeg, hip_eeg = load_neuronfile(eeg_path, srate=eeg_srate, channel=1, gain=200, inputrange=20)
    
    start = int((eeg_stimstart - 60) * eeg_srate)
    stop = int((eeg_stimstart) * eeg_srate)
    
    seizure_lfp = hip_eeg[start:stop]
    
    #### Make frequency bands
    test = [1]
    while test[-1]<40:
        test.append(test[-1] * 2)
        
    filter_low = test[:-1]
    filter_high = test[1:]
    
    ### LOAD spike data
    path_phy_data = phy_path
    unit_srate = unit_srate
    data = phy_data(path_phy_data, unit_srate)
    templates = templates
    
    test_results = np.zeros((len(filter_low), len(templates)))
    test_results_mean = np.zeros((len(filter_low), len(templates)))
    
    counter = 0
    for template in templates:
        spikes = data.extract_spike_trains(template)
        
        rel_lfp_start = seizure_start + 10
        n_spikes = len(spike_lfp_radians(seizure_lfp, eeg_srate, rel_lfp_start, spikes))
        spike_radians = np.zeros((n_spikes, len(filter_low)))
        
        for i in range(len(filter_low)):
            filtered_lfp = butter_bandpass_filter(seizure_lfp, filter_low[i], filter_high[i], eeg_srate, order=1)
            spike_radians[:,i] = spike_lfp_radians(filtered_lfp, eeg_srate, rel_lfp_start, spikes)

        for filter_setting in range(len(filter_low)):
            test_results[filter_setting, counter] = rayleightest(spike_radians[:,filter_setting])
            test_results_mean[filter_setting, counter] = circmean(spike_radians[:,filter_setting])
        counter = counter + 1

    return test_results, test_results_mean
