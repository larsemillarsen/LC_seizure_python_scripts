# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 12:36:46 2022

@author: User
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, r'C:\Users\User\OneDrive - UGent\python_functions')

from vies.parse.neuron import load_neuronfile
from vies.lfp.filter import butter_bandpass_filter
from scipy.signal import find_peaks
from vies.parse.phy import phy_data
#from vies.spike_analysis.extract_frequency_bins import extract_frequency_bins
from scipy.stats import ks_2samp


#eeg_path = r'E:\LCSeizureData\LCSeizureNIDAQ\LL023\Seizures\4_Seizure300_2201218A0004.mat'
#eeg_srate = 10000
#eeg_stimstart = 60
#phy_path = r'E:\LCSeizureData\LCSeizure_OpenEphys\LL023\2020-12-18_13-43-19\Analysis108\LL023\LL023-final.GUI'
#eeg_stimstart = 60
#unit_srate = 30000
#templates = [4]
#seizure_start = 2061.56

def popspike_coupling(eeg_path, eeg_srate, eeg_stimstart, phy_path, unit_srate, templates, seizure_start, window=1.0):

    filter_low = 200
    filter_high = 4500
    unit_srate = unit_srate
    eeg_srate = eeg_srate
    templates = templates
    
    time_eeg, hip_eeg = load_neuronfile(eeg_path, srate=eeg_srate, channel=1, gain=200, inputrange=20)
    filtered_lfp = butter_bandpass_filter(hip_eeg, filter_low, filter_high, eeg_srate, order=1)
    
    sd_signal = np.std(filtered_lfp)
    peaks = find_peaks(-filtered_lfp, height=5*sd_signal, distance=0.01*eeg_srate)
    peaks_seconds = peaks[0] / eeg_srate
         
    peaks_seconds = peaks_seconds - eeg_stimstart
    peaks_seconds = peaks_seconds[peaks_seconds>10]
    peaks_seconds = peaks_seconds[peaks_seconds<60]   
    path_phy_data = phy_path
    
    data = phy_data(path_phy_data, unit_srate)
    
    width_histogram = window # in seconds
    test_outcome = np.zeros(len(templates))
    counter = 0
    
    dict_spikes = []
    for t in templates:
        spike_times = data.extract_spike_trains(t) - seizure_start # aligns spikes with LFP trace
        
        for i in range(len(peaks_seconds)):
            #print(i)
            if i == 0:
                #print('test')
                rel_spikes = spike_times - peaks_seconds[i]
                
                rel_spikes = rel_spikes[rel_spikes > -width_histogram/2]
                rel_spikes = rel_spikes[rel_spikes < width_histogram/2]
                
                spikes = rel_spikes
            else:
                rel_spikes = spike_times - peaks_seconds[i]
                
                rel_spikes = rel_spikes[rel_spikes >= -width_histogram/2]
                rel_spikes = rel_spikes[rel_spikes <= width_histogram/2]
                
                spikes = np.concatenate((spikes, rel_spikes), axis=0)       
        
        dict_spikes.append(spikes)
        
        
        #uniform_dist=np.random.uniform(low=-width_histogram/2, high=width_histogram/2, size=len(spikes))    
        uniform_dist = np.linspace(start = -width_histogram/2, stop=width_histogram/2, num=len(spikes))
        
        test = ks_2samp(spikes, uniform_dist)
        test_outcome[counter] = test[1]
    
        counter = counter + 1



    #print(histogram[0])
    return test_outcome, dict_spikes

### LL023

#def popspike_coupling(eeg_path, eeg_srate, phy_path, unit_srate, templates, seizure_start):
#eeg_path = r'E:\LCSeizureData\LCSeizureNIDAQ\LL023\Seizures\4_Seizure300_2201218A0004.mat'
#eeg_srate = 10000
#eeg_stimstart = 60
#phy_path = r'E:\LCSeizureData\LCSeizure_OpenEphys\LL023\2020-12-18_13-43-19\Analysis108\LL023\LL023-final.GUI'
#eeg_stimstart = 60
#unit_srate = 30000
#templates = [4]
#seizure_start = 2061.56

#test_function, histogram = popspike_coupling(eeg_path, eeg_srate, eeg_stimstart, phy_path, unit_srate, templates, seizure_start)
