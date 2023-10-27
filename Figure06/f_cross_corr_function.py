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

sys.path.insert(1, r'E:\OneDrive - UGent\python_functions')


from vies.lfp.filter import butter_bandpass_filter, butter_bandstop_filter
from vies.parse.neuron import load_neuronfile
from vies.lfp.lfp_analysis import spectrogram


from vies.parse.phy import phy_data

from vies.spike.extract_frequency_bins import extract_frequency_bins
from scipy.stats import pearsonr, spearmanr
from scipy.signal import correlate, correlation_lags
from vies.lfp.filter import movingaverage

from sklearn.cluster import KMeans


def eegpower_spikerate_crosscorr(eeg_path, eeg_srate, eeg_downsample_rate, eeg_stimstart, phy_path, unit_srate, spike_templates, seizure_start):
  
    eeg_path = eeg_path
    eeg_srate = eeg_srate
    down_samplerate_eeg = eeg_downsample_rate
    stim_start = eeg_stimstart
    
    phy_path = phy_path
    unit_srate = unit_srate
    spike_templates = spike_templates
    seizure_start = seizure_start
    
    ## Load + process EEG
    time_eeg, hip_eeg = load_neuronfile(eeg_path, srate=eeg_srate, channel=1, gain=500, inputrange=20)
    hip_eeg = butter_bandpass_filter(hip_eeg, 2, int(down_samplerate_eeg/2), eeg_srate, order=2)
    hip_eeg = butter_bandstop_filter(hip_eeg, 48, 52, eeg_srate, order=2) ## notch?
    hip_eeg = signal.resample(hip_eeg, int(len(hip_eeg)/int(eeg_srate/down_samplerate_eeg)))
    
    spectrotime, frequencies, spectraldata = spectrogram(hip_eeg, down_samplerate_eeg, windowlength=1, overlap=0.5, highfreq=1000)
    
    start_eeg = stim_start + 10
    stop_eeg = start_eeg + 50
    index_start_eeg = np.where(spectrotime == start_eeg)[0][0] + 1
    index_stop_eeg = np.where(spectrotime == stop_eeg)[0][0] + 1
    
    spectraldata = spectraldata[:, index_start_eeg:index_stop_eeg] 
    spectrotime = spectrotime[index_start_eeg:index_stop_eeg]
    
    new_time = np.zeros(int(np.ceil(len(spectrotime) / 2)))
    new_spectraldata = np.zeros((np.shape(spectraldata)[0], int(np.ceil(len(spectrotime) / 2))))
    for i in range(len(new_time)):
        if i == 0:
            new_time[i] = spectrotime[i]
            new_spectraldata[:,i] = spectraldata[:,i]
            
        elif i == len(new_time) - 1:
            start = int(i*2 - 1)
            new_time[i] = np.mean(spectrotime[start:])
            new_spectraldata[:,i] = np.mean(spectraldata[:, start:], axis=1)
        else:
            start = int(i*2 - 1)
            stop = int(start+3)
            new_time[i] = np.mean(spectrotime[start:stop])
            new_spectraldata[:,i] = np.mean(spectraldata[:, start:stop], axis=1)
    
    for i in range(np.shape(new_spectraldata)[0]):
        new_spectraldata[i,:] = movingaverage(new_spectraldata[i,:], 5)
    ## Load + process unit data
        
    unit_data = phy_data(phy_path, unit_srate)
    seizure_start = seizure_start + 10
    unit_start = seizure_start
    unit_stop = seizure_start + 50
    
    unit_bins = np.zeros((len(new_time), len(spike_templates)))
    for i in range(len(spike_templates)):
        spikes = unit_data.extract_spike_trains(spike_templates[i])
        unit_bins[:,i] = extract_frequency_bins(spikes, 1, unit_start, unit_stop)[1]
        unit_bins[:,i] = movingaverage(unit_bins[:,i], 5)
    
    corrs_all_freqs = np.zeros(((np.shape(new_spectraldata)[1]), np.shape(new_spectraldata)[0], len(spike_templates)))
    #corrs_spearman = np.zeros((np.shape(new_spectraldata)[0], len(spike_templates)))
    corrs_total = np.zeros(((np.shape(new_spectraldata)[1]), len(spike_templates)))
    
    total_amplitude = np.sum(np.sqrt(new_spectraldata), axis = 0)
    for i in range(len(spike_templates)):
        a = (total_amplitude - np.mean(total_amplitude)) / (np.std(total_amplitude) * len(total_amplitude))
        b = (unit_bins[:,i] - np.mean(unit_bins[:,i])) / (np.std(unit_bins[:,i]))
        corrs_total[:,i] = correlate(a, b, mode='same')
        for k in range(np.shape(new_spectraldata)[0]):
            a = (np.sqrt(new_spectraldata[k,:]) - np.mean(np.sqrt(new_spectraldata[k,:]))) / (np.std(np.sqrt(new_spectraldata[k,:])) * len(np.sqrt(new_spectraldata[k,:])))
            b = (unit_bins[:,i] - np.mean(unit_bins[:,i])) / (np.std(unit_bins[:,i]))
            corrs_all_freqs[:,k,i] = correlate(a, b, mode='same')
            #corrs_spearman[k,i] = spearmanr(np.sqrt(new_spectraldata[k,:]), unit_bins[:,i])[0] ** 2

    return corrs_all_freqs, corrs_total
