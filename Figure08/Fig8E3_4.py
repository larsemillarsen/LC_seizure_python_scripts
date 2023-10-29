# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 09:42:44 2022

@author: User
"""

#import os
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(1, r'C:\Users\scaestec\OneDrive - UGent\PhD\Python\Phyton_functions')

from vies.parse.pyphotometry import import_ppd
from vies.parse.neuron import load_neuronfile
from vies.lfp.lfp_analysis import spectrogram
#from vies.lfp.filter import movingaverage
from scipy.signal import correlate, correlation_lags

directory = r'E:\Manuscript_analysis_files\data\GRABne'
savepath_grp = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure08\output\Fig8E_3.png'
savepath_boxplot = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure08\output\Fig8E_4.png'

sc025 = {}
sc038 = {}
sc040 = {}
sc041 = {}
sc044 = {} 
sc045 = {}
sc046 = {}
sc047 = {}
sc048 = {}
sc054 = {}


number_of_animals = 10
lowcut = 0.001
highcut = 2

sc025['file'] = directory + '\SC025_awake-2022-06-29-112116.ppd'
sc038['file'] = directory + '\SC038_seiz-2023-02-01-132816.ppd'
sc040['file'] = directory + '\SC040_thresh-2023-02-02-120825.ppd'
sc041['file'] = directory + '\SC041_seiz-2023-02-02-152751.ppd'
sc044['file'] = directory + '\SC044_seiz-2023-02-03-115555.ppd'
sc045['file'] = directory + '\SC045_seiz-2023-02-03-135432.ppd'
sc046['file'] = directory + '\SC046_thresh-2023-02-06-121340.ppd'
sc047['file'] = directory + '\SC047_seiz-2023-02-06-152843.ppd'
sc048['file'] = directory + '\SC048_thresh-2023-02-07-120254.ppd'
sc054['file'] = directory + '\SC054_thresh-2023-02-07-142631.ppd'

sc025['data'] = import_ppd(sc025['file'], low_pass=1, f_type='bandpass')
sc038['data'] = import_ppd(sc038['file'], low_pass=1, f_type='bandpass')
sc040['data'] = import_ppd(sc040['file'], low_pass=1, f_type='bandpass')
sc041['data'] = import_ppd(sc041['file'], low_pass=1, f_type='bandpass')
sc044['data'] = import_ppd(sc044['file'], low_pass=1, f_type='bandpass')
sc045['data'] = import_ppd(sc045['file'], low_pass=1, f_type='bandpass')
sc046['data'] = import_ppd(sc046['file'], low_pass=1, f_type='bandpass')
sc047['data'] = import_ppd(sc047['file'], low_pass=1, f_type='bandpass')
sc048['data'] = import_ppd(sc048['file'], low_pass=1, f_type='bandpass')
sc054['data'] = import_ppd(sc054['file'], low_pass=1, f_type='bandpass')

sc025['seizure_start'] = np.where(np.diff(sc025['data']['digital_1'])==1)[0][-1]
sc038['seizure_start'] = np.where(np.diff(sc038['data']['digital_1'])==1)[0][-1]
sc040['seizure_start'] = np.where(np.diff(sc040['data']['digital_1'])==1)[0][-1]
sc041['seizure_start'] = np.where(np.diff(sc041['data']['digital_1'])==1)[0][-1]
sc044['seizure_start'] = np.where(np.diff(sc044['data']['digital_1'])==1)[0][-1]
sc045['seizure_start'] = np.where(np.diff(sc045['data']['digital_1'])==1)[0][-1]
sc046['seizure_start'] = np.where(np.diff(sc046['data']['digital_1'])==1)[0][-1]
sc047['seizure_start'] = np.where(np.diff(sc047['data']['digital_1'])==1)[0][-1]
sc048['seizure_start'] = np.where(np.diff(sc048['data']['digital_1'])==1)[0][-1]
sc054['seizure_start'] = np.where(np.diff(sc054['data']['digital_1'])==1)[0][-1]

sc025['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC025\SC025_Seizure1_awake220629A0011.mat'
sc038['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC038\SC038_seizure1230201A0009.mat'
sc040['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC040\SC040_threshold50230202A0001.mat'
sc041['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC041\SC041_seizure1230202A0004.mat'
sc044['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC044\SC044_seizure1230203A0001.mat'
sc045['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC045\SC045_seizure230203A0003.mat'
sc046['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC046\SC046_threshold100230206A0001.mat'
sc047['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC047\SC047_seizure1230206A0005.mat'
sc048['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC048\SC048_threshold100230207A0001.mat'
sc054['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC054\SC054_threshold50230207A0002.mat'

sc025['eeg_data'] = load_neuronfile(sc025['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc038['eeg_data'] = load_neuronfile(sc038['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc040['eeg_data'] = load_neuronfile(sc040['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc041['eeg_data'] = load_neuronfile(sc041['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc044['eeg_data'] = load_neuronfile(sc044['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc045['eeg_data'] = load_neuronfile(sc045['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc046['eeg_data'] = load_neuronfile(sc046['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc047['eeg_data'] = load_neuronfile(sc047['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc048['eeg_data'] = load_neuronfile(sc048['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc054['eeg_data'] = load_neuronfile(sc054['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)

start = 60*2
stop = start + 50*2

sc025['eeg_spectrogram'] = np.sqrt(spectrogram(sc025['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc038['eeg_spectrogram'] = np.sqrt(spectrogram(sc038['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc040['eeg_spectrogram'] = np.sqrt(spectrogram(sc040['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc041['eeg_spectrogram'] = np.sqrt(spectrogram(sc041['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc044['eeg_spectrogram'] = np.sqrt(spectrogram(sc044['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc045['eeg_spectrogram'] = np.sqrt(spectrogram(sc045['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc046['eeg_spectrogram'] = np.sqrt(spectrogram(sc046['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc047['eeg_spectrogram'] = np.sqrt(spectrogram(sc047['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc048['eeg_spectrogram'] = np.sqrt(spectrogram(sc048['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc054['eeg_spectrogram'] = np.sqrt(spectrogram(sc054['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])

eeg_amplitude_data = np.zeros((50, number_of_animals))

for i in range(np.shape(eeg_amplitude_data)[0]):
    if i == 0:
        eeg_amplitude_data[i,0] = np.sum(sc025['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,1] = np.sum(sc038['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,2] = np.sum(sc040['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,3] = np.sum(sc041['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,4] = np.sum(sc044['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,5] = np.sum(sc045['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,6] = np.sum(sc046['eeg_spectrogram'][:,i], axis=0)   
        eeg_amplitude_data[i,7] = np.sum(sc047['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,8] = np.sum(sc048['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,9] = np.sum(sc054['eeg_spectrogram'][:,i], axis=0)
               
        
    elif i == np.shape(eeg_amplitude_data)[0] - 1:
        start = int(i*2 - 1)
        eeg_amplitude_data[i,0] = np.sum(sc025['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,1] = np.sum(sc038['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,2] = np.sum(sc040['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,3] = np.sum(sc041['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,4] = np.sum(sc044['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,5] = np.sum(sc045['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,6] = np.sum(sc046['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,7] = np.sum(sc047['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,8] = np.sum(sc048['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,9] = np.sum(sc054['eeg_spectrogram'][:,start], axis=0)
              
    else:
        start = int(i*2 - 1)
        stop = int(start+3)
        eeg_amplitude_data[i,0] = np.sum(np.mean(sc025['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,1] = np.sum(np.mean(sc038['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,2] = np.sum(np.mean(sc040['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,3] = np.sum(np.mean(sc041['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,4] = np.sum(np.mean(sc044['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,5] = np.sum(np.mean(sc045['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,6] = np.sum(np.mean(sc046['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,7] = np.sum(np.mean(sc047['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,8] = np.sum(np.mean(sc048['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,9] = np.sum(np.mean(sc054['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        

#for i in range(np.shape(eeg_amplitude_data)[1]):
#    eeg_amplitude_data[:,i] = movingaverage(eeg_amplitude_data[:,i], 5)


scaled_eeg_amplitude_data = np.zeros((np.shape(eeg_amplitude_data)[0], np.shape(eeg_amplitude_data)[1]))
for i in range(np.shape(scaled_eeg_amplitude_data)[1]):
    scaled_eeg_amplitude_data[:, i] = (eeg_amplitude_data[:, i] - np.min(eeg_amplitude_data[:, i])) / (np.max(eeg_amplitude_data[:, i], axis=0) - np.min(eeg_amplitude_data[:, i])) 
    

sz_grab_data = np.zeros((120*130, number_of_animals))


sz_grab_data[:,0] = sc025['data']['analog_1'][int(sc025['seizure_start'] - 130*60):int(sc025['seizure_start'] + 130*60)]
sz_grab_data[:,1] = sc038['data']['analog_1'][int(sc038['seizure_start'] - 130*60):int(sc038['seizure_start'] + 130*60)]
sz_grab_data[:,2] = sc040['data']['analog_1'][int(sc040['seizure_start'] - 130*60):int(sc040['seizure_start'] + 130*60)]
sz_grab_data[:,3] = sc041['data']['analog_1'][int(sc041['seizure_start'] - 130*60):int(sc041['seizure_start'] + 130*60)]
sz_grab_data[:,4] = sc044['data']['analog_1'][int(sc044['seizure_start'] - 130*60):int(sc044['seizure_start'] + 130*60)]
sz_grab_data[:,5] = sc045['data']['analog_1'][int(sc045['seizure_start'] - 130*60):int(sc045['seizure_start'] + 130*60)]
sz_grab_data[:,6] = sc046['data']['analog_1'][int(sc046['seizure_start'] - 130*60):int(sc046['seizure_start'] + 130*60)]
sz_grab_data[:,7] = sc047['data']['analog_1'][int(sc047['seizure_start'] - 130*60):int(sc047['seizure_start'] + 130*60)]
sz_grab_data[:,8] = sc048['data']['analog_1'][int(sc048['seizure_start'] - 130*60):int(sc048['seizure_start'] + 130*60)]
sz_grab_data[:,9] = sc054['data']['analog_1'][int(sc054['seizure_start'] - 130*60):int(sc054['seizure_start'] + 130*60)]

baseline = int(60 * 130)
sz_start = int(70 * 130)
base_means = np.mean(sz_grab_data[0:baseline], 0)

sz_grab_data[:,0] = sc025['data']['analog_1_filt'][int(sc025['seizure_start'] - 130*60):int(sc025['seizure_start'] + 130*60)] + base_means[0]
sz_grab_data[:,1] = sc038['data']['analog_1_filt'][int(sc038['seizure_start'] - 130*60):int(sc038['seizure_start'] + 130*60)] + base_means[1]
sz_grab_data[:,2] = sc040['data']['analog_1_filt'][int(sc040['seizure_start'] - 130*60):int(sc040['seizure_start'] + 130*60)] + base_means[2]
sz_grab_data[:,3] = sc041['data']['analog_1_filt'][int(sc041['seizure_start'] - 130*60):int(sc041['seizure_start'] + 130*60)] + base_means[3]
sz_grab_data[:,4] = sc044['data']['analog_1_filt'][int(sc044['seizure_start'] - 130*60):int(sc044['seizure_start'] + 130*60)] + base_means[4]
sz_grab_data[:,5] = sc045['data']['analog_1_filt'][int(sc045['seizure_start'] - 130*60):int(sc045['seizure_start'] + 130*60)] + base_means[5]
sz_grab_data[:,6] = sc046['data']['analog_1_filt'][int(sc046['seizure_start'] - 130*60):int(sc046['seizure_start'] + 130*60)] + base_means[6]
sz_grab_data[:,7] = sc047['data']['analog_1_filt'][int(sc047['seizure_start'] - 130*60):int(sc047['seizure_start'] + 130*60)] + base_means[7]
sz_grab_data[:,8] = sc048['data']['analog_1_filt'][int(sc048['seizure_start'] - 130*60):int(sc048['seizure_start'] + 130*60)] + base_means[8]
sz_grab_data[:,9] = sc054['data']['analog_1_filt'][int(sc054['seizure_start'] - 130*60):int(sc054['seizure_start'] + 130*60)] + base_means[9]

sz_grab_data_normalized = np.zeros((120*130, number_of_animals))
sz_grab_data_normalized_1s = np.zeros((120, number_of_animals))

for i in range(np.shape(sz_grab_data_normalized)[1]):
    baseline = np.mean(sz_grab_data[0:int(50*130),i], axis = 0)
    sz_grab_data_normalized[:,i] = ((sz_grab_data[:,i] - baseline) / baseline) * 100
    
for i in range(np.shape(sz_grab_data_normalized_1s)[0]):
    start = int(i*130)
    stop = int(start + 130)
    sz_grab_data_normalized_1s[i,:] = np.mean(sz_grab_data_normalized[start:stop,:], axis=0) 

#for i in range(np.shape(sz_grab_data_normalized_1s)[1]):
#    sz_grab_data_normalized_1s[:,i] = movingaverage(sz_grab_data_normalized_1s[:,i], 5)
   
sz_grab_data_normalized_1s = sz_grab_data_normalized_1s[70:,:]

for i in range(np.shape(sz_grab_data_normalized_1s)[1]):
    sz_grab_data_normalized_1s[:,i] = (sz_grab_data_normalized_1s[:,i] - np.min(sz_grab_data_normalized_1s[:,i])) / (np.max(sz_grab_data_normalized_1s[:,i]) - np.min(sz_grab_data_normalized_1s[:,i]))



        
correlations = np.zeros((50,10))
test = correlate(sz_grab_data_normalized_1s[:,0], scaled_eeg_amplitude_data[:,0], mode='same')
for i in range(np.shape(sz_grab_data_normalized_1s)[1]):
    a = (scaled_eeg_amplitude_data[:,i] - np.mean(scaled_eeg_amplitude_data[:,i], axis=0)) / (np.std(scaled_eeg_amplitude_data[:,i],axis = 0) * len(scaled_eeg_amplitude_data[:,i]))
    b = (sz_grab_data_normalized_1s[:,i] - np.mean(sz_grab_data_normalized_1s[:,i], axis=0)) / (np.std(sz_grab_data_normalized_1s[:,i]))
    #plt.plot(b)
    correlations[:,i] = correlate(a, b, mode='same')
    

#plt.plot(scaled_eeg_amplitude_data[:,6])
#plt.plot(sz_grab_data_normalized_1s[:,6])

x = correlation_lags(50, 50, mode='same')

lags = np.zeros((10))
peak_r = np.zeros((10))
for i in range(np.shape(correlations)[1]):
    #plt.plot(x, correlations[:,i])
    index = np.where(correlations[:,i] == np.max(correlations[:,i]))[0]
    lags[i] = x[index]
    peak_r[i] = np.max(correlations[:,i])
    

mean_corr = np.mean(correlations, axis=1)
std_corr = np.std(correlations, axis=1)
x_lag = correlation_lags(50, 50, mode='same')

### GRP PLOTS

x = np.linspace(0.65, 0.82, number_of_animals)
meanprops = {'linestyle': 'solid',
            'color':'black',
            'linewidth':4}

medianprops  = {'linestyle': 'solid',
            'color':'red',
            'linewidth':4}

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(3,3))
      
ax[0].scatter(x, peak_r)
ax[0].boxplot(peak_r, showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)
ax[0].set_xlim(0.6 ,1.15)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_linewidth(2)
ax[0].set_ylabel('xcorr',fontsize=16, fontweight='bold')
ax[0].set_xticks([])
ax[0].set_xticklabels([])
ax[0].set_yticks([0.4, 0.6, 0.8, 1.0])
ax[0].set_yticklabels([0.4, 0.6, 0.8, 1.0], fontsize=12)
#ax[0].text(0.83, 115, '*', fontsize = 25, fontweight='bold')

ax[1].scatter(x, lags)
ax[1].boxplot(lags, showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)
ax[1].set_xlim(0.6 ,1.15)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_linewidth(2)
ax[1].set_ylabel('Lag (s)',fontsize=16, fontweight='bold')
#ax[1].hlines(1.96, 0.6, 1.4, colors='r', linestyles='--')
ax[1].set_xticks([])
ax[1].set_xticklabels([])
ax[1].set_yticks([10, 5, 0, -5, -10, -15])
ax[1].set_yticklabels([10, 5, 0, -5, -10, -15], fontsize=12)
plt.tight_layout()
plt.savefig(savepath_boxplot, dpi=600)
#plt.close()

    

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(3,3))


ax.plot(x_lag, mean_corr, color='k')
ax.fill_between(x_lag, mean_corr-std_corr, mean_corr+std_corr, color='gray')
ax.vlines(0,-1,1, colors='r', linestyles='--', zorder=3)
ax.set_xticks([-25, 0, 25])
ax.set_xticklabels([-25, 0, 25], fontsize = 12)
ax.set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
ax.set_yticklabels([-1.0, -0.5, 0, 0.5, 1.0], fontsize = 12)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.set_xlabel('Lag (s)',fontsize=14, fontweight='bold')
ax.set_ylabel('xcorr',fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig(savepath_grp, dpi=600)



