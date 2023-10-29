# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 09:42:44 2022

@author: User
"""

#import os
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(1, r'E:\Manuscript_analysis_files')

from vies.parse.pyphotometry import import_ppd, open_pyphotometry_csv
from vies.parse.neuron import load_neuronfile
from vies.lfp.lfp_analysis import spectrogram
#from vies.lfp.filter import movingaverage
from scipy.signal import correlate, correlation_lags

directory = r'E:\Manuscript_analysis_files\data\GRABne'
savepath_grp = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure08\output\Fig8D_3.png'
savepath_boxplot = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure08\output\Fig8D_4.png'

sc016 = {}
sc024 = {}
sc025 = {}
sc026 = {}
#sc027 = {}
sc028 = {} 
sc030 = {}
sc031 = {}
#sc032 = {}
sc033 = {}
sc034 = {}
#sc036 = {}
#sc037 = {}
sc038 = {}
sc039 = {}
sc040 = {}
sc041 = {}
#sc043 = {}
sc044 = {}
sc045 = {}
sc046 = {}
sc047 = {}
sc048 = {}
sc049 = {}


number_of_animals = 19
lowcut = 0.001
highcut = 2

sc016['file'] = directory + '\SC016_seiz1-2022-03-24-123537.csv'
sc024['file'] = directory + '\sc024-2022-06-14-142418.csv'
sc025['file'] = directory + '\sc025-2022-06-15-152403.csv'
sc026['file'] = directory + '\sc026-2022-06-16-125500.csv'
#sc027['file'] = directory + '\sc027-2022-06-17-140853.csv'
sc028['file'] = directory + '\sc028-2022-06-20-150958.csv'
sc030['file'] = directory + '\sc030-2022-06-21-145138.csv'
sc031['file'] = directory + '\sc031-2022-06-22-162759.ppd'
#sc032['file'] = directory + '\SC032-2022-06-23-141205.csv'
sc033['file'] = directory + '\SC033-2022-06-24-153503.ppd'
sc034['file'] = directory + '\SC034_acute-2022-11-08-155546.ppd'
#sc036['file'] = directory + '\SC036_acute-2022-11-17-161742.ppd'
#sc037['file'] = directory + '\SC037_acute-2022-11-18-124841.ppd'
sc038['file'] = directory + '\SC038_acute-2022-11-21-130818.ppd'
sc039['file'] = directory + '\SC039_acute-2022-11-25-114242.ppd'
sc040['file'] = directory + '\SC040_acute-2022-11-29-115358.ppd'
sc041['file'] = directory + '\SC041_acute-2022-12-02-113303.ppd'
#sc043['file'] = directory + '\SC043_acute-2022-11-22-135339.ppd'
sc044['file'] = directory + '\SC044_acute-2022-11-23-121522.ppd'
sc045['file'] = directory + '\SC045_acute-2022-11-24-125435.ppd'
sc046['file'] = directory + '\SC046_acute-2022-11-28-112930.ppd'
sc047['file'] = directory + '\SC047_acute-2022-11-30-114328.ppd'
sc048['file'] = directory + '\SC048_acute-2022-12-01-113127.ppd'
sc049['file'] = directory + '\SC049_acute-2022-12-05-115326.ppd'

sc016['data'] = open_pyphotometry_csv(sc016['file'], [1], [1], srate=130, delimiter=',', low_pass=1, f_type='bandpass')
sc024['data'] = open_pyphotometry_csv(sc024['file'], [1], [1], srate=130, delimiter=',', low_pass=1, f_type='bandpass')
sc025['data'] = open_pyphotometry_csv(sc025['file'], [1], [1], srate=130, delimiter=',', low_pass=1, f_type='bandpass')
sc026['data'] = open_pyphotometry_csv(sc026['file'], [1], [1], srate=130, delimiter=',', low_pass=1, f_type='bandpass')
#sc027['data'] = open_pyphotometry_csv(sc027['file'], [1], [1], srate=130, delimiter=',', low_pass=1, f_type='bandpass')
sc028['data'] = open_pyphotometry_csv(sc028['file'], [1], [1], srate=130, delimiter=',', low_pass=1, f_type='bandpass')
sc030['data'] = open_pyphotometry_csv(sc030['file'], [1], [1], srate=130, delimiter=',', low_pass=1, f_type='bandpass')
sc031['data'] = import_ppd(sc031['file'], low_pass=1, f_type='bandpass')
#sc032['data'] = open_pyphotometry_csv(sc032['file'], [1], [1], srate=130, delimiter=',', low_pass=1, f_type='bandpass')
sc033['data'] = import_ppd(sc033['file'], low_pass=1, f_type='bandpass')
sc034['data'] = import_ppd(sc034['file'], low_pass=1, f_type='bandpass')
#sc036['data'] = import_ppd(sc036['file'], low_pass=1, f_type='bandpass')
#sc037['data'] = import_ppd(sc037['file'], low_pass=1, f_type='bandpass')
sc038['data'] = import_ppd(sc038['file'], low_pass=1, f_type='bandpass')
sc039['data'] = import_ppd(sc039['file'], low_pass=1, f_type='bandpass')
sc040['data'] = import_ppd(sc040['file'], low_pass=1, f_type='bandpass')
sc041['data'] = import_ppd(sc041['file'], low_pass=1, f_type='bandpass')
#sc043['data'] = import_ppd(sc043['file'], low_pass=1, f_type='bandpass')
sc044['data'] = import_ppd(sc044['file'], low_pass=1, f_type='bandpass')
sc045['data'] = import_ppd(sc045['file'], low_pass=1, f_type='bandpass')
sc046['data'] = import_ppd(sc046['file'], low_pass=1, f_type='bandpass')
sc047['data'] = import_ppd(sc047['file'], low_pass=1, f_type='bandpass')
sc048['data'] = import_ppd(sc048['file'], low_pass=1, f_type='bandpass')
sc049['data'] = import_ppd(sc049['file'], low_pass=1, f_type='bandpass')


sc016['seizure_start'] = 60*130
sc024['seizure_start'] = np.where(np.diff(sc024['data'][2])==1)[0][-1]
sc025['seizure_start'] = np.where(np.diff(sc025['data'][2])==1)[0][-2]
sc026['seizure_start'] = np.where(np.diff(sc026['data'][2])==1)[0][-1]
sc028['seizure_start'] = np.where(np.diff(sc028['data'][2])==1)[0][-1]
sc030['seizure_start'] = np.where(np.diff(sc030['data'][2])==1)[0][-1]
sc031['seizure_start'] = np.where(np.diff(sc031['data']['digital_1'])==1)[0][-1]
#sc032['seizure_start'] = np.where(np.diff(sc030['data'][2])==1)[0][-1]
sc033['seizure_start'] = np.where(np.diff(sc033['data']['digital_1'])==1)[0][-1]
sc034['seizure_start'] = np.where(np.diff(sc034['data']['digital_1'])==1)[0][-1]
#sc036['seizure_start'] = np.where(np.diff(sc031['data']['digital_1'])==1)[0][-1]
#sc037['seizure_start'] = np.where(np.diff(sc031['data']['digital_1'])==1)[0][-1]
sc038['seizure_start'] = np.where(np.diff(sc038['data']['digital_1'])==1)[0][-1]
sc039['seizure_start'] = np.where(np.diff(sc039['data']['digital_1'])==1)[0][-1]
sc040['seizure_start'] = np.where(np.diff(sc040['data']['digital_1'])==1)[0][-1]
sc041['seizure_start'] = np.where(np.diff(sc041['data']['digital_1'])==1)[0][-1]
#sc043['seizure_start'] = np.where(np.diff(sc031['data']['digital_1'])==1)[0][-1]
sc044['seizure_start'] = np.where(np.diff(sc044['data']['digital_1'])==1)[0][-1]
sc045['seizure_start'] = np.where(np.diff(sc045['data']['digital_1'])==1)[0][-1]
sc046['seizure_start'] = np.where(np.diff(sc046['data']['digital_1'])==1)[0][-1]
sc047['seizure_start'] = np.where(np.diff(sc047['data']['digital_1'])==1)[0][-1]
sc048['seizure_start'] = np.where(np.diff(sc048['data']['digital_1'])==1)[0][-1]
sc049['seizure_start'] = np.where(np.diff(sc049['data']['digital_1'])==1)[0][-1]

sc016['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC016\Seizure1220324A0003.mat'
sc024['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC024\SC024_Seizure1220614A0005.mat'
sc025['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC025\SC025_Seizure1220615A0006.mat'
sc026['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC026\SC026_Seizure1220616A0005.mat'
sc028['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC028\SC028_Seizure1220620A0006.mat'
sc030['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC030\SC030_Seizure1220621A0005.mat'
sc031['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC031\SC031_Seizure1220622A0006.mat'
#sc032['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC032\SC032_Seizure1220623A0000.mat'
sc033['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC033\SC033_Seizure1220624A0006.mat'
sc034['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC034\SC034_seizure1221108A0007.mat'
#sc036['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC036\SC036_seizure1221117A0004.mat'
#sc037['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC037\SC037_seizure221118A0005.mat'
sc038['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC038\SC038_seizure1221121A0007.mat'
sc039['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC039\SC039_seizure1221125A0006.mat'
sc040['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC040\SC040_seizure221129A0007.mat'
sc041['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC041\SC041_seizure1221202A0007.mat'
sc044['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC044\SC044_seizure1221123A0008.mat'
sc045['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC045\SC045_Seizure1221124A0006.mat'
sc046['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC046\SC046_seizure1221128A0006.mat'
sc047['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC047\SC047_seizure1221130A0007.mat'
sc048['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC048\SC048_seizure1221201A0007.mat'
sc049['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC049\SC049_seizure1221205A0007.mat'


sc016['eeg_data'] = load_neuronfile(sc016['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc024['eeg_data'] = load_neuronfile(sc024['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc025['eeg_data'] = load_neuronfile(sc025['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc026['eeg_data'] = load_neuronfile(sc026['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc028['eeg_data'] = load_neuronfile(sc028['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc030['eeg_data'] = load_neuronfile(sc030['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc031['eeg_data'] = load_neuronfile(sc031['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
#sc032['eeg_data'] = load_neuronfile(sc032['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc033['eeg_data'] = load_neuronfile(sc033['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc034['eeg_data'] = load_neuronfile(sc034['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
#sc036['eeg_data'] = load_neuronfile(sc036['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
#sc037['eeg_data'] = load_neuronfile(sc037['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc038['eeg_data'] = load_neuronfile(sc038['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc039['eeg_data'] = load_neuronfile(sc039['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc040['eeg_data'] = load_neuronfile(sc040['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc041['eeg_data'] = load_neuronfile(sc041['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc044['eeg_data'] = load_neuronfile(sc044['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc045['eeg_data'] = load_neuronfile(sc045['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc046['eeg_data'] = load_neuronfile(sc046['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc047['eeg_data'] = load_neuronfile(sc047['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc048['eeg_data'] = load_neuronfile(sc048['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)
sc049['eeg_data'] = load_neuronfile(sc049['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)

start = 60*2
stop = start + 50*2
sc016['eeg_spectrogram'] = np.sqrt(spectrogram(sc016['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,int(start+10*2):int(stop+10*2)])
sc024['eeg_spectrogram'] = np.sqrt(spectrogram(sc024['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc025['eeg_spectrogram'] = np.sqrt(spectrogram(sc025['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc026['eeg_spectrogram'] = np.sqrt(spectrogram(sc026['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc028['eeg_spectrogram'] = np.sqrt(spectrogram(sc028['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc030['eeg_spectrogram'] = np.sqrt(spectrogram(sc030['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc031['eeg_spectrogram'] = np.sqrt(spectrogram(sc031['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
#sc032['eeg_spectrogram'] = np.sqrt(spectrogram(sc032['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc033['eeg_spectrogram'] = np.sqrt(spectrogram(sc033['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc034['eeg_spectrogram'] = np.sqrt(spectrogram(sc034['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
#sc036['eeg_spectrogram'] = np.sqrt(spectrogram(sc036['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
#sc037['eeg_spectrogram'] = np.sqrt(spectrogram(sc037['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc038['eeg_spectrogram'] = np.sqrt(spectrogram(sc038['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc039['eeg_spectrogram'] = np.sqrt(spectrogram(sc039['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc040['eeg_spectrogram'] = np.sqrt(spectrogram(sc040['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc041['eeg_spectrogram'] = np.sqrt(spectrogram(sc041['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc044['eeg_spectrogram'] = np.sqrt(spectrogram(sc044['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc045['eeg_spectrogram'] = np.sqrt(spectrogram(sc045['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc046['eeg_spectrogram'] = np.sqrt(spectrogram(sc046['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc047['eeg_spectrogram'] = np.sqrt(spectrogram(sc047['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc048['eeg_spectrogram'] = np.sqrt(spectrogram(sc048['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])
sc049['eeg_spectrogram'] = np.sqrt(spectrogram(sc049['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)[2][:,start:stop])

eeg_amplitude_data = np.zeros((50, number_of_animals))

for i in range(np.shape(eeg_amplitude_data)[0]):
    if i == 0:
        eeg_amplitude_data[i,0] = np.sum(sc016['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,1] = np.sum(sc024['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,2] = np.sum(sc025['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,3] = np.sum(sc026['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,4] = np.sum(sc028['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,5] = np.sum(sc030['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,6] = np.sum(sc031['eeg_spectrogram'][:,i], axis=0)
        #eeg_amplitude_data[i,7] = np.sum(sc032['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,7] = np.sum(sc033['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,8] = np.sum(sc034['eeg_spectrogram'][:,i], axis=0)
        #eeg_amplitude_data[i,10] = np.sum(sc036['eeg_spectrogram'][:,i], axis=0)
        #eeg_amplitude_data[i,11] = np.sum(sc037['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,9] = np.sum(sc038['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,10] = np.sum(sc039['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,11] = np.sum(sc040['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,12] = np.sum(sc041['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,13] = np.sum(sc044['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,14] = np.sum(sc045['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,15] = np.sum(sc046['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,16] = np.sum(sc047['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,17] = np.sum(sc048['eeg_spectrogram'][:,i], axis=0)
        eeg_amplitude_data[i,18] = np.sum(sc049['eeg_spectrogram'][:,i], axis=0)
        
        
    elif i == np.shape(eeg_amplitude_data)[0] - 1:
        start = int(i*2 - 1)
        eeg_amplitude_data[i,0] = np.sum(sc016['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,1] = np.sum(sc024['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,2] = np.sum(sc025['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,3] = np.sum(sc026['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,4] = np.sum(sc028['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,5] = np.sum(sc030['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,6] = np.sum(sc031['eeg_spectrogram'][:,start], axis=0)
        #eeg_amplitude_data[i,7] = np.sum(sc032['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,7] = np.sum(sc033['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,8] = np.sum(sc034['eeg_spectrogram'][:,start], axis=0)
        #eeg_amplitude_data[i,10] = np.sum(sc036['eeg_spectrogram'][:,start], axis=0)
        #eeg_amplitude_data[i,11] = np.sum(sc037['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,9] = np.sum(sc038['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,10] = np.sum(sc039['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,11] = np.sum(sc040['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,12] = np.sum(sc041['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,13] = np.sum(sc044['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,14] = np.sum(sc045['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,15] = np.sum(sc046['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,16] = np.sum(sc047['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,17] = np.sum(sc048['eeg_spectrogram'][:,start], axis=0)
        eeg_amplitude_data[i,18] = np.sum(sc049['eeg_spectrogram'][:,start], axis=0)
        
    else:
        start = int(i*2 - 1)
        stop = int(start+3)
        eeg_amplitude_data[i,0] = np.sum(np.mean(sc016['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,1] = np.sum(np.mean(sc024['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,2] = np.sum(np.mean(sc025['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,3] = np.sum(np.mean(sc026['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,4] = np.sum(np.mean(sc028['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,5] = np.sum(np.mean(sc030['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,6] = np.sum(np.mean(sc031['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        #eeg_amplitude_data[i,7] = np.sum(np.mean(sc032['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,7] = np.sum(np.mean(sc033['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,8] = np.sum(np.mean(sc034['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        #eeg_amplitude_data[i,10] = np.sum(np.mean(sc036['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        #eeg_amplitude_data[i,11] = np.sum(np.mean(sc037['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,9] = np.sum(np.mean(sc038['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,10] = np.sum(np.mean(sc039['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,11] = np.sum(np.mean(sc040['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,12] = np.sum(np.mean(sc041['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,13] = np.sum(np.mean(sc044['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,14] = np.sum(np.mean(sc045['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,15] = np.sum(np.mean(sc046['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,16] = np.sum(np.mean(sc047['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,17] = np.sum(np.mean(sc048['eeg_spectrogram'][:,start:stop], axis=1), axis=0)
        eeg_amplitude_data[i,18] = np.sum(np.mean(sc049['eeg_spectrogram'][:,start:stop], axis=1), axis=0)


#for i in range(np.shape(eeg_amplitude_data)[1]):
#    eeg_amplitude_data[:,i] = movingaverage(eeg_amplitude_data[:,i], 5)


scaled_eeg_amplitude_data = np.zeros((np.shape(eeg_amplitude_data)[0], np.shape(eeg_amplitude_data)[1]))
for i in range(np.shape(scaled_eeg_amplitude_data)[1]):
    scaled_eeg_amplitude_data[:, i] = (eeg_amplitude_data[:, i] - np.min(eeg_amplitude_data[:, i])) / (np.max(eeg_amplitude_data[:, i], axis=0) - np.min(eeg_amplitude_data[:, i])) 
    

sz_grab_data = np.zeros((120*130, number_of_animals))

sz_grab_data[:,0] = sc016['data'][1][int(sc016['seizure_start'] - 130*60):int(sc016['seizure_start'] + 130*60)]
sz_grab_data[:,1] = sc024['data'][1][int(sc024['seizure_start'] - 130*60):int(sc024['seizure_start'] + 130*60)]
sz_grab_data[:,2] = sc025['data'][1][int(sc025['seizure_start'] - 130*60):int(sc025['seizure_start'] + 130*60)]
sz_grab_data[:,3] = sc026['data'][1][int(sc026['seizure_start'] - 130*60):int(sc026['seizure_start'] + 130*60)]
sz_grab_data[:,4] = sc028['data'][1][int(sc028['seizure_start'] - 130*60):int(sc028['seizure_start'] + 130*60)]
sz_grab_data[:,5] = sc030['data'][1][int(sc030['seizure_start'] - 130*60):int(sc030['seizure_start'] + 130*60)]
sz_grab_data[:,6] = sc031['data']['analog_1'][int(sc031['seizure_start'] - 130*60):int(sc031['seizure_start'] + 130*60)]
#sz_grab_data[:,7] = sc032['data'][1][int(sc032['seizure_start'] - 130*60):int(sc032['seizure_start'] + 130*60)]
sz_grab_data[:,7] = sc033['data']['analog_1'][int(sc033['seizure_start'] - 130*60):int(sc033['seizure_start'] + 130*60)]
sz_grab_data[:,8] = sc034['data']['analog_1'][int(sc034['seizure_start'] - 130*60):int(sc034['seizure_start'] + 130*60)]
#sz_grab_data[:,10] = sc036['data']['analog_1'][int(sc036['seizure_start'] - 130*60):int(sc036['seizure_start'] + 130*60)]
#sz_grab_data[:,11] = sc037['data']['analog_1'][int(sc037['seizure_start'] - 130*60):int(sc037['seizure_start'] + 130*60)]
sz_grab_data[:,9] = sc038['data']['analog_1'][int(sc038['seizure_start'] - 130*60):int(sc038['seizure_start'] + 130*60)]
sz_grab_data[:,10] = sc039['data']['analog_1'][int(sc039['seizure_start'] - 130*60):int(sc039['seizure_start'] + 130*60)]
sz_grab_data[:,11] = sc040['data']['analog_1'][int(sc040['seizure_start'] - 130*60):int(sc040['seizure_start'] + 130*60)]
sz_grab_data[:,12] = sc041['data']['analog_1'][int(sc041['seizure_start'] - 130*60):int(sc041['seizure_start'] + 130*60)]
sz_grab_data[:,13] = sc044['data']['analog_1'][int(sc044['seizure_start'] - 130*60):int(sc044['seizure_start'] + 130*60)]
sz_grab_data[:,14] = sc045['data']['analog_1'][int(sc045['seizure_start'] - 130*60):int(sc045['seizure_start'] + 130*60)]
sz_grab_data[:,15] = sc046['data']['analog_1'][int(sc046['seizure_start'] - 130*60):int(sc046['seizure_start'] + 130*60)]
sz_grab_data[:,16] = sc047['data']['analog_1'][int(sc047['seizure_start'] - 130*60):int(sc047['seizure_start'] + 130*60)]
sz_grab_data[:,17] = sc048['data']['analog_1'][int(sc048['seizure_start'] - 130*60):int(sc048['seizure_start'] + 130*60)]
sz_grab_data[:,18] = sc049['data']['analog_1'][int(sc049['seizure_start'] - 130*60):int(sc049['seizure_start'] + 130*60)]

baseline = int(60 * 130)
sz_start = int(70 * 130)
base_means = np.mean(sz_grab_data[0:baseline], 0)

sz_grab_data[:,0] = sc016['data'][3][int(sc016['seizure_start'] - 130*60):int(sc016['seizure_start'] + 130*60)] + base_means[0]
sz_grab_data[:,1] = sc024['data'][3][int(sc024['seizure_start'] - 130*60):int(sc024['seizure_start'] + 130*60)] + base_means[1]
sz_grab_data[:,2] = sc025['data'][3][int(sc025['seizure_start'] - 130*60):int(sc025['seizure_start'] + 130*60)] + base_means[2]
sz_grab_data[:,3] = sc026['data'][3][int(sc026['seizure_start'] - 130*60):int(sc026['seizure_start'] + 130*60)] + base_means[3]
sz_grab_data[:,4] = sc028['data'][3][int(sc028['seizure_start'] - 130*60):int(sc028['seizure_start'] + 130*60)] + base_means[4]
sz_grab_data[:,5] = sc030['data'][3][int(sc030['seizure_start'] - 130*60):int(sc030['seizure_start'] + 130*60)] + base_means[5]
sz_grab_data[:,6] = sc031['data']['analog_1_filt'][int(sc031['seizure_start'] - 130*60):int(sc031['seizure_start'] + 130*60)] + base_means[6]
#sz_grab_data[:,7] = sc032['data'][1][int(sc032['seizure_start'] - 130*60):int(sc032['seizure_start'] + 130*60)] + base_means[7]
sz_grab_data[:,7] = sc033['data']['analog_1_filt'][int(sc033['seizure_start'] - 130*60):int(sc033['seizure_start'] + 130*60)] + base_means[7]
sz_grab_data[:,8] = sc034['data']['analog_1_filt'][int(sc034['seizure_start'] - 130*60):int(sc034['seizure_start'] + 130*60)] + base_means[8]
sz_grab_data[:,9] = sc038['data']['analog_1_filt'][int(sc038['seizure_start'] - 130*60):int(sc038['seizure_start'] + 130*60)] + base_means[9]
sz_grab_data[:,10] = sc039['data']['analog_1_filt'][int(sc039['seizure_start'] - 130*60):int(sc039['seizure_start'] + 130*60)] + base_means[10]
sz_grab_data[:,11] = sc040['data']['analog_1_filt'][int(sc040['seizure_start'] - 130*60):int(sc040['seizure_start'] + 130*60)] + base_means[11]
sz_grab_data[:,12] = sc041['data']['analog_1_filt'][int(sc041['seizure_start'] - 130*60):int(sc041['seizure_start'] + 130*60)] + base_means[12]
sz_grab_data[:,13] = sc044['data']['analog_1_filt'][int(sc044['seizure_start'] - 130*60):int(sc044['seizure_start'] + 130*60)] + base_means[13]
sz_grab_data[:,14] = sc045['data']['analog_1_filt'][int(sc045['seizure_start'] - 130*60):int(sc045['seizure_start'] + 130*60)] + base_means[14]
sz_grab_data[:,15] = sc046['data']['analog_1_filt'][int(sc046['seizure_start'] - 130*60):int(sc046['seizure_start'] + 130*60)] + base_means[15]
sz_grab_data[:,16] = sc047['data']['analog_1_filt'][int(sc047['seizure_start'] - 130*60):int(sc047['seizure_start'] + 130*60)] + base_means[16]
sz_grab_data[:,17] = sc048['data']['analog_1_filt'][int(sc048['seizure_start'] - 130*60):int(sc048['seizure_start'] + 130*60)] + base_means[17]
sz_grab_data[:,18] = sc049['data']['analog_1_filt'][int(sc049['seizure_start'] - 130*60):int(sc049['seizure_start'] + 130*60)] + base_means[18]

sz_grab_data_normalized = np.zeros((120*130, number_of_animals))
sz_grab_data_normalized_1s = np.zeros((120, number_of_animals))

for i in range(np.shape(sz_grab_data_normalized)[1]):
    baseline = np.mean(sz_grab_data[0:int(50*130),i], axis = 0)
    sz_grab_data_normalized[:,i] = ((sz_grab_data[:,i] - baseline) / baseline) * 100
    
for i in range(np.shape(sz_grab_data_normalized_1s)[0]):
    start = int(i*130)
    stop = int(start + 130)
    sz_grab_data_normalized_1s[i,:] = np.mean(sz_grab_data_normalized[start:stop,:], axis=0) 


   
sz_grab_data_normalized_1s = sz_grab_data_normalized_1s[70:,:]

for i in range(np.shape(sz_grab_data_normalized_1s)[1]):
    sz_grab_data_normalized_1s[:,i] = (sz_grab_data_normalized_1s[:,i] - np.min(sz_grab_data_normalized_1s[:,i])) / (np.max(sz_grab_data_normalized_1s[:,i]) - np.min(sz_grab_data_normalized_1s[:,i]))



        
correlations = np.zeros((50,19))
test = correlate(sz_grab_data_normalized_1s[:,0], scaled_eeg_amplitude_data[:,0], mode='same')
for i in range(np.shape(sz_grab_data_normalized_1s)[1]):
    a = (scaled_eeg_amplitude_data[:,i] - np.mean(scaled_eeg_amplitude_data[:,i], axis=0)) / (np.std(scaled_eeg_amplitude_data[:,i],axis = 0) * len(scaled_eeg_amplitude_data[:,i]))
    b = (sz_grab_data_normalized_1s[:,i] - np.mean(sz_grab_data_normalized_1s[:,i], axis=0)) / (np.std(sz_grab_data_normalized_1s[:,i]))
    #plt.plot(b)
    correlations[:,i] = correlate(a, b, mode='same')
    


x = correlation_lags(50, 50, mode='same')

lags = np.zeros((19))
peak_r = np.zeros((19))
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
ax[1].set_yticks([2, 0, -2, -4, -6, -8, -10, -12])
ax[1].set_yticklabels([2, 0, -2, -4, -6, -8, -10, -12], fontsize=12)
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



