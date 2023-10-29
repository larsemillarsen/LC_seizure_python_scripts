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
#import pandas as pd
#from scipy.signal import butter, filtfilt
from scipy.stats import ttest_1samp
import scipy

directory = r'E:\Manuscript_analysis_files\data\GRABne'
savepath_grp = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure08\output\Fig8D_1_2.png'
savepath1 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure08\output\Fig8A_C.png'

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

number_of_animals = 19

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

change_percent = np.zeros(number_of_animals)
change_z = np.zeros(number_of_animals)
for i in range(number_of_animals):
    mean_base = np.mean(sz_grab_data[0:baseline,i])
    std_base = np.std(sz_grab_data[0:baseline,i])
    peak = np.max(sz_grab_data[sz_start:,i])
    z = (sz_grab_data[:,i] - mean_base) / std_base
    #plt.plot(z)
    change_z[i] = np.max(z[sz_start:])
    change_percent[i] = ((peak - mean_base) / mean_base) * 100
  

p_change = ttest_1samp(change_percent, 0)

p_values = scipy.stats.norm.sf(abs(change_z))*2




### GRP PLOTS

x = np.linspace(0.7, 0.75, number_of_animals)
meanprops = {'linestyle': 'solid',
            'color':'black',
            'linewidth':2}

medianprops  = {'linestyle': '--',
            'color':'black',
            'linewidth':2}

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6,3))
ax[0].scatter(x, change_percent)
ax[0].boxplot(change_percent, showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)
ax[0].set_xlim(0.6 ,1.4)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_linewidth(2)
ax[0].set_ylabel('$GRAB_{NE2m}$ $\Delta$F/$F_{0}$',fontsize=16, fontweight='bold')
ax[0].set_xticks([])
ax[0].set_xticklabels([])
ax[0].set_yticks([0,20,40,60,80,100, 120])
ax[0].set_yticklabels([0,20,40,60,80,100, 120], fontsize=12)
ax[0].text(0.83, 115, '*', fontsize = 25, fontweight='bold')
ax[1].scatter(x, change_z)
ax[1].boxplot(change_z, showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)
ax[1].set_xlim(0.6 ,1.4)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_linewidth(2)
ax[1].set_ylabel('$GRAB_{NE2m}$ Z-score',fontsize=16, fontweight='bold')
ax[1].hlines(1.96, 0.6, 1.4, colors='r', linestyles='--')
ax[1].set_xticks([])
ax[1].set_xticklabels([])
ax[1].set_yticks([0, 100, 200, 300, 400, 500])
ax[1].set_yticklabels([0, 100, 200, 300, 400, 500], fontsize=12)
plt.tight_layout()
plt.savefig(savepath_grp, dpi=600)

#### EXAMPLE PLOT

sc025['eeg_path'] = r'E:\Manuscript_analysis_files\data\lfp_data\SC025\SC025_Seizure1220615A0006.mat'
sc025['eeg_data'] = load_neuronfile(sc025['eeg_path'], srate = 10000, channel=1, gain=100, inputrange=20)

sc025['eeg_spectrogram'] = spectrogram(sc025['eeg_data'][1], srate=10000, windowlength=1, overlap=0.5,highfreq=100)


example_FoverF0 = (sz_grab_data[:,2] - np.mean(sz_grab_data[0:baseline,2])) / np.mean(sz_grab_data[0:baseline,2]) * 100
example_t = np.arange(0, 120, 1/130)

example_z = (sz_grab_data[:,1] - np.mean(sz_grab_data[0:baseline,1])) / np.std(sz_grab_data[0:baseline,1])
z_bins = np.zeros(int(120/5))
for i in range(len(z_bins)):
    start = i*(5*130)
    stop = start + (5*130) - 1
    z_bins[i] = np.mean(example_z[start:stop])
    
z_time = np.arange(0, 120, 5) + 5/2

t_index = np.where(example_FoverF0 == np.max(example_FoverF0))[0][0]
t = example_t[t_index]

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(3, figsize=(9,6))

ax[0].plot(sc025['eeg_data'][0]-50,sc025['eeg_data'][1],'b', label='Baseline (pre-stimulation)', linewidth=0.3)
ax[0].vlines(0, -15, 15, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[0].vlines(10, -15, 15, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[1].contourf(sc025['eeg_spectrogram'][0]-50, sc025['eeg_spectrogram'][1],np.log(sc025['eeg_spectrogram'][2]),300,cmap='jet')
ax[1].vlines(0, 0, 100, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[1].vlines(10, 0, 100, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[2].plot(example_t-60,example_FoverF0,'g')
ax[2].scatter(t-60, np.max(example_FoverF0), marker='*', color='k', zorder=3, s=100)   
ax[2].vlines(0, np.min(example_FoverF0), np.max(example_FoverF0), colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[2].vlines(10, np.min(example_FoverF0), np.max(example_FoverF0), colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_linewidth(2)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_linewidth(2)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_linewidth(2)
ax[2].spines['left'].set_linewidth(2)
ax[0].set_ylim([-15,10])
ax[0].set_xlim([-10,50])
ax[1].set_xlim([-10,50])
ax[2].set_xlim([-10,50])
ax[0].set_xticklabels([],[])
ax[0].tick_params(axis='x', which='both', bottom=False)
ax[0].tick_params(axis='x', which='both', bottom=False)
ax[1].set_xticklabels([],[])
ax[1].tick_params(axis='x', which='both', bottom=False)
ax[0].yaxis.set_tick_params(width=2)
ax[1].yaxis.set_tick_params(width=2)
ax[2].yaxis.set_tick_params(width=2)
ax[0].tick_params(axis="y", labelsize=14)
ax[1].tick_params(axis="y", labelsize=14)
ax[2].tick_params(axis="y", labelsize=14)
ax[0].set_ylabel('Hip LFP\nAmplitude (mV)', fontsize=14, fontweight='bold')
ax[1].set_ylabel('Hip LFP\nFrequency (Hz)', fontsize=14, fontweight='bold')
ax[2].set_ylabel('$GRAB_{NE2m}$\n$\Delta$F/$F_{0}$', fontsize=14, fontweight='bold')
ax[2].set_xlabel('Time (s)', fontsize=16, fontweight='bold')
ax[2].text(25, 100, '*Z-score = ' + str(np.around(change_z[2], 1)), fontsize = 14, fontweight='bold')
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.tight_layout()
plt.savefig(savepath1, dpi=600)
#plt.close()


mean_change = np.mean(change_percent)
std_change = np.std(change_percent)


print(mean_change)
print(std_change)