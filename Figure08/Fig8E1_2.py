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
from vies.parse.pyphotometry import import_ppd
from vies.parse.neuron import load_neuronfile
from vies.lfp.lfp_analysis import spectrogram
from scipy.stats import ttest_1samp
import scipy

directory = r'E:\Manuscript_analysis_files\data\GRABne'
savepath_grp = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure08\output\Fig8E_1_2.png'

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

number_of_animals = 10

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
ax[0].set_yticks([0,20,40,60])
ax[0].set_yticklabels([0,20,40,60], fontsize=12)
ax[0].text(0.83, 60, '*', fontsize = 25, fontweight='bold')

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
ax[1].set_yticks([0, 100, 200, 300, 400])
ax[1].set_yticklabels([0, 100, 200, 300, 400], fontsize=12)
plt.tight_layout()
plt.savefig(savepath_grp, dpi=600)




