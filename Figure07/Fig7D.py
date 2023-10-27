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

from statsmodels.stats.multitest import multipletests

sys.path.insert(1, r'E:\Manuscript_analysis_files')

from f_preictal import spike_lfpphase_coupling

### GENERAL PARAMS
unit_srate = 30000


### LL001
eeg_path = r'E:\Data\LCSeizureData\LCSeizureNIDAQ\13_LL001\Seizure\Seizure1.mat'
eeg_srate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
LL001_spike_templates = [5, 6, 48, 59, 91, 114, 125, 149, 164]
seizure_start = 3782.72

pLL001, meanLL001 = spike_lfpphase_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL001_spike_templates,
                             seizure_start
                             )

### LL004
eeg_path = r'E:\Data\LCSeizureData\LCSeizureNIDAQ\11_LL004\Seizure1.mat'
eeg_srate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL004-final.GUI'
LL004_spike_templates = [51, 52]
seizure_start = 1161.98

pLL004, meanLL004 = spike_lfpphase_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL004_spike_templates,
                             seizure_start
                             )


### LL005
eeg_path = r'E:\Data\LCSeizureData\LCSeizureNIDAQ\14_LL005\Seizure1.mat'
eeg_srate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'
LL005_spike_templates = [95, 120, 136]
seizure_start = 770.52

pLL005, meanLL005 = spike_lfpphase_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL005_spike_templates,
                             seizure_start
                             )

### LL007
eeg_path = r'E:\Data\LCSeizureData\LCSeizureNIDAQ\16_LL007\Seizure\Seizure0.mat'
eeg_srate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'
LL007_spike_templates = [81, 157]
seizure_start = 1193.42

pLL007, meanLL007 = spike_lfpphase_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL007_spike_templates,
                             seizure_start
                             )


### LL011contra
eeg_path = r'E:\Data\LCSeizureData\LCSeizureNIDAQ\LL011\Seizures_contra\5_seizure1_nolight201026A0018.mat'
eeg_srate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL011_contra-final.GUI'
LL0011contra_spike_templates = [107]
seizure_start = 542.118

pLL011contra, meanLL011contra = spike_lfpphase_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL0011contra_spike_templates,
                             seizure_start
                             )


### LL012ipsi
eeg_path = r'E:\Data\LCSeizureData\LCSeizureNIDAQ\LL012\Seizures - ipsi\2_seizure1_nolight201028A0012.mat'
eeg_srate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL012_ipsi-final.GUI'
LL0012ipsi_spike_templates = [122, 137] ## 112 bugs out somehow
seizure_start = 729.362

pLL012ipsi, meanLL012ipsi = spike_lfpphase_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL0012ipsi_spike_templates,
                             seizure_start
                             )

### LL012contra
eeg_path = r'E:\Data\LCSeizureData\LCSeizureNIDAQ\LL012\5_cl_seizure1_nolight201028A0020.mat'
eeg_srate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL012_contra-final.GUI'
LL0012contra_spike_templates = [318, 329]
seizure_start = 124.806

pLL012contra, meanLL012contra = spike_lfpphase_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL0012contra_spike_templates,
                             seizure_start
                             )


### LL023
eeg_path = r'E:\Data\LCSeizureData\LCSeizureNIDAQ\LL023\Seizures\4_Seizure300_2201218A0004.mat'
eeg_srate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL023-final.GUI'
LL023_spike_templates = [4, 9, 12, 33, 48, 49]
seizure_start = 2061.56

pLL023, meanLL023 = spike_lfpphase_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL023_spike_templates,
                             seizure_start
                             )

all_neurons = np.column_stack((pLL001, pLL004, pLL005, pLL007, pLL011contra, pLL012ipsi, pLL012contra, pLL023))
all_neurons_mean = np.column_stack((meanLL001, pLL004, meanLL005, meanLL007, meanLL011contra, meanLL012ipsi, meanLL012contra, meanLL023))

all_neurons_corrected = all_neurons

for i in range(np.shape(all_neurons)[0]):
    all_neurons_corrected[i,:] = multipletests(all_neurons[i,:], alpha=0.05, method='bonferroni')[1]
#all_neurons = all_neurons * np.shape(all_neurons)[0]


plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
### boxplot with dots
ngroup = np.shape(all_neurons)[1]

clevels = np.zeros((np.shape(all_neurons)[0], np.shape(all_neurons)[1]))
for i in range(np.shape(all_neurons)[0]):
    clevels[i,:] = np.linspace(i+0.3, i+0.7, ngroup)

fig, ax = plt.subplots(1, figsize=(5,5))

meanprops = {'linestyle': 'solid',
            'color':'black',
            'linewidth':2}

medianprops  = {'linestyle': '--',
            'color':'black',
            'linewidth':2}

#ax.boxplot(all_neurons[0,:], showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)
for i in range(np.shape(all_neurons)[0]):
    ax.scatter(clevels[i,:], all_neurons[i,:])

ax.hlines(0.05,np.min(clevels), np.max(clevels), colors='r', linestyles='--')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_linewidth(2)

ax.set_ylim(np.min(all_neurons)-(np.min(all_neurons)/2), 2)
ax.set_yscale('log')
#ax.set_scale
ax.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
ax.set_xticklabels(['1-2 Hz', '2-4 Hz', '4-8 Hz', '8-16 Hz', '16-32 Hz', '32-64 Hz'], fontsize=12, fontweight='bold', rotation=75) 
ax.tick_params(axis='y') 


ax.set_xlim(-0.4, 6.6)
ax.set_ylim(1e-6, 1.1)

ax.text(3.0, 3, 'n = ' + str(np.shape(all_neurons)[1]) + ' neurons', fontsize = 14)
#ax.set_xlabel('Hip LFP Frequency (Hz)',fontsize=16, fontweight='bold')
ax.set_ylabel('p-value',fontsize=16, fontweight='bold')
#ax.set_xticklabels([],[])
ax.tick_params(axis='x', which='both', bottom=False)



plt.tight_layout()
savepath2 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure07\output\Fig7D.png'
plt.savefig(os.path.abspath(savepath2), dpi=600)






