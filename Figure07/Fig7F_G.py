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

from f_lfpspike_unitcoupling import popspike_coupling
from scipy.stats import ks_2samp

### GENERAL PARAMS
unit_srate = 30000

### LL001
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL001\Seizure1.mat'
eeg_srate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
LL001_spike_templates = [5, 6, 48, 59, 91, 114, 125, 149, 164]
seizure_start = 3782.72

pLL001, hLL001 = popspike_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL001_spike_templates,
                             seizure_start
                             )

### LL004
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL004\Seizure1.mat'
eeg_srate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL004-final.GUI'
LL004_spike_templates = [51, 52]
seizure_start = 1161.98

pLL004, hLL004 = popspike_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL004_spike_templates,
                             seizure_start, 
                             window=1.0
                             )


### LL005
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL005\Seizure1.mat'
eeg_srate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'
LL005_spike_templates = [95, 120, 136]
seizure_start = 770.52

pLL005, hLL005 = popspike_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL005_spike_templates,
                             seizure_start, 
                             window=1.0
                             )

### LL007
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL007\Seizure0.mat'
eeg_srate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'
LL007_spike_templates = [81, 157]
seizure_start = 1193.42

pLL007, hLL007 = popspike_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL007_spike_templates,
                             seizure_start, 
                             window=1.0
                             )


### LL011contra
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL011\Seizures_contra\5_seizure1_nolight201026A0018.mat'
eeg_srate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL011_contra-final.GUI'
LL0011contra_spike_templates = [107]
seizure_start = 542.118

pLL011contra, hLL011contra = popspike_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL0011contra_spike_templates,
                             seizure_start, 
                             window=1.0
                             )


### LL012ipsi
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL012\Seizures - ipsi\2_seizure1_nolight201028A0012.mat'
eeg_srate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL012_ipsi-final.GUI'
LL0012ipsi_spike_templates = [122, 137] ## 112 bugs out somehow
seizure_start = 729.362

pLL012ipsi, hLL012ipsi = popspike_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL0012ipsi_spike_templates,
                             seizure_start, 
                             window=1.0
                             )

### LL012contra
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL012\Seizures - contra\5_cl_seizure1_nolight201028A0020.mat'
eeg_srate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL012_contra-final.GUI'
LL0012contra_spike_templates = [318, 329]
seizure_start = 124.806

pLL012contra, hLL012contra = popspike_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL0012contra_spike_templates,
                             seizure_start, 
                             window=1.0
                             )


### LL023
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL023\4_Seizure300_2201218A0004.mat'
eeg_srate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL023-final.GUI'
LL023_spike_templates = [4, 9, 12, 33, 48, 49]
seizure_start = 2061.56

pLL023, hLL023 = popspike_coupling(eeg_path, eeg_srate,
                             eeg_stimstart,
                             phy_path,
                             unit_srate,
                             LL023_spike_templates,
                             seizure_start, 
                             window=1.0
                             )



all_neurons = np.concatenate((pLL001, pLL004, pLL005, pLL007, pLL011contra, pLL012ipsi, pLL012contra, pLL023))

all_neurons_corrected = multipletests(all_neurons, alpha=0.05, method='bonferroni')[1]

window = 1 # in seconds
bins = 10
bin_width = window/bins


histograms = np.zeros((bins, len(all_neurons)))

counter = 0
for i in range(len(hLL001)):
    histograms[:, counter] = np.histogram(hLL001[i], bins=bins, range=(-window/2, window/2))[0]
    
    counter = counter + 1

for i in range(len(hLL004)):
    histograms[:, counter] = np.histogram(hLL004[i], bins=bins, range=(-window/2, window/2))[0]
    
    counter = counter + 1

for i in range(len(hLL005)):
    histograms[:, counter] = np.histogram(hLL005[i], bins=bins, range=(-window/2, window/2))[0]
    
    counter = counter + 1

for i in range(len(hLL007)):
    histograms[:, counter] = np.histogram(hLL007[i], bins=bins, range=(-window/2, window/2))[0]
    
    counter = counter + 1   

for i in range(len(hLL011contra)):
    histograms[:, counter] = np.histogram(hLL011contra[i], bins=bins, range=(-window/2, window/2))[0]
    
    counter = counter + 1

for i in range(len(hLL012contra)):
    histograms[:, counter] = np.histogram(hLL012contra[i], bins=bins, range=(-window/2, window/2))[0]
    
    counter = counter + 1   

for i in range(len(hLL012ipsi)):
    histograms[:, counter] = np.histogram(hLL012ipsi[i], bins=bins, range=(-window/2, window/2))[0]
    
    counter = counter + 1
    
for i in range(len(hLL023)):
    histograms[:, counter] = np.histogram(hLL023[i], bins=bins, range=(-window/2, window/2))[0]
    
    counter = counter + 1

x = np.histogram(hLL023[0], bins=bins, range=(-window/2, window/2))[1]
x = x[:-1] + bin_width/2





clevels = np.linspace(0.3, 0.7, len(all_neurons)) 

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

fig, ax = plt.subplots(1, figsize=(3,5))
  
ax.scatter(clevels, all_neurons)

ax.hlines(0.05, 0, 1.0, colors='r', linestyles='--')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_linewidth(2)

ax.set_ylim(np.min(all_neurons)-(np.min(all_neurons)/2), 2)
ax.set_yscale('log')
ax.set_xticks([])
ax.set_xticklabels([]) 
ax.tick_params(axis='y') 

ax.set_xlim(-0.4, 1.6)

ax.text(-0.6, 3, 'n = ' + str(np.shape(all_neurons)[0]) + ' neurons', fontsize = 14)
ax.set_ylabel('p-value',fontsize=16, fontweight='bold')
ax.tick_params(axis='x', which='both', bottom=False)


plt.tight_layout()
savepath2 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure07\output\Fig7G.png'
plt.savefig(os.path.abspath(savepath2), dpi=600)




uni_dist = np.linspace(np.min(hLL007[0]), np.max(hLL007[0]), len(hLL007[0]))
p = ks_2samp(hLL007[0], uni_dist)[1]
            
x = np.arange(-500,500, 100) + 50
y = np.histogram(hLL007[0], bins = 10)[0] / np.sum(np.histogram(hLL007[0], bins = 10)[0])
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(4,4))
ax.bar(x, y, width = 100*0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.text(250, np.max(y), 'p = ' + str(np.around(p, 3)) + '\n' +
    'p* = ' + str(np.around(all_neurons_corrected[14], 3)), fontsize = 12)
ax.set_yticks([0, 0.05, 0.1, 0.15, 0.2])
ax.set_yticklabels([0, 0.05, 0.10, 0.15, 0.20]) 
#ax.set_xticks([110, 120, 130, 140, 150])
#ax.set_xticklabels([-10, 0, 10, 20, 30]) 
ax.set_xlabel('Time (ms)',fontsize=14, fontweight='bold')
ax.set_ylabel('Spike probability',fontsize=14, fontweight='bold')
savepath2 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure07\output\Fig7F.png'
plt.savefig(os.path.abspath(savepath2), dpi=600)
