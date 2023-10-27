# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 12:15:00 2022

@author: llarsen
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(1, r'E:\OneDrive - UGent')

from f_cross_corr_function import eegpower_spikerate_crosscorr
from scipy.signal import correlation_lags
from sklearn.cluster import DBSCAN
from matplotlib.colors import ListedColormap
from scipy.stats import ttest_ind


### GENERAL PARAMS
eeg_downsample_rate = 2000
unit_totalpowerate = 30000

savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure06\output\Fig6C.png'


### LL001
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL001\Seizure1.mat'
eeg_totalpowerate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
LL001_spike_templates = [5, 6, 10, 48, 59, 91, 92, 101, 110, 114, 118, 125, 149, 164, 165, 176, 178, 180, 197]
seizure_start = 3782.72

LL001_all_freqs, LL001_totalpower = eegpower_spikerate_crosscorr(eeg_path, eeg_totalpowerate,
                             eeg_downsample_rate,
                             eeg_stimstart,
                             phy_path,
                             unit_totalpowerate,
                             LL001_spike_templates,
                             seizure_start
                             )

### LL004
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL004\Seizure1.mat'
eeg_totalpowerate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL004-final.GUI'
LL004_spike_templates = [18, 30, 40, 51, 52]
seizure_start = 1161.98

LL004_all_freqs, LL004_totalpower = eegpower_spikerate_crosscorr(eeg_path, eeg_totalpowerate,
                             eeg_downsample_rate,
                             eeg_stimstart,
                             phy_path,
                             unit_totalpowerate,
                             LL004_spike_templates,
                             seizure_start
                             )


### LL005
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL005\Seizure1.mat'
eeg_totalpowerate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'
LL005_spike_templates = [17, 20, 26, 36, 38, 40, 47, 53, 61, 63, 77, 84, 95, 100, 112, 120, 124, 128, 136, 146]
seizure_start = 770.52

LL005_all_freqs, LL005_totalpower = eegpower_spikerate_crosscorr(eeg_path, eeg_totalpowerate,
                             eeg_downsample_rate,
                             eeg_stimstart,
                             phy_path,
                             unit_totalpowerate,
                             LL005_spike_templates,
                             seizure_start
                             )

### LL007
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL007\Seizure0.mat'
eeg_totalpowerate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'
LL007_spike_templates = [81, 107, 112, 135, 157, 159, 181]
seizure_start = 1193.42

LL007_all_freqs, LL007_totalpower = eegpower_spikerate_crosscorr(eeg_path, eeg_totalpowerate,
                             eeg_downsample_rate,
                             eeg_stimstart,
                             phy_path,
                             unit_totalpowerate,
                             LL007_spike_templates,
                             seizure_start
                             )

### LL009
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL009\5_Seizure1_nolight201019A0000.mat'
eeg_totalpowerate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL009-final.GUI'
LL009_spike_templates = [47]
seizure_start = 4408.35

LL009_all_freqs, LL009_totalpower = eegpower_spikerate_crosscorr(eeg_path, eeg_totalpowerate,
                             eeg_downsample_rate,
                             eeg_stimstart,
                             phy_path,
                             unit_totalpowerate,
                             LL009_spike_templates,
                             seizure_start
                             )


### LL011contra
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL011\Seizures_contra\5_seizure1_nolight201026A0018.mat'
eeg_totalpowerate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL011_contra-final.GUI'
LL0011contra_spike_templates = [15, 33, 95, 104, 107, 140]
seizure_start = 542.118

LL0011contra_all_freqs, LL0011contra_totalpower = eegpower_spikerate_crosscorr(eeg_path, eeg_totalpowerate,
                             eeg_downsample_rate,
                             eeg_stimstart,
                             phy_path,
                             unit_totalpowerate,
                             LL0011contra_spike_templates,
                             seizure_start
                             )


### LL012ipsi
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL012\Seizures - ipsi\2_seizure1_nolight201028A0012.mat'
eeg_totalpowerate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL012_ipsi-final.GUI'
LL0012ipsi_spike_templates = [137] ## 112 bugs out somehow
seizure_start = 729.362

LL0012ipsi_all_freqs, LL0012ipsi_totalpower = eegpower_spikerate_crosscorr(eeg_path, eeg_totalpowerate,
                             eeg_downsample_rate,
                             eeg_stimstart,
                             phy_path,
                             unit_totalpowerate,
                             LL0012ipsi_spike_templates,
                             seizure_start
                             )

### LL012contra
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL012\Seizures - contra\5_cl_seizure1_nolight201028A0020.mat'
eeg_totalpowerate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL012_contra-final.GUI'
LL0012contra_spike_templates = [12, 65, 156, 166, 170, 200, 204, 218, 317, 318, 329, 348, 354]
seizure_start = 124.806

LL0012contra_all_freqs, LL0012contra_totalpower = eegpower_spikerate_crosscorr(eeg_path, eeg_totalpowerate,
                             eeg_downsample_rate,
                             eeg_stimstart,
                             phy_path,
                             unit_totalpowerate,
                             LL0012contra_spike_templates,
                             seizure_start
                             )


### LL017
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL017\05_Seizure5201208A0017.mat'
eeg_totalpowerate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL017-final.GUI'
LL017_spike_templates = [276, 293, 300, 317]
seizure_start = 3621.83

LL017_all_freqs, LL017_totalpower = eegpower_spikerate_crosscorr(eeg_path, eeg_totalpowerate,
                             eeg_downsample_rate,
                             eeg_stimstart,
                             phy_path,
                             unit_totalpowerate,
                             LL017_spike_templates,
                             seizure_start
                             )

### LL023
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL023\4_Seizure300_2201218A0004.mat'
eeg_totalpowerate = 10000
eeg_stimstart = 60
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL023-final.GUI'
LL023_spike_templates = [4, 5, 9, 12, 33, 48, 49]
seizure_start = 2061.56

LL023_all_freqs, LL023_totalpower = eegpower_spikerate_crosscorr(eeg_path, eeg_totalpowerate,
                             eeg_downsample_rate,
                             eeg_stimstart,
                             phy_path,
                             unit_totalpowerate,
                             LL023_spike_templates,
                             seizure_start
                             )

all_corr_total = np.column_stack((LL001_totalpower, LL004_totalpower, LL005_totalpower, LL007_totalpower, LL009_totalpower, LL0011contra_totalpower, LL0012ipsi_totalpower, LL0012contra_totalpower, LL017_totalpower, LL023_totalpower))

indices = np.zeros(np.shape(all_corr_total)[1])
for i in range(np.shape(all_corr_total)[1]):
    indices[i] = int(np.where(np.max(np.abs(all_corr_total[:,i])) == np.abs(all_corr_total[:,i]))[0])

lags = correlation_lags(50,50, mode='same')

peak_corrs = np.zeros(np.shape(all_corr_total)[1])
lag_peak_corrs = np.zeros(np.shape(all_corr_total)[1])

for i in range(np.shape(all_corr_total)[1]):
    peak_corrs[i] = all_corr_total[int(indices[i]), i]
    lag_peak_corrs[i] = lags[int(indices[i])]



def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))


data = np.column_stack((NormalizeData(lag_peak_corrs), NormalizeData(peak_corrs)))

clustering = DBSCAN(eps=0.15, min_samples=8).fit(data)
data = np.column_stack((data, clustering.labels_))

data_real = np.column_stack((lag_peak_corrs, peak_corrs, data[:,2]))

inhibited = 0
excited = 0
decorrelated = 0

for i in range(len(data[:,2])):
    if data[i,2] == -1:
        decorrelated = decorrelated + 1
    elif data[i,2] == 0:
        excited = excited + 1
    elif data[i,2] == 1:
        inhibited = inhibited + 1
    
inhibited = np.zeros((inhibited, 2))
excited = np.zeros((excited, 2))
decorrelated = np.zeros((decorrelated, 2))

counter0 = 0
counter1 = 0
counter2= 0
for i in range(len(data[:,2])):
    if data[i,2] == -1:
        decorrelated[counter0,0] = lag_peak_corrs[i]
        decorrelated[counter0,1]  = peak_corrs[i]   
        counter0 = counter0 + 1
    elif data[i,2] == 0:
        excited[counter1,0]  = lag_peak_corrs[i]
        excited[counter1,1]  = peak_corrs[i]        
        counter1 = counter1 + 1
    elif data[i,2] == 1:
        inhibited[counter2,0]  = lag_peak_corrs[i]
        inhibited[counter2,1]  = peak_corrs[i]
        counter2 = counter2 + 1

x_excited = np.linspace(1.3+0.5, 1.6+0.5, len(excited[:,0]))
x_inhibited = np.linspace(0.3, 0.6, len(inhibited[:,0]))

p_lag = ttest_ind(inhibited[:,0], excited[:,0])
p_xcorr = ttest_ind(inhibited[:,1], excited[:,1])
means_inhibited = np.mean(inhibited, axis=0)
means_excited = np.mean(excited, axis=0)
std_inhibited = np.std(inhibited, axis=0)
std_excited = np.std(excited, axis=0)


col = ListedColormap(['g', 'r', 'b'])

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(3, figsize=(4,6))

ax[0].scatter(lag_peak_corrs, peak_corrs, c=data[:,2], cmap=col)

ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_linewidth(2)
ax[0].spines['left'].set_linewidth(2)

ax[0].set_xticklabels([],[])
ax[0].tick_params(axis='x', which='both', bottom=False)

ax[0].set_xticks([-50, -25, 0, 25, 50])
ax[0].set_xticklabels([-50, -25, 0, 25, 50], fontsize = 12)


ax[0].set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
ax[0].set_yticklabels([-1.0, -0.5, 0, 0.5, 1.0], fontsize = 12)

ax[0].set_xlabel('Lag (s)',fontsize=14, fontweight='bold')
ax[0].set_ylabel('xcorr',fontsize=14, fontweight='bold')

ax[0].set_ylim(-1,1)
ax[0].set_xlim(-25,25)

ax[1].boxplot(inhibited[:,0], positions = [0], showmeans=True, meanline=True, widths = 0.2, 
           boxprops= dict(linewidth=2.0, color='black'), 
           whiskerprops=dict(linestyle='-',linewidth=2.0, color='black'),
           capprops=dict(linestyle='-',linewidth=2.0, color='black'), 
           medianprops=dict(linestyle='--',linewidth=2.0, color='black'),
           meanprops=dict(linestyle='-',linewidth=2.0, color='black'))
ax[1].boxplot(excited[:,0], positions = [1+0.5], showmeans=True, meanline=True, widths = 0.2, 
           boxprops= dict(linewidth=2.0, color='black'), 
           whiskerprops=dict(linestyle='-',linewidth=2.0, color='black'),
           capprops=dict(linestyle='-',linewidth=2.0, color='black'), 
           medianprops=dict(linestyle='--',linewidth=2.0, color='black'),
           meanprops=dict(linestyle='-',linewidth=2.0, color='black'))

ax[1].scatter(x_inhibited, inhibited[:,0], c='b')
ax[1].scatter(x_excited, excited[:,0], c='r')
ax[1].hlines(0, 0,2.1, linestyle='--', linewidth=1.5, color='k',zorder=3)

ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_linewidth(2)

ax[1].set_xticks([0.25, 1.5+0.25])
ax[1].set_xticklabels(['Inhibited', 'Excited'], fontsize = 14)
ax[1].set_yticks([-10,-5,0])
ax[1].set_yticklabels([-10,-5,0], fontsize = 12)
ax[1].tick_params(axis='x', which='both', bottom=False)
ax[1].set_ylabel('Lag (s)',fontsize=14, fontweight='bold')
ax[1].set_ylim(-10, 2)
#ax.boxplot(decorrelated[:,0], positions = [2])

ax[2].boxplot(np.sqrt(inhibited[:,1]**2), positions = [0], showmeans=True, meanline=True, widths = 0.2, 
           boxprops= dict(linewidth=2.0, color='black'), 
           whiskerprops=dict(linestyle='-',linewidth=2.0, color='black'),
           capprops=dict(linestyle='-',linewidth=2.0, color='black'), 
           medianprops=dict(linestyle='--',linewidth=2.0, color='black'),
           meanprops=dict(linestyle='-',linewidth=2.0, color='black'))
ax[2].boxplot(np.sqrt(excited[:,1]**2), positions = [1+0.5], showmeans=True, meanline=True, widths = 0.2, 
           boxprops= dict(linewidth=2.0, color='black'), 
           whiskerprops=dict(linestyle='-',linewidth=2.0, color='black'),
           capprops=dict(linestyle='-',linewidth=2.0, color='black'), 
           medianprops=dict(linestyle='--',linewidth=2.0, color='black'),
           meanprops=dict(linestyle='-',linewidth=2.0, color='black'))

ax[2].scatter(x_inhibited, np.abs(inhibited[:,1]), c='b')
ax[2].scatter(x_excited, np.abs(excited[:,1]), c='r')

ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[2].spines['left'].set_linewidth(2)   

ax[2].set_xticks([0.25, 1.5+0.25])
ax[2].set_xticklabels(['Inhibited', 'Excited'], fontsize = 14)
ax[2].tick_params(axis='x', which='both', bottom=False)
ax[2].set_yticks([0.4, 0.6, 0.8, 1.0])
ax[2].set_yticklabels([0.4, 0.6, 0.8, 1.0], fontsize = 12)
ax[2].set_ylabel('Absolute xcorr',fontsize=14, fontweight='bold')


plt.tight_layout()
plt.savefig(os.path.abspath(savepath), dpi=600)

### compute correlograms

#for i in range(np.shape(all_corr_total)[1]):
#    if data[i,2] == -1:
#        plt.plot(all_corr_total[:,i], c='k')
#    elif data[i,2] == 0:
#        plt.plot(all_corr_total[:,i], c='b')
#    elif data[i,2] == 1:
#        plt.plot(all_corr_total[:,i], c='r')




