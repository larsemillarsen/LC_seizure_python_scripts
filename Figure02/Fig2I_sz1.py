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
import pandas as pd

sys.path.insert(1, r'E:\OneDrive - UGent\python_functions')

from vies.parse.phy import phy_data

path = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\Fig2H\output_data'
savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\output\fig2I_sz1.png'

def spike_train_extract(phy_path, spike_templates, seizure_start):
    
    unit_srate = 30000
    unit_data = phy_data(phy_path, unit_srate)
    
    spikes = [None] * len(spike_templates)
    for i in range(len(spike_templates)):
        spike_train = unit_data.extract_spike_trains(spike_templates[i]) - seizure_start
        spike_train = spike_train[spike_train>-30]
        spike_train = spike_train[spike_train<50]
        spikes[i] = spike_train
        
    return spikes


files = os.listdir(path)

df = np.load(path + '\\' + files[0])
for i in range(len(files)):
    if i == 0:
        df = np.load(path + '\\' + files[0])
    else:
        df = np.row_stack((df, np.load(path + '\\' + files[i])))
        
df = df[:,:-1].astype('float64')
        
headers = ['template ID', 'effect_sz1', 'effect_sz2', 'effect_sz3', 'z_sz1', 'z_sz2','z_sz3', 'p_sz1', 'p_sz2', 'p_sz3']
df2 = pd.DataFrame(data=df, columns=headers)


### LL001
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
LL001_spike_templates = [5, 6, 10, 48, 59, 75, 91, 92, 101, 110, 114, 118, 125, 149, 164, 165, 176, 178, 180, 197]
seizure_start = 3782.72

LL001_all_spikes = spike_train_extract(
                             phy_path,
                             LL001_spike_templates,
                             seizure_start
                             )

### LL004

phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL004-final.GUI'
LL004_spike_templates = [18, 30, 40, 51, 52]
seizure_start = 1161.98

LL004_all_spikes = spike_train_extract(
                             phy_path,
                             LL004_spike_templates,
                             seizure_start
                             )


### LL005
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'
LL005_spike_templates = [17, 20, 26, 36, 38, 40, 47, 53, 61, 63, 77, 84, 95, 100, 112, 120, 124, 128, 136, 146]
seizure_start = 770.52

LL005_all_spikes = spike_train_extract(
                             phy_path,
                             LL005_spike_templates,
                             seizure_start
                             )

### LL007
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'
LL007_spike_templates = [81, 107, 112, 135, 157, 159, 181]
seizure_start = 1193.42

LL007_all_spikes = spike_train_extract(
                             phy_path,
                             LL007_spike_templates,
                             seizure_start
                             )

### LL009
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL009-final.GUI'
LL009_spike_templates = [28, 47]
seizure_start = 4408.35

LL009_all_spikes = spike_train_extract(
                             phy_path,
                             LL009_spike_templates,
                             seizure_start
                             )

### LL011ipsi
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL011_ipsi-final.GUI'
LL011ipsi_spike_templates = [16, 28]
seizure_start = 2650.96

LL011ipsi_all_spikes = spike_train_extract(
                             phy_path,
                             LL011ipsi_spike_templates,
                             seizure_start
                             )

### LL011contra
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL011_contra-final.GUI'
LL011contra_spike_templates = [15, 31, 33, 95, 101, 104, 107, 140, 150]
seizure_start = 542.118

LL011contra_all_spikes = spike_train_extract(
                             phy_path,
                             LL011contra_spike_templates,
                             seizure_start
                             )


### LL012ipsi

phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL012_ipsi-final.GUI'
LL012ipsi_spike_templates = [122, 137] ## 112 bugs out somehow
seizure_start = 729.362

LL012ipsi_all_spikes = spike_train_extract(
                             phy_path,
                             LL012ipsi_spike_templates,
                             seizure_start
                             )

### LL012contra
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL012_contra-final.GUI'
LL0012contra_spike_templates = [12, 65, 156, 166, 170, 200, 204, 218, 273, 310, 315, 317, 318, 329, 348, 354]
seizure_start = 124.806

LL012contra_all_spikes = spike_train_extract(
                             phy_path,
                             LL0012contra_spike_templates,
                             seizure_start
                             )


### LL017
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL017-final.GUI'
LL017_spike_templates = [276, 278, 293, 300, 317]
seizure_start = 3621.83

LL017_all_spikes = spike_train_extract(
                             phy_path,
                             LL017_spike_templates,
                             seizure_start
                             )

### LL023
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL023-final.GUI'
LL023_spike_templates = [4, 5, 9, 12, 17, 33, 48, 49, 52]
seizure_start = 2061.56

LL023_all_spikes = spike_train_extract(
                             phy_path,
                             LL023_spike_templates,
                             seizure_start
                             )


spikes = LL001_all_spikes + LL004_all_spikes + LL005_all_spikes + LL007_all_spikes + LL009_all_spikes + LL011contra_all_spikes+ LL011ipsi_all_spikes + LL012contra_all_spikes+ LL012ipsi_all_spikes + LL017_all_spikes + LL023_all_spikes

change = np.zeros((3,3))
effect_1 = []
effect_2 = []
effect_3 = []
x1=np.ones((np.shape(df)[0]))
x2=np.ones((np.shape(df)[0])) * 20
x3=np.ones((np.shape(df)[0])) * 40
for i in range(np.shape(df)[0]):
    effect_1.append('empty')
    effect_2.append('empty')
    effect_3.append('empty')

## SZ1

p_threshold = 0.05
counter_inhibited = 0
counter_excited = 0
counter_nochange = 0
for y in range(np.shape(df)[0]):
    if df[y,4] < 0 and df[y,7] < p_threshold:
        counter_inhibited = counter_inhibited + 1
        effect_1[y] = 'inhibited (p<0.05)'
    elif df[y,4] > 0 and df[y,7] < p_threshold:
        counter_excited = counter_excited + 1
        effect_1[y] = 'excited (p<0.05)'
    elif df[y,7] > p_threshold:
        counter_nochange = counter_nochange + 1
        effect_1[y] = 'no change (p>0.05)'

    change[0,0] = counter_inhibited
    change[1,0] = counter_excited
    change[2,0] = counter_nochange


##SZ2
counter_inhibited = 0
counter_excited = 0
counter_nochange = 0
for y in range(np.shape(df)[0]):
    if df[y,5] < 0 and df[y,8] < p_threshold:
        counter_inhibited = counter_inhibited + 1
        effect_2[y] = 'inhibited (p<0.05)'
    elif df[y,5] > 0 and df[y,8] < p_threshold:
        counter_excited = counter_excited + 1
        effect_2[y] = 'excited (p<0.05)'
    elif df[y,8] > p_threshold:
        counter_nochange = counter_nochange + 1
        effect_2[y] = 'no change (p>0.05)'

    change[0,1] = counter_inhibited
    change[1,1] = counter_excited
    change[2,1] = counter_nochange 
    
##SZ3
counter_inhibited = 0
counter_excited = 0
counter_nochange = 0
for y in range(np.shape(df)[0]):
    if df[y,6] < 0 and df[y,9] < p_threshold:
        counter_inhibited = counter_inhibited + 1
        effect_3[y] = 'inhibited (p<0.05)'
    elif df[y,6] > 0 and df[y,9] < p_threshold:
        counter_excited = counter_excited + 1
        effect_3[y] = 'excited (p<0.05)'
    elif df[y,9] > p_threshold:
        counter_nochange = counter_nochange + 1
        effect_3[y] = 'no change (p>0.05)'

    change[0,2] = counter_inhibited
    change[1,2] = counter_excited
    change[2,2] = counter_nochange



effect = np.array(df2['effect_sz1'].iloc[:97])
effect_label = effect_1
indices = list(np.argsort(effect).astype(int))

reshuffled_spikes = [spikes[i] for i in indices]
reshuffled_labels = [effect_label[i] for i in indices]


colors = [None]*len(reshuffled_labels)
for i in range(len(colors)):
    if reshuffled_labels[i] == 'inhibited (p<0.05)':
        colors[i] = 'b'
    elif reshuffled_labels[i] == 'excited (p<0.05)':
        colors[i] = 'r'
    elif reshuffled_labels[i] == 'no change (p>0.05)':
        colors[i] = 'g'


### PLOTTING
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(7.5,6))
ax.set_title('Seizure 1', fontsize=17, fontweight='bold')
ax.eventplot(reshuffled_spikes, color=colors, linelengths = 0.7)
ax.vlines(0, -0.5, 97, colors='r', linewidth=2, linestyles='dashed')
ax.vlines(10, -0.5, 97, colors='r', linewidth=2, linestyles='dashed')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(False)
ax.set_xlim(-10,40)
ax.set_ylim(-0.6,97)
ax.set_yticks([])
ax.set_xlabel('Time (s)', fontsize=17)
ax.set_ylabel('Neurons', fontsize=17)
plt.savefig(os.path.abspath(savepath), dpi=600)
plt.tight_layout()
#plt.close()