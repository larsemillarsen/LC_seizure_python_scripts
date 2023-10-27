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
from itertools import combinations

savepath_broad = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\output\4A_broad_example.png'
savepath_sharp = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\output\4B_sharp_example.png'


LL005 = {}
LL023 = {}

LL005['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL005_ranked.npy'
LL023['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL023_ranked.npy'

LL005 = np.load(LL005['path'], allow_pickle='TRUE').item()
LL023 = np.load(LL023['path'], allow_pickle='TRUE').item()

LL005['templates'] = [17, 20, 26, 36, 38, 40, 47, 53, 61, 63, 77, 84, 95, 100, 112, 120, 124, 128, 136, 146]
LL023['templates'] = [4, 5, 9, 12, 33, 48, 49]

LL005['pairs'] = iterations = list(combinations(LL005['templates'], 2))
LL023['pairs'] = iterations = list(combinations(LL023['templates'], 2))

broad_example = LL005['pairs'][170]
sharp_example = LL023['pairs'][14]

window = 400
bins = 40
bin_size = window/bins
x = np.arange(-window/2, window/2, bin_size) + bin_size/2
ones = np.ones((len(x)))


### BROAD PLOT EXAMPLE
fig, ax = plt.subplots(1, figsize=(4,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
      
ax.bar(x, LL005['correlograms'][170], width = 0.80*(bin_size))

ax.plot(x, LL005['mean'][170], color='gray')
ax.plot(x, LL005['local_std_99'][170], linestyle='--', color='gray')
ax.plot(x, LL005['local_std_1'][170], linestyle='--', color='gray')

ax.plot(x, LL005['global_std_1'][170]*ones, linestyle='dotted', color='blue')
ax.plot(x, LL005['global_std_99'][170]*ones, linestyle='dotted', color='blue')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

ax.set_yticks([0, 0.02, 0.04, 0.06]) 
ax.set_yticklabels([0, 0.02, 0.04, 0.06], fontsize=12, fontweight='bold') 
ax.set_xticks([-200,0,200])
ax.set_xticklabels([-200,0,200], fontsize=12, fontweight='bold') 

ax.set_xlabel('Time (ms)', fontsize=14, fontweight='bold')
ax.set_ylabel('Normalized spike probability', fontsize=14, fontweight='bold')
ax.set_ylim((0, 0.06))
plt.tight_layout()       
plt.savefig(savepath_broad, dpi=600)

plt.close()

### BROAD PLOT EXAMPLE

window = 6
bins = 12
bin_size = window/bins
x = np.arange(-window/2, window/2, bin_size) + bin_size/2
ones = np.ones((len(x)))

fig, ax = plt.subplots(1, figsize=(4,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
      
ax.bar(x, LL023['correlograms'][14], width = 0.80*(bin_size))

ax.plot(x, LL023['mean'][14], color='gray')
ax.plot(x, LL023['local_std_99'][14], linestyle='--', color='gray')
ax.plot(x, LL023['local_std_1'][14], linestyle='--', color='gray')

ax.plot(x, LL023['global_std_1'][14]*ones, linestyle='dotted', color='blue')
ax.plot(x, LL023['global_std_99'][14]*ones, linestyle='dotted', color='blue')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

ax.set_yticks([0, 0.2, 0.4, 0.6]) 
ax.set_yticklabels([0, 0.2, 0.4, 0.6], fontsize=12, fontweight='bold') 
ax.set_xticks([-3,0,3])
ax.set_xticklabels([-3,0,3], fontsize=12, fontweight='bold') 

ax.set_xlabel('Time (ms)', fontsize=14, fontweight='bold')
ax.set_ylabel('Normalized spike probability', fontsize=14, fontweight='bold')
ax.set_ylim((0, 0.6))
plt.tight_layout()       
plt.savefig(savepath_sharp, dpi=600)



        