# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 12:53:20 2022

@author: User
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, r'C:\Users\User\OneDrive - UGent\python_functions')


file_broad = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\pair_type_coupling_broad.npy'
file_sharp = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\pair_type_coupling_sharp.npy'

file_broad = np.load(file_broad, allow_pickle='TRUE')
file_sharp = np.load(file_sharp, allow_pickle='TRUE')

savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\output\Fig4E.png'


### BROAD PLOT EXAMPLE
x_broad = np.arange(1.3, 5.3, 1)
x_sharp = np.arange(1.7, 5.7, 1)

fig, ax = plt.subplots(1, figsize=(4,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
      
ax.bar(x_broad, file_broad*100, width = 0.4, color='gray', label = 'Broad')
ax.bar(x_sharp, file_sharp*100, width = 0.4, color='darkgray', label = 'Sharp')


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

ax.set_xticks([1.5, 2.5, 3.5, 4.5])
ax.set_xticklabels(['excited', 'inhibited', 'mixed', 'other']) 

ax.set_xlabel('Types of pairs', fontsize=14, fontweight='bold')
ax.set_ylabel('% coupled pairs', fontsize=14, fontweight='bold')
ax.set_ylim((0, 40))
ax.legend(frameon=False, loc='upper right')
plt.tight_layout()       
plt.savefig(savepath, dpi=600)