# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 12:15:00 2022

@author: llarsen
"""

import os
import numpy as np
import matplotlib.pyplot as plt


loaddir = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\input\pinch'
savedir = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\output'

files = os.listdir(loaddir)

data = np.load(loaddir + '\\' + files[0], allow_pickle=True).item()

for i in range(len(files)):
    if i == 0:
        pass
    else:
        new_data = np.load(loaddir + '\\' + files[i], allow_pickle=True).item()
        if len(new_data['LC_excited']) > 0:
            data['LC_excited'] = np.column_stack((data['LC_excited'], new_data['LC_excited']))
        if len(new_data['LC_inhibited']) > 0:
            data['LC_inhibited'] = np.column_stack((data['LC_inhibited'], new_data['LC_inhibited']))
        if len(new_data['LC_nochange']) > 0:
            data['LC_nochange'] = np.column_stack((data['LC_nochange'], new_data['LC_nochange']))
        if len(new_data['no_LC']) > 0:
            data['no_LC'] = np.column_stack((data['no_LC'], new_data['no_LC']))
            
t = np.arange(-20,6)
    
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(6,4))
ax.plot(t, np.mean(data['LC_excited'], axis = 1), color='r', label = 'Excited LC neurons (n=' + str(np.shape(data['LC_excited'])[1]) + ')')
ax.fill_between(t, np.mean(data['LC_excited'], axis=1), np.mean(data['LC_excited'], axis = 1) + (np.std(data['LC_excited'], axis = 1)/np.sqrt(len(data['LC_excited'])))*1.96, color='r', alpha=0.3)
ax.fill_between(t, np.mean(data['LC_excited'], axis=1), np.mean(data['LC_excited'], axis = 1) - (np.std(data['LC_excited'], axis = 1)/np.sqrt(len(data['LC_excited'])))*1.96, color='r', alpha=0.3)

ax.plot(t, np.mean(data['LC_inhibited'], axis = 1), color='b', label = 'Inhibited LC neurons (n=' + str(np.shape(data['LC_inhibited'])[1]) + ')')
ax.fill_between(t, np.mean(data['LC_inhibited'], axis=1), np.mean(data['LC_inhibited'], axis = 1) + (np.std(data['LC_inhibited'], axis = 1)/np.sqrt(len(data['LC_inhibited'])))*1.96, color='b', alpha=0.3)
ax.fill_between(t, np.mean(data['LC_inhibited'], axis=1), np.mean(data['LC_inhibited'], axis = 1) - (np.std(data['LC_inhibited'], axis = 1)/np.sqrt(len(data['LC_inhibited'])))*1.96, color='b', alpha=0.3)

ax.plot(t, np.mean(data['LC_nochange'], axis = 1), color='g', label = 'No change LC neurons (n=' + str(np.shape(data['LC_nochange'])[1]) + ')')
ax.fill_between(t, np.mean(data['LC_nochange'], axis=1), np.mean(data['LC_nochange'], axis = 1) + (np.std(data['LC_nochange'], axis = 1)/np.sqrt(len(data['LC_nochange'])))*1.96, color='g', alpha=0.3)
ax.fill_between(t, np.mean(data['LC_nochange'], axis=1), np.mean(data['LC_nochange'], axis = 1) - (np.std(data['LC_nochange'], axis = 1)/np.sqrt(len(data['LC_nochange'])))*1.96, color='g', alpha=0.3)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.set_xlim(-20, 5)
ax.set_xlabel('Time (s)', fontsize=14)
ax.set_ylabel('Spike Rate (Hz)', fontsize=14)
plt.legend(frameon=False, loc='upper left')
plt.tight_layout()
save_filename = '\e_pinch.png'
savepath = savedir + save_filename
plt.savefig(os.path.abspath(savepath), dpi=600)





''' FIGURE NOT USED IN MANUSCRIPT -- NON LC NEURONS
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(6,4))
ax.plot(t, np.mean(data['no_LC'], axis = 1), color='k', label = 'Non LC neurons (n=' + str(np.shape(data['no_LC'])[1]) + ')')
ax.fill_between(t, np.mean(data['no_LC'], axis=1), np.mean(data['no_LC'], axis = 1) + (np.std(data['no_LC'], axis = 1)/np.sqrt(len(data['no_LC'])))*1.96, color='k', alpha=0.3)
ax.fill_between(t, np.mean(data['no_LC'], axis=1), np.mean(data['no_LC'], axis = 1) - (np.std(data['no_LC'], axis = 1)/np.sqrt(len(data['no_LC'])))*1.96, color='k', alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.set_xlim(-20, 5)
ax.set_xlabel('Time (s)', fontsize=14)
ax.set_ylabel('Spike Rate (Hz)', fontsize=14)
plt.legend(frameon=False, loc='upper left')
plt.tight_layout()
#save_filename = r'\no_LC_pinch.png'
#savepath = savedir + save_filename
#plt.savefig(os.path.abspath(savepath), dpi=300)
#plt.close()
'''