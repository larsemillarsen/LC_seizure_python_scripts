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

from scipy.stats import pearsonr, spearmanr

sys.path.insert(1, r'E:\OneDrive - UGent\python_functions')

#path = r'E:\OneDrive - UGent\88_LCseizureproject\results\results_seizure_final_mod.xlsx'
path = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\Fig2H\output_data'

savepath_peak_asymmetry = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\output\unused_asymmetry.png'
savepath_mean_spike_rate = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\output\unused_mean_spike_rate.png'
savepath_median_ISI = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\output\unused_median_ISI.png'
savepath_box_spike_width = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\output\a_spike_width.png'
savepath_box_spike_asymmetry = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\output\b_spike_asymmetry.png'
savepath_box_spike_rate_mean = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\output\c_spike_rate_mean.png'
savepath_box_spike_ISI_median = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\output\d_spike_ISI_median.png'

files = os.listdir(path)

order_files = [0, 1, 2, 3, 4, 6, 5, 8, 7, 9, 10]
df = np.load(path + '\\' + files[0])
for i in order_files:
    if i == 0:
        df = np.load(path + '\\' + files[0])
    else:
        df = np.row_stack((df, np.load(path + '\\' + files[i])))
        
df = df[:,:-1].astype('float64')
        


change = np.zeros((3,3))
effect_1 = [None]*97
p_threshold = 0.05
counter_inhibited = 0
counter_excited = 0
counter_nochange = 0
for y in range(np.shape(df)[0]):
    if df[y,4] < 0 and df[y,7] < p_threshold:
        counter_inhibited = counter_inhibited + 1
        effect_1[y] = 'inhibited'
    elif df[y,4] > 0 and df[y,7] < p_threshold:
        counter_excited = counter_excited + 1
        effect_1[y] = 'excited'
    elif df[y,7] > p_threshold:
        counter_nochange = counter_nochange + 1
        effect_1[y] = 'no change'

    change[0,0] = counter_inhibited
    change[1,0] = counter_excited
    change[2,0] = counter_nochange

headers = ['template ID', 'effect_sz1', 'effect_sz2', 'effect_sz3', 'z_sz1', 'z_sz2','z_sz3', 'p_sz1', 'p_sz2', 'p_sz3']
seizure_responses = pd.DataFrame(data=df, columns=headers)
seizure_responses['effect'] = effect_1


firing_characteristics = np.load(r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\intermediate_data\firing_characteristics.npy', allow_pickle='TRUE').item()
firing_characteristics['median_ISI'] = np.load(r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\intermediate_data\median_ISI.npy')



waveform_data = np.load(r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\intermediate_data\width_asymmetry.npy', allow_pickle='TRUE').item()

seizure_responses['mean_firing'] = firing_characteristics['mean_firing']
seizure_responses['median_ISI'] = firing_characteristics['median_ISI']
seizure_responses['spike_width'] = waveform_data['widths']
seizure_responses['peak_asymmetry'] = waveform_data['peak_asymmetry']

d_exicted = seizure_responses[seizure_responses['effect']=='excited']
d_inhibited = seizure_responses[seizure_responses['effect']=='inhibited']
d_nochange = seizure_responses[seizure_responses['effect']=='no change']

### Peak Asymmetry

fig, ax = plt.subplots(1, figsize=(4,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
      
ax.scatter(seizure_responses['spike_width'][seizure_responses['effect']=='excited'],  seizure_responses['peak_asymmetry'][seizure_responses['effect']=='excited'], c='red', s=15)
ax.scatter(seizure_responses['spike_width'][seizure_responses['effect']=='inhibited'],  seizure_responses['peak_asymmetry'][seizure_responses['effect']=='inhibited'], c='blue', s=15)
ax.scatter(seizure_responses['spike_width'][seizure_responses['effect']=='no change'],  seizure_responses['peak_asymmetry'][seizure_responses['effect']=='no change'], c='green', s=15)
#ax.scatter(x_sharp, fraction_sh[0:4,0], width = 0.4)


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

ax.set_yticks([-1.0, -0.5, 0, 0.5, 1.0]) 
ax.set_yticklabels([-1.0, -0.5, 0, 0.5, 1.0], fontsize=12, fontweight='bold') 
ax.set_xticks([0, 0.5, 1, 1.5, 2])
ax.set_xticklabels([0, 0.5, 1, 1.5, 2], fontsize=12, fontweight='bold') 

ax.set_xlabel('Spike Width (ms)', fontsize=14, fontweight='bold')
ax.set_ylabel('Peak asymmetry (B-A)/(B+A)', fontsize=14, fontweight='bold')
ax.set_ylim((-1, 1))
ax.set_xlim((0, 2))
plt.tight_layout()       
plt.savefig(savepath_peak_asymmetry, dpi=300)
plt.close()

### mean_spikerate
fig, ax = plt.subplots(1, figsize=(4,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
      
#ax.scatter(waveform_data['widths'],  firing_characteristics['mean_firing'])
ax.scatter(seizure_responses['spike_width'][seizure_responses['effect']=='excited'],  seizure_responses['mean_firing'][seizure_responses['effect']=='excited'], c='red', s=15)
ax.scatter(seizure_responses['spike_width'][seizure_responses['effect']=='inhibited'],  seizure_responses['mean_firing'][seizure_responses['effect']=='inhibited'], c='blue', s=15)
ax.scatter(seizure_responses['spike_width'][seizure_responses['effect']=='no change'],  seizure_responses['mean_firing'][seizure_responses['effect']=='no change'], c='green', s=15)

#ax.scatter(x_sharp, fraction_sh[0:4,0], width = 0.4)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

ax.set_yticks([0, 2, 4, 6, 8]) 
ax.set_yticklabels([0, 2, 4, 6, 8], fontsize=12, fontweight='bold') 
ax.set_xticks([0, 0.5, 1, 1.5, 2])
ax.set_xticklabels([0, 0.5, 1, 1.5, 2], fontsize=12, fontweight='bold') 

ax.set_xlabel('Spike Width (ms)', fontsize=14, fontweight='bold')
ax.set_ylabel('Mean spike rate (Hz)', fontsize=14, fontweight='bold')
ax.set_ylim((0, 8))
ax.set_xlim((0, 2))
plt.tight_layout()       
plt.savefig(savepath_mean_spike_rate, dpi=300)
plt.close()

### median ISI
fig, ax = plt.subplots(1, figsize=(4,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
      
#ax.scatter(waveform_data['widths'],  firing_characteristics['median_ISI'])
ax.scatter(seizure_responses['spike_width'][seizure_responses['effect']=='excited'],  seizure_responses['median_ISI'][seizure_responses['effect']=='excited'], c='red', s=15)
ax.scatter(seizure_responses['spike_width'][seizure_responses['effect']=='inhibited'],  seizure_responses['median_ISI'][seizure_responses['effect']=='inhibited'], c='blue', s=15)
ax.scatter(seizure_responses['spike_width'][seizure_responses['effect']=='no change'],  seizure_responses['median_ISI'][seizure_responses['effect']=='no change'], c='green', s=15)

#ax.scatter(x_sharp, fraction_sh[0:4,0], width = 0.4)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0, 2.5]) 
ax.set_yticklabels([0, 0.5, 1.0, 1.5, 2.0, 2.5], fontsize=12, fontweight='bold') 
ax.set_xticks([0, 0.5, 1, 1.5, 2])
ax.set_xticklabels([0, 0.5, 1, 1.5, 2], fontsize=12, fontweight='bold') 

ax.set_xlabel('Spike Width (ms)', fontsize=14, fontweight='bold')
ax.set_ylabel('Median spike interval (s)', fontsize=14, fontweight='bold')
ax.set_ylim((0, 2.5))
ax.set_xlim((0, 2))
plt.tight_layout()       
plt.savefig(savepath_median_ISI, dpi=300)
plt.close()



### spike_width boxplot
meanprops = {'linestyle': 'solid',
            'color':'black',
            'linewidth':2}

medianprops  = {'linestyle': '--',
            'color':'black',
            'linewidth':2}


fig, ax = plt.subplots(1, figsize=(3,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

x_excited = np.linspace(1.7, 1.75, len(seizure_responses['spike_width'][seizure_responses['effect']=='excited']))
x_inhibited = np.linspace(0.7, 0.75, len(seizure_responses['spike_width'][seizure_responses['effect']=='inhibited']))
x_no_change = np.linspace(2.7, 2.75, len(seizure_responses['spike_width'][seizure_responses['effect']=='no change']))      

box_plot_data = [seizure_responses['spike_width'][seizure_responses['effect']=='inhibited'],
                 seizure_responses['spike_width'][seizure_responses['effect']=='excited'],
                 seizure_responses['spike_width'][seizure_responses['effect']=='no change']
                ]

#ax.scatter(waveform_data['widths'],  firing_characteristics['median_ISI'])
ax.scatter(x_excited, seizure_responses['spike_width'][seizure_responses['effect']=='excited'], c='red', s=15)
ax.scatter(x_inhibited, seizure_responses['spike_width'][seizure_responses['effect']=='inhibited'], c='blue', s=15)
ax.scatter(x_no_change, seizure_responses['spike_width'][seizure_responses['effect']=='no change'], c='green', s=15)

#ax.boxplot(seizure_responses['spike_width'][seizure_responses['effect']=='excited'], positions=[1], showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)
ax.boxplot(box_plot_data, positions=[1, 2, 3], showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)

#ax.scatter(x_sharp, fraction_sh[0:4,0], width = 0.4)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

ax.set_yticks([0.5, 1.0, 1.5]) 
ax.set_yticklabels([0.5, 1.0, 1.5], fontsize=12, fontweight='bold') 
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['inhibted', 'excited', 'no change'], fontsize=14, fontweight='bold', rotation=45) 

#ax.set_xlabel('Spike Width (ms)', fontsize=14, fontweight='bold')
ax.set_ylabel('Spike Width (ms)', fontsize=14, fontweight='bold')
#ax.set_ylim((0, 2))
#ax.set_xlim((0, 2))
plt.tight_layout()       
plt.savefig(savepath_box_spike_width, dpi=300)
plt.close()



### spike_asymmetry boxplot
meanprops = {'linestyle': 'solid',
            'color':'black',
            'linewidth':2}

medianprops  = {'linestyle': '--',
            'color':'black',
            'linewidth':2}


fig, ax = plt.subplots(1, figsize=(3,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

x_excited = np.linspace(1.7, 1.75, len(seizure_responses['peak_asymmetry'][seizure_responses['effect']=='excited']))
x_inhibited = np.linspace(0.7, 0.75, len(seizure_responses['peak_asymmetry'][seizure_responses['effect']=='inhibited']))
x_no_change = np.linspace(2.7, 2.75, len(seizure_responses['peak_asymmetry'][seizure_responses['effect']=='no change']))      

box_plot_data = [seizure_responses['peak_asymmetry'][seizure_responses['effect']=='inhibited'],
                 seizure_responses['peak_asymmetry'][seizure_responses['effect']=='excited'],
                 seizure_responses['peak_asymmetry'][seizure_responses['effect']=='no change']
                ]

#ax.scatter(waveform_data['widths'],  firing_characteristics['median_ISI'])
ax.scatter(x_excited, seizure_responses['peak_asymmetry'][seizure_responses['effect']=='excited'], c='red', s=15)
ax.scatter(x_inhibited, seizure_responses['peak_asymmetry'][seizure_responses['effect']=='inhibited'], c='blue', s=15)
ax.scatter(x_no_change, seizure_responses['peak_asymmetry'][seizure_responses['effect']=='no change'], c='green', s=15)

#ax.boxplot(seizure_responses['peak_asymmetry'][seizure_responses['effect']=='excited'], positions=[1], showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)
ax.boxplot(box_plot_data, positions=[1, 2, 3], showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)

#ax.scatter(x_sharp, fraction_sh[0:4,0], width = 0.4)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

ax.set_yticks([-1.0, -0.5, 0]) 
ax.set_yticklabels([-1.0, -0.5, 0], fontsize=12, fontweight='bold') 
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['inhibted', 'excited', 'no change'], fontsize=14, fontweight='bold', rotation=45) 

#ax.set_xlabel('Spike Width (ms)', fontsize=14, fontweight='bold')
ax.set_ylabel('Peak asymmetry', fontsize=14, fontweight='bold')
ax.set_ylim((-1, 0.2))
#ax.set_xlim((0, 2))
plt.tight_layout()       
plt.savefig(savepath_box_spike_asymmetry, dpi=300)
plt.close()



### spike_rate_mean boxplot
meanprops = {'linestyle': 'solid',
            'color':'black',
            'linewidth':2}

medianprops  = {'linestyle': '--',
            'color':'black',
            'linewidth':2}


fig, ax = plt.subplots(1, figsize=(3,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

x_excited = np.linspace(1.7, 1.75, len(seizure_responses['mean_firing'][seizure_responses['effect']=='excited']))
x_inhibited = np.linspace(0.7, 0.75, len(seizure_responses['mean_firing'][seizure_responses['effect']=='inhibited']))
x_no_change = np.linspace(2.7, 2.75, len(seizure_responses['mean_firing'][seizure_responses['effect']=='no change']))      

box_plot_data = [seizure_responses['mean_firing'][seizure_responses['effect']=='inhibited'],
                 seizure_responses['mean_firing'][seizure_responses['effect']=='excited'],
                 seizure_responses['mean_firing'][seizure_responses['effect']=='no change']
                ]

#ax.scatter(waveform_data['widths'],  firing_characteristics['median_ISI'])
ax.scatter(x_excited, seizure_responses['mean_firing'][seizure_responses['effect']=='excited'], c='red', s=15)
ax.scatter(x_inhibited, seizure_responses['mean_firing'][seizure_responses['effect']=='inhibited'], c='blue', s=15)
ax.scatter(x_no_change, seizure_responses['mean_firing'][seizure_responses['effect']=='no change'], c='green', s=15)

#ax.boxplot(seizure_responses['mean_firing'][seizure_responses['effect']=='excited'], positions=[1], showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)
ax.boxplot(box_plot_data, positions=[1, 2, 3], showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)

#ax.scatter(x_sharp, fraction_sh[0:4,0], width = 0.4)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

ax.set_yticks([0, 2, 4, 6, 8]) 
ax.set_yticklabels([0, 2, 4, 6, 8], fontsize=12, fontweight='bold') 
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['inhibted', 'excited', 'no change'], fontsize=14, fontweight='bold', rotation=45) 

#ax.set_xlabel('Spike Width (ms)', fontsize=14, fontweight='bold')
ax.set_ylabel('Mean spike rate (Hz)', fontsize=14, fontweight='bold')
#ax.set_ylim((-1, 0.2))
#ax.set_xlim((0, 2))
plt.tight_layout()       
plt.savefig(savepath_box_spike_rate_mean, dpi=300)
plt.close()




### spike_ISI_median boxplot
meanprops = {'linestyle': 'solid',
            'color':'black',
            'linewidth':2}

medianprops  = {'linestyle': '--',
            'color':'black',
            'linewidth':2}


fig, ax = plt.subplots(1, figsize=(3,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

x_excited = np.linspace(1.7, 1.75, len(seizure_responses['median_ISI'][seizure_responses['effect']=='excited']))
x_inhibited = np.linspace(0.7, 0.75, len(seizure_responses['median_ISI'][seizure_responses['effect']=='inhibited']))
x_no_change = np.linspace(2.7, 2.75, len(seizure_responses['median_ISI'][seizure_responses['effect']=='no change']))      

box_plot_data = [seizure_responses['median_ISI'][seizure_responses['effect']=='inhibited'],
                 seizure_responses['median_ISI'][seizure_responses['effect']=='excited'],
                 seizure_responses['median_ISI'][seizure_responses['effect']=='no change']
                ]

#ax.scatter(waveform_data['widths'],  firing_characteristics['median_ISI'])
ax.scatter(x_excited, seizure_responses['median_ISI'][seizure_responses['effect']=='excited'], c='red', s=15)
ax.scatter(x_inhibited, seizure_responses['median_ISI'][seizure_responses['effect']=='inhibited'], c='blue', s=15)
ax.scatter(x_no_change, seizure_responses['median_ISI'][seizure_responses['effect']=='no change'], c='green', s=15)

#ax.boxplot(seizure_responses['median_ISI'][seizure_responses['effect']=='excited'], positions=[1], showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)
ax.boxplot(box_plot_data, positions=[1, 2, 3], showmeans=True, meanline=True, meanprops=meanprops, medianprops=medianprops)

#ax.scatter(x_sharp, fraction_sh[0:4,0], width = 0.4)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0]) 
ax.set_yticklabels([0, 0.5, 1.0, 1.5, 2.0], fontsize=12, fontweight='bold') 
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['inhibted', 'excited', 'no change'], fontsize=14, fontweight='bold', rotation=45) 

#ax.set_xlabel('Spike Width (ms)', fontsize=14, fontweight='bold')
ax.set_ylabel('Median spike interval (s)', fontsize=14, fontweight='bold')
#ax.set_ylim((-1, 0.2))
#ax.set_xlim((0, 2))
plt.tight_layout()       
plt.savefig(savepath_box_spike_ISI_median, dpi=300)
plt.close()





