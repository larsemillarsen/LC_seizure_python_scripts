# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 17:56:31 2022

@author: User
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, r'E:\Manuscript_analysis_files')

from vies.spike.autocorrelogram import spike_train_correlogram
from vies.parse.phy import phy_data

from itertools import combinations
#import random
#from scipy.optimize import curve_fit
import time

unit_srate = 30000

LL001 = {}
#LL004 = {}
#LL005 = {}
#LL007 = {}
#LL009 = {}
#LL011_ipsi = {}
#LL011_contra = {}
#LL012_ipsi = {}
#LL012_contra = {}
#LL017 = {}
#LL023 = {}

LL001['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
#LL004['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL004-final.GUI'
#LL005['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'
#LL007['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'
#LL009['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL009-final.GUI'
#LL011_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL011_ipsi-final.GUI'
#LL011_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL01_contra-final.GUI'
#LL012_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_ipsi-final.GUI'
#LL012_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_contra-final.GUI'
#LL017['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL017-final.GUI'
#LL023['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL023-final.GUI'

LL001['templates'] = [5, 6, 48, 59, 91, 114, 125, 149, 164]
#LL004['templates'] = [51, 52]
#LL005['templates'] = [95, 120, 136]
#LL007['templates'] = [81, 157]
#LL009['templates'] = [47]
#LL011_contra['templates'] = [107]
#LL012_ipsi['templates'] = [122, 137]
#LL012_contra['templates'] = [318, 329]
#LL017['templates'] = [276, 293, 300, 317]
#LL023['templates'] = [4, 9, 12, 33, 48, 49]

LL001['seizure_start'] = [3782.72, 4570.89, 5191.40, 5803.33]
#LL004['seizure_start'] = [1161.98, 1877.26, 2482.50, 3301.26]
#LL005['seizure_start'] = [770.52, 1664.53, 2309.58, 2986.12]
#LL007['seizure_start'] = [1193.42, 1947.62, 2539.26, 3392.5, 4197.42]
#LL009['seizure_start'] = [4408.35, 6246.86, 7525.43, 9029.11]
#LL011_contra['seizure_start'] = [542.118, 1793.91, 7718.44]
#LL012_ipsi['seizure_start'] = [729.362]
#LL012_contra['seizure_start'] = [124.806, 1392.35, 5306.96, 5306.96]
#LL017['seizure_start'] = [964.981, 2297.39, 3621.83, 5063.79, 6353.71]
#LL023['seizure_start'] = [112.027, 2061.56, 2753.56, 3454.29, 5612.5, 6857.6]

iterations = list(combinations(range(len(LL001['templates'])), 2))

LL001['data'] = phy_data(LL001['phy_path'], unit_srate)

bins = 12
window = 6
bin_size = window/bins

x = np.arange(-window/2, window/2, bin_size) + bin_size/2

counter = 0
n_shuffles = 200
jitter = 1 # +- in ms

data_LL001 = {}
data_LL001['peak_latency'] = np.zeros((len(iterations)))
data_LL001['significant_rank'] = np.zeros((len(iterations)))
data_LL001['significant_parametric'] = np.zeros((len(iterations)))

data_LL001['correlograms'] = []
data_LL001['jittered_correlograms'] = []

data_LL001['median'] = []
data_LL001['local_1'] = []
data_LL001['local_99'] = []
data_LL001['global_1'] = []
data_LL001['global_99'] = []

data_LL001['mean'] = []
data_LL001['std'] = []
data_LL001['local_std_1'] = []
data_LL001['local_std_99'] = []
data_LL001['global_std_1'] = []
data_LL001['global_std_99'] = []

for i in range(len(iterations)):

    start_time = time.time()    
    template_1 = iterations[i][0]
    template_2 = iterations[i][1]
    
    spikes_1 = LL001['data'].extract_spike_trains(LL001['templates'][template_1])
    spikes_2 = LL001['data'].extract_spike_trains(LL001['templates'][template_2])

    count = 0
    for times in LL001['seizure_start']:
        if count == 0:
            seizure_index_1 = np.where(np.logical_and(spikes_1>times+10, spikes_1<times+60))[0]
            seizure_index_2 = np.where(np.logical_and(spikes_2>times+10, spikes_2<times+60))[0]
        else:
            seizure_index_1 = np.append(seizure_index_1, np.where(np.logical_and(spikes_1>times+10, spikes_1<times+60))[0])
            seizure_index_2 = np.append(seizure_index_2, np.where(np.logical_and(spikes_2>times+10, spikes_2<times+60))[0])
        count = count +1
            
    spikes_1 = spikes_1[seizure_index_1]   
    spikes_2 = spikes_2[seizure_index_2]
    
    spikes_1 = spikes_1 * 1000
    spikes_2 = spikes_2 * 1000    
    
    correllogram = spike_train_correlogram(spikes_1, spikes_2, bins=bins, window=window)
        
    jittered_correlogram = np.zeros((bins, n_shuffles))
    
    for k in range(n_shuffles):
        
        spikes_1 = LL001['data'].extract_spike_trains(LL001['templates'][template_1])
        spikes_2 = LL001['data'].extract_spike_trains(LL001['templates'][template_2])
        
        count = 0
        for times in LL001['seizure_start']:
            if count == 0:
                seizure_index_1 = np.where(np.logical_and(spikes_1>times+10, spikes_1<times+60))[0]
                seizure_index_2 = np.where(np.logical_and(spikes_2>times+10, spikes_2<times+60))[0]
            else:
                seizure_index_1 = np.append(seizure_index_1, np.where(np.logical_and(spikes_1>times+10, spikes_1<times+60))[0])
                seizure_index_2 = np.append(seizure_index_2, np.where(np.logical_and(spikes_2>times+10, spikes_2<times+60))[0])
            count = count +1
                
        spikes_1 = spikes_1[seizure_index_1]   
        spikes_2 = spikes_2[seizure_index_2]
            
        spikes_1 = spikes_1 * 1000
        spikes_2 = spikes_2 * 1000
        
        jit = np.around(np.random.uniform(-jitter, jitter, size=len(spikes_2)), 2)
        spikes_2 = spikes_2 + jit
        spike_train_correlogram(spikes_1, spikes_2, bins=bins, window=window)
        
        jittered_correlogram[:,k] = spike_train_correlogram(spikes_1, spikes_2, bins=bins, window=window)
        
    min_jitter = np.min(jittered_correlogram, axis=0)
    max_jitter = np.max(jittered_correlogram, axis=0)
    

    mean_jitter = np.mean(jittered_correlogram, axis=1)
    std_jitter = np.std(jittered_correlogram, axis=1)
    local_std_1 = np.mean(jittered_correlogram, axis=1) - 3*np.std(jittered_correlogram, axis=1)
    local_std_max_99 = np.mean(jittered_correlogram, axis=1) + 3*np.std(jittered_correlogram, axis=1)
    global_std_1 = np.mean(min_jitter) - 3*np.std(min_jitter)
    global_std_99 = np.mean(max_jitter) + 3*np.std(max_jitter)
    
    median_jitter = np.median(jittered_correlogram, axis=1)    
    local_99 = np.percentile(jittered_correlogram, 99, axis=1)
    local_1 = np.percentile(jittered_correlogram, 1, axis=1)
    global_1 = np.percentile(min_jitter, 1)
    global_99 = np.percentile(max_jitter, 99)

    data_LL001['correlograms'].append(correllogram)
    data_LL001['jittered_correlograms'].append(jittered_correlogram)
    
    data_LL001['median'].append(median_jitter)
    data_LL001['local_1'].append(local_1)
    data_LL001['local_99'].append(local_99)
    data_LL001['global_1'].append(global_1)
    data_LL001['global_99'].append(global_99)

    data_LL001['mean'].append(mean_jitter)
    data_LL001['std'].append(std_jitter)
    data_LL001['local_std_1'].append(local_std_1)
    data_LL001['local_std_99'].append(local_std_max_99)
    data_LL001['global_std_1'].append(global_std_1)
    data_LL001['global_std_99'].append(global_std_99)
    
    index_max = np.where((correllogram == np.max(correllogram)))[0][0]
    data_LL001['peak_latency'][i] = x[index_max]
    
    for latency in range(bins):
        if correllogram[latency] > local_99[latency] and correllogram[latency] > global_99:
            data_LL001['significant_rank'][i] = 1
    
    for latency in range(bins):
        if correllogram[latency] > local_std_max_99[latency] and correllogram[latency] > global_std_99:
            data_LL001['significant_parametric'][i] = 1
    
    ones = np.ones((len(correllogram)))
    
    plt.bar(x, correllogram, width = 0.9*(bin_size))

    plt.bar(x, correllogram, width=0.9*bin_size)
    plt.plot(x, mean_jitter, color='gray')
    plt.plot(x, local_std_max_99, linestyle='--', color='gray')
    plt.plot(x, local_std_1, linestyle='--', color='gray')
    
    plt.plot(x, global_std_99*ones, linestyle='dotted', color='blue')
    plt.plot(x, global_std_1*ones, linestyle='dotted', color='blue')
        
    plt.ylim((0, 1))
    
    savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz'

    if os.path.isdir(savepath):
        #print("Saving Directory Exists")
        pass
    else:
        print("Created Saving Directory")
        os.makedirs(savepath)
      
    #suppress sanity plots          
    #plt.savefig(savepath + '\\' + str(LL001['templates'][template_1]) + '_' + str(LL001['templates'][template_2]) + '_correlogram_ranked.png' , dpi=300)
    
    plt.close()
    end_time = time.time()
    print(str(i+1) + ' of ' + str(len(iterations)) + ', this iteration took ' + str(int(np.round(end_time - start_time))) + ' seconds')

np.save(r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL001_ranked.npy', data_LL001)
