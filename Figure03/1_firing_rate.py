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

sys.path.insert(1, r'E:\Manuscript_analysis_files')
from vies.parse.phy import phy_data

unit_srate = 30000

savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\intermediate_data\firing_characteristics.npy'

LL001 = {}
LL004 = {}
LL005 = {}
LL007 = {}
LL009 = {}
LL011_ipsi = {}
LL011_contra = {}
LL012_ipsi = {}
LL012_contra = {}
LL017 = {}
LL023 = {}

LL001['templates'] = [5, 6, 10, 48, 59, 75, 91, 92, 101, 110, 114, 118, 125, 149, 164, 165, 176, 178, 180, 197]
LL004['templates'] = [18, 30, 40, 51, 52]
LL005['templates'] = [17, 20, 26, 36, 38, 40, 47, 53, 61, 63, 77, 84, 95, 100, 112, 120, 124, 128, 136, 146]
LL007['templates'] = [81, 107, 112, 135, 157, 159, 181]
LL009['templates'] = [28, 47]
LL011_ipsi['templates'] = [16, 28]
LL011_contra['templates'] = [15, 31, 33, 95, 101, 104, 107, 140, 150]
LL012_ipsi['templates'] = [122, 137]
LL012_contra['templates'] = [12, 65, 156, 166, 170, 200, 204, 218, 273, 310, 315, 317, 318, 329, 348, 354]
LL017['templates'] = [276, 278, 293, 300, 317]
LL023['templates'] = [4, 5, 9, 12, 17, 33, 48, 49, 52]

LL001['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
LL004['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL004-final.GUI'
LL005['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'
LL007['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'
LL009['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL009-final.GUI'
LL011_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL011_ipsi-final.GUI'
LL011_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL011_contra-final.GUI'
LL012_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_ipsi-final.GUI'
LL012_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_contra-final.GUI'
LL017['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL017-final.GUI'
LL023['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL023-final.GUI'

LL001['seizure_start'] = np.array([3782.72, 4570.89, 5191.40, 5803.33])
LL004['seizure_start'] = np.array([1161.98, 1877.26, 2482.50, 3301.26])
LL005['seizure_start'] = np.array([770.52, 1664.53, 2309.58, 2986.12])
LL007['seizure_start'] = np.array([1193.42, 1947.62, 2539.26, 3392.5, 4197.42])
LL009['seizure_start'] = np.array([4408.35, 5048.82, 6246.86, 6872.27, 7525.43, 8147.10, 9029.11, 9649.13])
LL011_ipsi['seizure_start'] = np.array([2650.96, 3337.68, 3887.28])
LL011_contra['seizure_start'] = np.array([542.118, 1182.99, 1793.91, 2643.96, 7106.03, 7718.44])
LL012_ipsi['seizure_start'] = np.array([729.362, 1421.73])
LL012_contra['seizure_start'] = np.array([124.806, 742.641, 1392.35, 2008.74, 5306.96, 5925.43])
LL017['seizure_start'] = np.array([964.981, 1626.73, 2297.39, 2966.71, 3621.83, 4365.61, 5063.79, 5712.88, 6353.71, 7003.33])
LL023['seizure_start'] = np.array([112.027, 2061.56, 2753.56, 3454.29, 4110.26, 4878.67, 5612.5, 6238, 6857.6, 7479.57])

LL001['data'] = phy_data(LL001['phy_path'], unit_srate)
LL004['data'] = phy_data(LL004['phy_path'], unit_srate)
LL005['data'] = phy_data(LL005['phy_path'], unit_srate)
LL007['data'] = phy_data(LL007['phy_path'], unit_srate)
LL009['data'] = phy_data(LL009['phy_path'], unit_srate)
LL011_ipsi['data'] = phy_data(LL011_ipsi['phy_path'], unit_srate)
LL011_contra['data'] = phy_data(LL011_contra['phy_path'], unit_srate)
LL012_ipsi['data'] = phy_data(LL012_ipsi['phy_path'], unit_srate)
LL012_contra['data'] = phy_data(LL012_contra['phy_path'], unit_srate)
LL017['data'] = phy_data(LL017['phy_path'], unit_srate)
LL023['data'] = phy_data(LL023['phy_path'], unit_srate)

LL001['spike_freq_bins'] = LL001['data'].extract_frequency_bins(60, 0, 10000, smooth_method = 'none', templates = LL001['templates'])
LL004['spike_freq_bins'] = LL004['data'].extract_frequency_bins(60, 0, 10000, templates = LL004['templates'])
LL005['spike_freq_bins'] = LL005['data'].extract_frequency_bins(60, 0, 10000, templates = LL005['templates'])
LL007['spike_freq_bins'] = LL007['data'].extract_frequency_bins(60, 0, 10000, templates = LL007['templates'])
LL009['spike_freq_bins'] = LL009['data'].extract_frequency_bins(60, 0, 10000, templates = LL009['templates'])
LL011_ipsi['spike_freq_bins'] = LL011_ipsi['data'].extract_frequency_bins(60, 0, 10000, templates = LL011_ipsi['templates'])
LL011_contra['spike_freq_bins'] = LL011_contra['data'].extract_frequency_bins(60, 0, 10000, templates = LL011_contra['templates'])
LL012_ipsi['spike_freq_bins'] = LL012_ipsi['data'].extract_frequency_bins(60, 0, 10000, templates = LL012_ipsi['templates'])
LL012_contra['spike_freq_bins'] = LL012_contra['data'].extract_frequency_bins(60, 0, 10000, templates = LL012_contra['templates'])
LL017['spike_freq_bins'] = LL017['data'].extract_frequency_bins(60, 0, 10000, templates = LL017['templates'])
LL023['spike_freq_bins'] = LL023['data'].extract_frequency_bins(60, 0, 10000, templates = LL023['templates'])

for i in range(np.shape(LL001['spike_freq_bins'][1])[0]):
    LL001['spike_freq_bins'][1][i,:][LL001['spike_freq_bins'][1][i,:] == 0] = np.nan
    #LL001['spike_freq_bins'][1][i,:][LL001['spike_freq_bins'][1][i,:] > 6] = np.nan
    LL004['spike_freq_bins'][1][i,:][LL004['spike_freq_bins'][1][i,:] == 0] = np.nan
    #LL004['spike_freq_bins'][1][i,:][LL004['spike_freq_bins'][1][i,:] > 6] = np.nan
    LL005['spike_freq_bins'][1][i,:][LL005['spike_freq_bins'][1][i,:] == 0] = np.nan
    #LL005['spike_freq_bins'][1][i,:][LL005['spike_freq_bins'][1][i,:] > 6] = np.nan
    LL007['spike_freq_bins'][1][i,:][LL007['spike_freq_bins'][1][i,:] == 0] = np.nan
    #LL007['spike_freq_bins'][1][i,:][LL007['spike_freq_bins'][1][i,:] > 6] = np.nan
    LL009['spike_freq_bins'][1][i,:][LL009['spike_freq_bins'][1][i,:] == 0] = np.nan    
    #LL009['spike_freq_bins'][1][i,:][LL009['spike_freq_bins'][1][i,:] > 6] = np.nan
    LL011_ipsi['spike_freq_bins'][1][i,:][LL011_ipsi['spike_freq_bins'][1][i,:] == 0] = np.nan
    #LL011_ipsi['spike_freq_bins'][1][i,:][LL011_ipsi['spike_freq_bins'][1][i,:] > 6] = np.nan
    LL011_contra['spike_freq_bins'][1][i,:][LL011_contra['spike_freq_bins'][1][i,:] == 0] = np.nan
    #LL011_contra['spike_freq_bins'][1][i,:][LL011_contra['spike_freq_bins'][1][i,:] > 6] = np.nan
    LL012_ipsi['spike_freq_bins'][1][i,:][LL012_ipsi['spike_freq_bins'][1][i,:] == 0] = np.nan
    #LL012_ipsi['spike_freq_bins'][1][i,:][LL012_ipsi['spike_freq_bins'][1][i,:] > 6] = np.nan
    LL012_contra['spike_freq_bins'][1][i,:][LL012_contra['spike_freq_bins'][1][i,:] == 0] = np.nan
    #LL012_contra['spike_freq_bins'][1][i,:][LL012_contra['spike_freq_bins'][1][i,:] > 6] = np.nan
    LL017['spike_freq_bins'][1][i,:][LL017['spike_freq_bins'][1][i,:] == 0] = np.nan
    #LL017['spike_freq_bins'][1][i,:][LL017['spike_freq_bins'][1][i,:] > 6] = np.nan
    LL023['spike_freq_bins'][1][i,:][LL023['spike_freq_bins'][1][i,:] == 0] = np.nan
    #LL023['spike_freq_bins'][1][i,:][LL023['spike_freq_bins'][1][i,:] > 6] = np.nan


LL001['seizure_start_60'] = LL001['seizure_start'] / 60
LL004['seizure_start_60'] = LL004['seizure_start'] / 60
LL005['seizure_start_60'] = LL005['seizure_start'] / 60
LL007['seizure_start_60'] = LL007['seizure_start'] / 60
LL009['seizure_start_60'] = LL009['seizure_start'] / 60
LL011_ipsi['seizure_start_60'] = LL011_ipsi['seizure_start'] / 60
LL011_contra['seizure_start_60'] = LL011_contra['seizure_start'] / 60
LL012_ipsi['seizure_start_60'] = LL012_ipsi['seizure_start'] / 60
LL012_contra['seizure_start_60'] = LL012_contra['seizure_start'] / 60
LL017['seizure_start_60'] = LL017['seizure_start'] / 60
LL023['seizure_start_60'] = LL023['seizure_start'] / 60

pre_window = 1
post_window = 2
for i in LL001['seizure_start_60']:
    start = int(np.floor(i) - pre_window)
    stop = start + post_window
    LL001['spike_freq_bins'][1][start:stop, :] = np.nan

for i in LL004['seizure_start_60']:
    start = int(np.floor(i) - pre_window)
    stop = start + post_window
    LL004['spike_freq_bins'][1][start:stop, :] = np.nan
    
for i in LL005['seizure_start_60']:
    start = int(np.floor(i) - pre_window)
    stop = start + post_window
    LL005['spike_freq_bins'][1][start:stop, :] = np.nan
    
for i in LL007['seizure_start_60']:
    start = int(np.floor(i) - pre_window)
    stop = start + post_window
    LL007['spike_freq_bins'][1][start:stop, :] = np.nan
    
for i in LL009['seizure_start_60']:
    start = int(np.floor(i) - pre_window)
    stop = start + post_window
    LL009['spike_freq_bins'][1][start:stop, :] = np.nan

for i in LL011_ipsi['seizure_start_60']:
    start = int(np.floor(i) - pre_window)
    stop = start + post_window
    LL011_ipsi['spike_freq_bins'][1][start:stop, :] = np.nan
    
for i in LL011_contra['seizure_start_60']:
    start = int(np.floor(i) - pre_window)
    stop = start + post_window
    LL011_contra['spike_freq_bins'][1][start:stop, :] = np.nan
    
for i in LL012_ipsi['seizure_start_60']:
    start = int(np.floor(i) - pre_window)
    stop = start + post_window
    LL012_ipsi['spike_freq_bins'][1][start:stop, :] = np.nan

for i in LL012_contra['seizure_start_60']:
    start = int(np.floor(i) - pre_window)
    stop = start + post_window
    LL012_contra['spike_freq_bins'][1][start:stop, :] = np.nan

for i in LL017['seizure_start_60']:
    start = int(np.floor(i) - pre_window)
    stop = start + post_window
    LL017['spike_freq_bins'][1][start:stop, :] = np.nan

for i in LL023['seizure_start_60']:
    start = int(np.floor(i) - pre_window)
    stop = start + post_window
    LL023['spike_freq_bins'][1][start:stop, :] = np.nan

all_neurons = np.column_stack((
    LL001['spike_freq_bins'][1],
    LL004['spike_freq_bins'][1],
    LL005['spike_freq_bins'][1],
    LL007['spike_freq_bins'][1],
    LL009['spike_freq_bins'][1],
    LL011_ipsi['spike_freq_bins'][1],
    LL011_contra['spike_freq_bins'][1],
    LL012_ipsi['spike_freq_bins'][1],
    LL012_contra['spike_freq_bins'][1],
    LL017['spike_freq_bins'][1],
    LL023['spike_freq_bins'][1]
    ))


#test = LL001['spike_freq_bins'][1][[LL001['spike_freq_bins'][1] == 0]]

#test = LL001['spike_freq_bins']

#test = test[1]
data = {}
data['mean_firing'] = np.nanmean(all_neurons, axis = 0)
data['grand_std'] = np.nanstd(data['mean_firing'])
data['grand_mean'] = np.nanmean(data['mean_firing'])
data['grand_median'] = np.nanmedian(data['mean_firing'])

np.save(savepath, data)





