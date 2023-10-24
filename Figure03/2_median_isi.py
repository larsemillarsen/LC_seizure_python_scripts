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

savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\median_ISI.npy'

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

LL001['spike_times'] = LL001['data'].extract_spike_trains(LL001['templates'])
LL004['spike_times'] = LL004['data'].extract_spike_trains(LL004['templates'])
LL005['spike_times'] = LL005['data'].extract_spike_trains(LL005['templates'])
LL007['spike_times'] = LL007['data'].extract_spike_trains(LL007['templates'])
LL009['spike_times'] = LL009['data'].extract_spike_trains(LL009['templates'])
LL011_ipsi['spike_times'] = LL011_ipsi['data'].extract_spike_trains(LL011_ipsi['templates'])
LL011_contra['spike_times'] = LL011_contra['data'].extract_spike_trains(LL011_contra['templates'])
LL012_ipsi['spike_times'] = LL012_ipsi['data'].extract_spike_trains(LL012_ipsi['templates'])
LL012_contra['spike_times'] = LL012_contra['data'].extract_spike_trains(LL012_contra['templates'])
LL017['spike_times'] = LL017['data'].extract_spike_trains(LL017['templates'])
LL023['spike_times'] = LL023['data'].extract_spike_trains(LL023['templates'])

for times in LL001['seizure_start']:
    for k in range(len(LL001['templates'])):
        LL001['spike_times'][k] = LL001['spike_times'][k][np.invert(np.logical_and(LL001['spike_times'][k]>times, LL001['spike_times'][k]<times+60))]   
        
for times in LL004['seizure_start']:
    for k in range(len(LL004['templates'])):
        LL004['spike_times'][k] = LL004['spike_times'][k][np.invert(np.logical_and(LL004['spike_times'][k]>times, LL004['spike_times'][k]<times+60))]   

for times in LL005['seizure_start']:
    for k in range(len(LL005['templates'])):
        LL005['spike_times'][k] = LL005['spike_times'][k][np.invert(np.logical_and(LL005['spike_times'][k]>times, LL005['spike_times'][k]<times+60))]   

for times in LL007['seizure_start']:
    for k in range(len(LL007['templates'])):
        LL007['spike_times'][k] = LL007['spike_times'][k][np.invert(np.logical_and(LL007['spike_times'][k]>times, LL007['spike_times'][k]<times+60))]   

for times in LL009['seizure_start']:
    for k in range(len(LL009['templates'])):
        LL009['spike_times'][k] = LL009['spike_times'][k][np.invert(np.logical_and(LL009['spike_times'][k]>times, LL009['spike_times'][k]<times+60))]   

for times in LL011_ipsi['seizure_start']:
    for k in range(len(LL011_ipsi['templates'])):
        LL011_ipsi['spike_times'][k] = LL011_ipsi['spike_times'][k][np.invert(np.logical_and(LL011_ipsi['spike_times'][k]>times, LL011_ipsi['spike_times'][k]<times+60))]   

for times in LL011_contra['seizure_start']:
    for k in range(len(LL011_contra['templates'])):
        LL011_contra['spike_times'][k] = LL011_contra['spike_times'][k][np.invert(np.logical_and(LL011_contra['spike_times'][k]>times, LL011_contra['spike_times'][k]<times+60))]   

for times in LL012_ipsi['seizure_start']:
    for k in range(len(LL012_ipsi['templates'])):
        LL012_ipsi['spike_times'][k] = LL012_ipsi['spike_times'][k][np.invert(np.logical_and(LL012_ipsi['spike_times'][k]>times, LL012_ipsi['spike_times'][k]<times+60))]   

for times in LL012_contra['seizure_start']:
    for k in range(len(LL012_contra['templates'])):
        LL012_contra['spike_times'][k] = LL012_contra['spike_times'][k][np.invert(np.logical_and(LL012_contra['spike_times'][k]>times, LL012_contra['spike_times'][k]<times+60))]   

for times in LL017['seizure_start']:
    for k in range(len(LL017['templates'])):
        LL017['spike_times'][k] = LL017['spike_times'][k][np.invert(np.logical_and(LL017['spike_times'][k]>times, LL017['spike_times'][k]<times+60))]   

for times in LL023['seizure_start']:
    for k in range(len(LL023['templates'])):
        LL023['spike_times'][k] = LL023['spike_times'][k][np.invert(np.logical_and(LL023['spike_times'][k]>times, LL023['spike_times'][k]<times+60))]   

LL001['ISI'] = [None] * len(LL001['templates'])
LL004['ISI'] = [None] * len(LL004['templates'])
LL005['ISI'] = [None] * len(LL005['templates'])
LL007['ISI'] = [None] * len(LL007['templates'])
LL009['ISI'] = [None] * len(LL009['templates'])
LL011_ipsi['ISI'] = [None] * len(LL011_ipsi['templates'])
LL011_contra['ISI'] = [None] * len(LL011_contra['templates'])
LL012_ipsi['ISI'] = [None] * len(LL012_ipsi['templates'])
LL012_contra['ISI'] = [None] * len(LL012_contra['templates'])
LL017['ISI'] = [None] * len(LL017['templates'])
LL023['ISI'] = [None] * len(LL023['templates'])

LL001['median_ISI'] = np.zeros((len(LL001['templates'])))
LL004['median_ISI'] = np.zeros((len(LL004['templates'])))
LL005['median_ISI'] = np.zeros((len(LL005['templates'])))
LL007['median_ISI'] = np.zeros((len(LL007['templates'])))
LL009['median_ISI'] = np.zeros((len(LL009['templates'])))
LL011_ipsi['median_ISI'] = np.zeros((len(LL011_ipsi['templates'])))
LL011_contra['median_ISI'] = np.zeros((len(LL011_contra['templates'])))
LL012_ipsi['median_ISI'] = np.zeros((len(LL012_ipsi['templates'])))
LL012_contra['median_ISI'] = np.zeros((len(LL012_contra['templates'])))
LL017['median_ISI'] = np.zeros((len(LL017['templates'])))
LL023['median_ISI'] = np.zeros((len(LL023['templates'])))

for i in range(len(LL001['ISI'])):
    LL001['ISI'][i] = np.diff(LL001['spike_times'][i])
    LL001['ISI'][i][LL001['ISI'][i]>60] = np.nan
    LL001['median_ISI'][i] = np.nanmedian(LL001['ISI'][i])

for i in range(len(LL004['ISI'])):
    LL004['ISI'][i] = np.diff(LL004['spike_times'][i])
    LL004['ISI'][i][LL004['ISI'][i]>60] = np.nan
    LL004['median_ISI'][i] = np.nanmedian(LL004['ISI'][i])

for i in range(len(LL005['ISI'])):
    LL005['ISI'][i] = np.diff(LL005['spike_times'][i])
    LL005['ISI'][i][LL005['ISI'][i]>60] = np.nan
    LL005['median_ISI'][i] = np.nanmedian(LL005['ISI'][i])
    
for i in range(len(LL007['ISI'])):
    LL007['ISI'][i] = np.diff(LL007['spike_times'][i])
    LL007['ISI'][i][LL007['ISI'][i]>60] = np.nan
    LL007['median_ISI'][i] = np.nanmedian(LL007['ISI'][i])
    
for i in range(len(LL009['ISI'])):
    LL009['ISI'][i] = np.diff(LL009['spike_times'][i])
    LL009['ISI'][i][LL009['ISI'][i]>60] = np.nan
    LL009['median_ISI'][i] = np.nanmedian(LL009['ISI'][i])
    
for i in range(len(LL011_ipsi['ISI'])):
    LL011_ipsi['ISI'][i] = np.diff(LL011_ipsi['spike_times'][i])
    LL011_ipsi['ISI'][i][LL011_ipsi['ISI'][i]>60] = np.nan
    LL011_ipsi['median_ISI'][i] = np.nanmedian(LL011_ipsi['ISI'][i])
    
for i in range(len(LL011_contra['ISI'])):
    LL011_contra['ISI'][i] = np.diff(LL011_contra['spike_times'][i])
    LL011_contra['ISI'][i][LL011_contra['ISI'][i]>60] = np.nan
    LL011_contra['median_ISI'][i] = np.nanmedian(LL011_contra['ISI'][i])
    
for i in range(len(LL012_ipsi['ISI'])):
    LL012_ipsi['ISI'][i] = np.diff(LL012_ipsi['spike_times'][i])
    LL012_ipsi['ISI'][i][LL012_ipsi['ISI'][i]>60] = np.nan
    LL012_ipsi['median_ISI'][i] = np.nanmedian(LL012_ipsi['ISI'][i])
    
for i in range(len(LL012_contra['ISI'])):
    LL012_contra['ISI'][i] = np.diff(LL012_contra['spike_times'][i])
    LL012_contra['ISI'][i][LL012_contra['ISI'][i]>60] = np.nan
    LL012_contra['median_ISI'][i] = np.nanmedian(LL012_contra['ISI'][i])

for i in range(len(LL017['ISI'])):
    LL017['ISI'][i] = np.diff(LL017['spike_times'][i])
    LL017['ISI'][i][LL017['ISI'][i]>60] = np.nan
    LL017['median_ISI'][i] = np.nanmedian(LL017['ISI'][i])
    
for i in range(len(LL023['ISI'])):
    LL023['ISI'][i] = np.diff(LL023['spike_times'][i])
    LL023['ISI'][i][LL023['ISI'][i]>60] = np.nan
    LL023['median_ISI'][i] = np.nanmedian(LL023['ISI'][i])
    

median_ISI_all = np.concatenate((LL001['median_ISI'],
                                 LL004['median_ISI'],
                                 LL005['median_ISI'],
                                 LL007['median_ISI'],
                                 LL009['median_ISI'],
                                 LL011_ipsi['median_ISI'],
                                 LL011_contra['median_ISI'],
                                 LL012_ipsi['median_ISI'],
                                 LL012_contra['median_ISI'],
                                 LL017['median_ISI'],
                                 LL023['median_ISI']
                                ))
    

np.save(savepath, median_ISI_all)
    
