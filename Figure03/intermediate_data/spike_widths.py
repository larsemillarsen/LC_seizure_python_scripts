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
from vies.spike.extract_waveforms import extract_waveforms
from scipy.signal import detrend
from scipy.signal import find_peaks

savepath_file = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure03\intermediate_data\width_asymmetry.npy'

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

LL001['spike_times'] = [None] * len(LL001['templates'])
for i in range(len(LL001['templates'])):
    LL001['spike_times'][i] = LL001['data'].extract_spike_trains(LL001['templates'][i])
    for times in LL001['seizure_start']:
        LL001['spike_times'][i] = LL001['spike_times'][i][np.invert(np.logical_and(LL001['spike_times'][i]>times, LL001['spike_times'][i]<times+60))] 
             
LL004['spike_times'] = [None] * len(LL004['templates'])
for i in range(len(LL004['templates'])):
    LL004['spike_times'][i] = LL004['data'].extract_spike_trains(LL004['templates'][i])
    for times in LL004['seizure_start']:
        LL004['spike_times'][i] = LL004['spike_times'][i][np.invert(np.logical_and(LL004['spike_times'][i]>times, LL004['spike_times'][i]<times+60))] 
             
LL005['spike_times'] = [None] * len(LL005['templates'])
for i in range(len(LL005['templates'])):
    LL005['spike_times'][i] = LL005['data'].extract_spike_trains(LL005['templates'][i])
    for times in LL005['seizure_start']:
        LL005['spike_times'][i] = LL005['spike_times'][i][np.invert(np.logical_and(LL005['spike_times'][i]>times, LL005['spike_times'][i]<times+60))] 
             
LL007['spike_times'] = [None] * len(LL007['templates'])
for i in range(len(LL007['templates'])):
    LL007['spike_times'][i] = LL007['data'].extract_spike_trains(LL007['templates'][i])
    for times in LL007['seizure_start']:
        LL007['spike_times'][i] = LL007['spike_times'][i][np.invert(np.logical_and(LL007['spike_times'][i]>times, LL007['spike_times'][i]<times+60))] 
             
LL009['spike_times'] = [None] * len(LL009['templates'])
for i in range(len(LL009['templates'])):
    LL009['spike_times'][i] = LL009['data'].extract_spike_trains(LL009['templates'][i])
    for times in LL009['seizure_start']:
        LL009['spike_times'][i] = LL009['spike_times'][i][np.invert(np.logical_and(LL009['spike_times'][i]>times, LL009['spike_times'][i]<times+60))] 
             
LL011_ipsi['spike_times'] = [None] * len(LL011_ipsi['templates'])
for i in range(len(LL011_ipsi['templates'])):
    LL011_ipsi['spike_times'][i] = LL011_ipsi['data'].extract_spike_trains(LL011_ipsi['templates'][i])
    for times in LL011_ipsi['seizure_start']:
        LL011_ipsi['spike_times'][i] = LL011_ipsi['spike_times'][i][np.invert(np.logical_and(LL011_ipsi['spike_times'][i]>times, LL011_ipsi['spike_times'][i]<times+60))] 
             
LL011_contra['spike_times'] = [None] * len(LL011_contra['templates'])
for i in range(len(LL011_contra['templates'])):
    LL011_contra['spike_times'][i] = LL011_contra['data'].extract_spike_trains(LL011_contra['templates'][i])
    for times in LL011_contra['seizure_start']:
        LL011_contra['spike_times'][i] = LL011_contra['spike_times'][i][np.invert(np.logical_and(LL011_contra['spike_times'][i]>times, LL011_contra['spike_times'][i]<times+60))] 
             
LL012_ipsi['spike_times'] = [None] * len(LL012_ipsi['templates'])
for i in range(len(LL012_ipsi['templates'])):
    LL012_ipsi['spike_times'][i] = LL012_ipsi['data'].extract_spike_trains(LL012_ipsi['templates'][i])
    for times in LL012_ipsi['seizure_start']:
        LL012_ipsi['spike_times'][i] = LL012_ipsi['spike_times'][i][np.invert(np.logical_and(LL012_ipsi['spike_times'][i]>times, LL012_ipsi['spike_times'][i]<times+60))] 
             
LL012_contra['spike_times'] = [None] * len(LL012_contra['templates'])
for i in range(len(LL012_contra['templates'])):
    LL012_contra['spike_times'][i] = LL012_contra['data'].extract_spike_trains(LL012_contra['templates'][i])
    for times in LL012_contra['seizure_start']:
        LL012_contra['spike_times'][i] = LL012_contra['spike_times'][i][np.invert(np.logical_and(LL012_contra['spike_times'][i]>times, LL012_contra['spike_times'][i]<times+60))] 
             
LL017['spike_times'] = [None] * len(LL017['templates'])
for i in range(len(LL017['templates'])):
    LL017['spike_times'][i] = LL017['data'].extract_spike_trains(LL017['templates'][i])
    for times in LL017['seizure_start']:
        LL017['spike_times'][i] = LL017['spike_times'][i][np.invert(np.logical_and(LL017['spike_times'][i]>times, LL017['spike_times'][i]<times+60))] 
             
LL023['spike_times'] = [None] * len(LL023['templates'])
for i in range(len(LL023['templates'])):
    LL023['spike_times'][i] = LL023['data'].extract_spike_trains(LL023['templates'][i])
    for times in LL023['seizure_start']:
        LL023['spike_times'][i] = LL023['spike_times'][i][np.invert(np.logical_and(LL023['spike_times'][i]>times, LL023['spike_times'][i]<times+60))] 
             
LL001['binary_path'] = r'E:\Data\LCSeizureData\LCSeizure_OpenEphys\13_LL001\analysis_108_400_4000\LL001.dat'
LL004['binary_path'] = r'E:\Data\LCSeizureData\LCSeizure_OpenEphys\11_LL004\analysis108_400_4000\LL004.dat'
LL005['binary_path'] = r'E:\Data\LCSeizureData\LCSeizure_OpenEphys\14_LL005\analysis_108_400_4000\LL005.dat'
LL007['binary_path'] = r'E:\Data\LCSeizureData\LCSeizure_OpenEphys\16_LL007\analysis_108_400_4000\LL007.dat'
LL009['binary_path'] = r'E:\Data\LCSeizureData\LCSeizure_OpenEphys\LL009\2_Ipsilateral\2020-10-19_18-20-39\Record Node 122\Analysis108\LL009.dat'
LL011_ipsi['binary_path'] = r'E:\Data\LCSeizureData\LCSeizure_OpenEphys\LL011\2_ipsilateral\2020-10-26_14-00-16\Record Node 131\analysis108\LL011_ipsi.dat'
LL011_contra['binary_path'] = r'E:\Data\LCSeizureData\LCSeizure_OpenEphys\LL011\3_contralateral\2020-10-26_16-45-38\Record Node 131\analysis108\LL011.dat'
LL012_ipsi['binary_path'] = r'E:\Data\LCSeizureData\LCSeizure_OpenEphys\LL012\1_ipsilateral\Record Node 138\analysis108\LL012_ipsi.dat'
LL012_contra['binary_path'] = r'E:\Data\LCSeizureData\LCSeizure_OpenEphys\LL012\3_contralateral\analysis108\LL012_contra.dat'
LL017['binary_path'] = r'E:\Data\LCSeizureData\LCSeizure_OpenEphys\LL017\Analysis108\LL017.dat'
LL023['binary_path'] = r'E:\Data\LCSeizureData\LCSeizure_OpenEphys\LL023\2020-12-18_13-43-19\Analysis108\LL023.dat'

### LL001
channels = LL001['data'].template_info
LL001['widths'] = np.zeros((len(LL001['templates'])))
LL001['peak_asymmetry'] = np.zeros((len(LL001['templates'])))
for i in range(np.shape(channels)[0]):
    if channels[i,1] in LL001['templates']:
        pass
    else:
        channels[i,:] = 0

channels = channels[~np.all(channels == 0, axis=1)] # removes all rows with 0s

#t = 0
for t in range(np.shape(channels)[0]):    
    
    if t == 0:
        waveforms = extract_waveforms(LL001['binary_path'], unit_srate, 32, LL001['spike_times'][t]*unit_srate, 5, 28)
    else:
        waveforms = extract_waveforms(LL001['binary_path'], unit_srate, 32, LL001['spike_times'][t]*unit_srate, 5, channels[t,0])
    
    for i in range(np.shape(waveforms)[1]):
        if np.argmin(waveforms[:,i], axis=0) < 70 or np.argmin(waveforms[:,i], axis=0) > 80:
            waveforms[:,i] = np.nan
            
    waveform_mean = np.nanmean(waveforms, axis=1)
    p1 = find_peaks(waveform_mean, height = np.mean(waveform_mean[0:55] + 6), distance=20, width = 4)
    n1 = find_peaks(waveform_mean*-1, height = np.mean(waveform_mean[0:55] + 100), distance=40, width = 4)
    
    plt.plot(waveform_mean)
    plt.scatter(p1[0], p1[1]['peak_heights'], c = 'r', marker = 'x')
    plt.scatter(n1[0], n1[1]['peak_heights']*-1, c = 'r', marker = 'x')
    
    savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\Wide_narrow\test_output\LL001'
    if os.path.isdir(savepath):
        #print("Saving Directory Exists")
        pass
    else:
        print("Created Saving Directory")
        os.makedirs(savepath)
                      
    #plt.savefig(savepath + '\\template' + str(LL001['templates'][t]) + '.png' , dpi=300)
    plt.close()
    
    a = p1[1]['peak_heights'][0]
    b = p1[1]['peak_heights'][1]
    
    LL001['widths'][t] = (p1[0][1] - p1[0][0]) / 30
    LL001['peak_asymmetry'][t] = (b-a)/(b+a)
    
    print(LL001['widths'][t], LL001['peak_asymmetry'][t])

### LL004
channels = LL004['data'].template_info
LL004['widths'] = np.zeros((len(LL004['templates'])))
LL004['peak_asymmetry'] = np.zeros((len(LL004['templates'])))
for i in range(np.shape(channels)[0]):
    if channels[i,1] in LL004['templates']:
        pass
    else:
        channels[i,:] = 0

#channels = channels[~np.all(channels == 0, axis=1)] # removes all rows with 0s

counter = 0
for t in range(np.shape(channels)[0]):    
    if counter < len(LL004['templates']):
        if channels[t,1] == LL004['templates'][counter]:
            waveforms = extract_waveforms(LL004['binary_path'], unit_srate, 32, LL004['spike_times'][counter]*unit_srate, 5, channels[t,0])
            
            for i in range(np.shape(waveforms)[1]):
                if np.argmin(waveforms[:,i], axis=0) < 70 or np.argmin(waveforms[:,i], axis=0) > 80:
                    waveforms[:,i] = np.nan
                    
            waveform_mean = np.nanmean(waveforms, axis=1)
            p1 = find_peaks(waveform_mean, height = np.mean(waveform_mean[0:55] + 6), distance=20, width = 4)
            n1 = find_peaks(waveform_mean*-1, height = np.mean(waveform_mean[0:55] + 100), distance=40, width = 4)
            
            plt.plot(waveform_mean)
            plt.scatter(p1[0], p1[1]['peak_heights'], c = 'r', marker = 'x')
            plt.scatter(n1[0], n1[1]['peak_heights']*-1, c = 'r', marker = 'x')
            
            savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\Wide_narrow\test_output\LL004'
            if os.path.isdir(savepath):
                #print("Saving Directory Exists")
                pass
            else:
                print("Created Saving Directory")
                os.makedirs(savepath)
                              
            #plt.savefig(savepath + '\\template' + str(LL004['templates'][counter]) + '.png' , dpi=300)
            plt.close()
            
            a = p1[1]['peak_heights'][0]
            b = p1[1]['peak_heights'][1]
            
            LL004['widths'][counter] = (p1[0][1] - p1[0][0]) / 30
            LL004['peak_asymmetry'][counter] = (b-a)/(b+a)
            
            print(LL004['widths'][counter], LL004['peak_asymmetry'][counter])
            counter = counter + 1
        else:
            pass
    else:
        pass


### LL005
channels = LL005['data'].template_info
LL005['widths'] = np.zeros((len(LL005['templates'])))
LL005['peak_asymmetry'] = np.zeros((len(LL005['templates'])))
for i in range(np.shape(channels)[0]):
    if channels[i,1] in LL005['templates']:
        pass
    else:
        channels[i,:] = 0

#channels = channels[~np.all(channels == 0, axis=1)] # removes all rows with 0s

counter = 0
for t in range(np.shape(channels)[0]):    
    if counter < len(LL005['templates']):
        if channels[t,1] == LL005['templates'][counter]:
            waveforms = extract_waveforms(LL005['binary_path'], unit_srate, 32, LL005['spike_times'][counter]*unit_srate, 5, channels[t,0])
            
            for i in range(np.shape(waveforms)[1]):
                if np.argmin(waveforms[:,i], axis=0) < 70 or np.argmin(waveforms[:,i], axis=0) > 80:
                    waveforms[:,i] = np.nan
                    
            waveform_mean = np.nanmean(waveforms, axis=1)
            p1 = find_peaks(waveform_mean, height = np.mean(waveform_mean[0:55] + 6), distance=20, width = 4)
            n1 = find_peaks(waveform_mean*-1, height = np.mean(waveform_mean[0:55] + 100), distance=40, width = 4)
            
            plt.plot(waveform_mean)
            plt.scatter(p1[0], p1[1]['peak_heights'], c = 'r', marker = 'x')
            plt.scatter(n1[0], n1[1]['peak_heights']*-1, c = 'r', marker = 'x')
            
            savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\Wide_narrow\test_output\LL005'
            if os.path.isdir(savepath):
                #print("Saving Directory Exists")
                pass
            else:
                print("Created Saving Directory")
                os.makedirs(savepath)
                              
            #plt.savefig(savepath + '\\template' + str(LL005['templates'][counter]) + '.png' , dpi=300)
            plt.close()
            
            a = p1[1]['peak_heights'][0]
            b = p1[1]['peak_heights'][1]
            
            LL005['widths'][counter] = (p1[0][1] - p1[0][0]) / 30
            LL005['peak_asymmetry'][counter] = (b-a)/(b+a)
            
            print(LL005['widths'][counter], LL005['peak_asymmetry'][counter])
            counter = counter + 1
        else:
            pass
    else:
        pass
    

### LL007
channels = LL007['data'].template_info
LL007['widths'] = np.zeros((len(LL007['templates'])))
LL007['peak_asymmetry'] = np.zeros((len(LL007['templates'])))
for i in range(np.shape(channels)[0]):
    if channels[i,1] in LL007['templates']:
        pass
    else:
        channels[i,:] = 0

#channels = channels[~np.all(channels == 0, axis=1)] # removes all rows with 0s

counter = 0
for t in range(np.shape(channels)[0]):    
    if counter < len(LL007['templates']):
        if channels[t,1] == LL007['templates'][counter]:
            waveforms = extract_waveforms(LL007['binary_path'], unit_srate, 32, LL007['spike_times'][counter]*unit_srate, 5, channels[t,0])
            
            for i in range(np.shape(waveforms)[1]):
                if np.argmin(waveforms[:,i], axis=0) < 70 or np.argmin(waveforms[:,i], axis=0) > 80:
                    waveforms[:,i] = np.nan
                    
            waveform_mean = np.nanmean(waveforms, axis=1)
            p1 = find_peaks(waveform_mean, height = np.mean(waveform_mean[0:55] + 6), distance=20, width = 4)
            n1 = find_peaks(waveform_mean*-1, height = np.mean(waveform_mean[0:55] + 100), distance=40, width = 4)
            
            plt.plot(waveform_mean)
            plt.scatter(p1[0], p1[1]['peak_heights'], c = 'r', marker = 'x')
            plt.scatter(n1[0], n1[1]['peak_heights']*-1, c = 'r', marker = 'x')
            
            savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\Wide_narrow\test_output\LL007'
            if os.path.isdir(savepath):
                #print("Saving Directory Exists")
                pass
            else:
                print("Created Saving Directory")
                os.makedirs(savepath)
                              
            #plt.savefig(savepath + '\\template' + str(LL007['templates'][counter]) + '.png' , dpi=300)
            plt.close()
            
            a = p1[1]['peak_heights'][0]
            b = p1[1]['peak_heights'][1]
            
            LL007['widths'][counter] = (p1[0][1] - p1[0][0]) / 30
            LL007['peak_asymmetry'][counter] = (b-a)/(b+a)
            
            print(LL007['widths'][counter], LL007['peak_asymmetry'][counter])
            counter = counter + 1
        else:
            pass
    else:
        pass
    

### LL009
channels = LL009['data'].template_info
LL009['widths'] = np.zeros((len(LL009['templates'])))
LL009['peak_asymmetry'] = np.zeros((len(LL009['templates'])))
for i in range(np.shape(channels)[0]):
    if channels[i,1] in LL009['templates']:
        pass
    else:
        channels[i,:] = 0

#channels = channels[~np.all(channels == 0, axis=1)] # removes all rows with 0s

counter = 0
for t in range(np.shape(channels)[0]):    
    if counter < len(LL009['templates']):
        if channels[t,1] == LL009['templates'][counter]:
            waveforms = extract_waveforms(LL009['binary_path'], unit_srate, 32, LL009['spike_times'][counter]*unit_srate, 5, channels[t,0])
            
            for i in range(np.shape(waveforms)[1]):
                if np.argmin(waveforms[:,i], axis=0) < 70 or np.argmin(waveforms[:,i], axis=0) > 80:
                    waveforms[:,i] = np.nan
                    
            waveform_mean = np.nanmean(waveforms, axis=1)
            p1 = find_peaks(waveform_mean, height = np.mean(waveform_mean[0:60] + 20), distance=20, width = 2)
            n1 = find_peaks(waveform_mean*-1, height = np.mean(waveform_mean[0:60] + 100), distance=40, width = 4)
            
            plt.plot(waveform_mean)
            plt.scatter(p1[0], p1[1]['peak_heights'], c = 'r', marker = 'x')
            plt.scatter(n1[0], n1[1]['peak_heights']*-1, c = 'r', marker = 'x')
            
            savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\Wide_narrow\test_output\LL009'
            if os.path.isdir(savepath):
                #print("Saving Directory Exists")
                pass
            else:
                print("Created Saving Directory")
                os.makedirs(savepath)
                              
            #plt.savefig(savepath + '\\template' + str(LL009['templates'][counter]) + '.png' , dpi=300)
            plt.close()
            
            a = p1[1]['peak_heights'][0]
            b = p1[1]['peak_heights'][1]
            
            LL009['widths'][counter] = (p1[0][1] - p1[0][0]) / 30
            LL009['peak_asymmetry'][counter] = (b-a)/(b+a)
            
            print(LL009['widths'][counter], LL009['peak_asymmetry'][counter])
            counter = counter + 1
        else:
            pass
    else:
        pass
    

### LL011_ipsi
channels = LL011_ipsi['data'].template_info
LL011_ipsi['widths'] = np.zeros((len(LL011_ipsi['templates'])))
LL011_ipsi['peak_asymmetry'] = np.zeros((len(LL011_ipsi['templates'])))
for i in range(np.shape(channels)[0]):
    if channels[i,1] in LL011_ipsi['templates']:
        pass
    else:
        channels[i,:] = 0

#channels = channels[~np.all(channels == 0, axis=1)] # removes all rows with 0s

counter = 0
for t in range(np.shape(channels)[0]):    
    if counter < len(LL011_ipsi['templates']):
        if channels[t,1] == LL011_ipsi['templates'][counter]:
            waveforms = extract_waveforms(LL011_ipsi['binary_path'], unit_srate, 32, LL011_ipsi['spike_times'][counter]*unit_srate, 5, channels[t,0])
            
            for i in range(np.shape(waveforms)[1]):
                if np.argmin(waveforms[:,i], axis=0) < 70 or np.argmin(waveforms[:,i], axis=0) > 80:
                    waveforms[:,i] = np.nan
                    
            waveform_mean = np.nanmean(waveforms, axis=1)
            p1 = find_peaks(waveform_mean, height = np.mean(waveform_mean[0:60]), distance=15, width = 4)
            n1 = find_peaks(waveform_mean*-1, height = np.mean(waveform_mean[0:60] + 100), distance=40, width = 4)
            
            plt.plot(waveform_mean)
            plt.scatter(p1[0], p1[1]['peak_heights'], c = 'r', marker = 'x')
            plt.scatter(n1[0], n1[1]['peak_heights']*-1, c = 'r', marker = 'x')
            
            #savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\Wide_narrow\test_output\LL011_ipsi'
            if os.path.isdir(savepath):
                #print("Saving Directory Exists")
                pass
            else:
                print("Created Saving Directory")
                os.makedirs(savepath)
                              
            plt.savefig(savepath + '\\template' + str(LL011_ipsi['templates'][counter]) + '.png' , dpi=300)
            plt.close()
            
            a = p1[1]['peak_heights'][0]
            b = p1[1]['peak_heights'][1]
            
            LL011_ipsi['widths'][counter] = (p1[0][1] - p1[0][0]) / 30
            LL011_ipsi['peak_asymmetry'][counter] = (b-a)/(b+a)
            
            print(LL011_ipsi['widths'][counter], LL011_ipsi['peak_asymmetry'][counter])
            counter = counter + 1
        else:
            pass
    else:
        pass
   
    
### LL011_contra
channels = LL011_contra['data'].template_info
LL011_contra['widths'] = np.zeros((len(LL011_contra['templates'])))
LL011_contra['peak_asymmetry'] = np.zeros((len(LL011_contra['templates'])))
for i in range(np.shape(channels)[0]):
    if channels[i,1] in LL011_contra['templates']:
        pass
    else:
        channels[i,:] = 0

#channels = channels[~np.all(channels == 0, axis=1)] # removes all rows with 0s

counter = 0
for t in range(np.shape(channels)[0]):    
    if counter < len(LL011_contra['templates']):
        if channels[t,1] == LL011_contra['templates'][counter]:
            waveforms = extract_waveforms(LL011_contra['binary_path'], unit_srate, 32, LL011_contra['spike_times'][counter]*unit_srate, 5, channels[t,0])
            
            for i in range(np.shape(waveforms)[1]):
                if np.argmin(waveforms[:,i], axis=0) < 70 or np.argmin(waveforms[:,i], axis=0) > 80:
                    waveforms[:,i] = np.nan
                    
            waveform_mean = np.nanmean(waveforms, axis=1)
            p1 = find_peaks(waveform_mean, height = np.mean(waveform_mean[0:60]+10), distance=15, width = 4)
            n1 = find_peaks(waveform_mean*-1, height = np.mean(waveform_mean[0:60] + 100), distance=40, width = 4)
            
            plt.plot(waveform_mean)
            plt.scatter(p1[0], p1[1]['peak_heights'], c = 'r', marker = 'x')
            plt.scatter(n1[0], n1[1]['peak_heights']*-1, c = 'r', marker = 'x')
            
            savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\Wide_narrow\test_output\LL011_contra'
            if os.path.isdir(savepath):
                #print("Saving Directory Exists")
                pass
            else:
                print("Created Saving Directory")
                os.makedirs(savepath)
                              
            #plt.savefig(savepath + '\\template' + str(LL011_contra['templates'][counter]) + '.png' , dpi=300)
            plt.close()
            
            a = p1[1]['peak_heights'][0]
            b = p1[1]['peak_heights'][1]
            
            LL011_contra['widths'][counter] = (p1[0][1] - p1[0][0]) / 30
            LL011_contra['peak_asymmetry'][counter] = (b-a)/(b+a)
            
            print(LL011_contra['widths'][counter], LL011_contra['peak_asymmetry'][counter])
            counter = counter + 1
        else:
            pass
    else:
        pass
    
    
### LL012_ipsi
channels = LL012_ipsi['data'].template_info
LL012_ipsi['widths'] = np.zeros((len(LL012_ipsi['templates'])))
LL012_ipsi['peak_asymmetry'] = np.zeros((len(LL012_ipsi['templates'])))
for i in range(np.shape(channels)[0]):
    if channels[i,1] in LL012_ipsi['templates']:
        pass
    else:
        channels[i,:] = 0

#channels = channels[~np.all(channels == 0, axis=1)] # removes all rows with 0s

counter = 0
for t in range(np.shape(channels)[0]):    
    if counter < len(LL012_ipsi['templates']):
        if channels[t,1] == LL012_ipsi['templates'][counter]:
            waveforms = extract_waveforms(LL012_ipsi['binary_path'], unit_srate, 32, LL012_ipsi['spike_times'][counter]*unit_srate, 5, channels[t,0])
            
            for i in range(np.shape(waveforms)[1]):
                if np.argmin(waveforms[:,i], axis=0) < 70 or np.argmin(waveforms[:,i], axis=0) > 80:
                    waveforms[:,i] = np.nan
                    
            waveform_mean = np.nanmean(waveforms, axis=1)
            p1 = find_peaks(waveform_mean, height = np.mean(waveform_mean[0:60]), distance=15, width = 4)
            n1 = find_peaks(waveform_mean*-1, height = np.mean(waveform_mean[0:60] + 100), distance=40, width = 4)
            
            plt.plot(waveform_mean)
            plt.scatter(p1[0], p1[1]['peak_heights'], c = 'r', marker = 'x')
            plt.scatter(n1[0], n1[1]['peak_heights']*-1, c = 'r', marker = 'x')
            
            savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\Wide_narrow\test_output\LL012_ipsi'
            if os.path.isdir(savepath):
                #print("Saving Directory Exists")
                pass
            else:
                print("Created Saving Directory")
                os.makedirs(savepath)
                              
            #plt.savefig(savepath + '\\template' + str(LL012_ipsi['templates'][counter]) + '.png' , dpi=300)
            plt.close()
            
            a = p1[1]['peak_heights'][0]
            b = p1[1]['peak_heights'][1]
            
            LL012_ipsi['widths'][counter] = (p1[0][1] - p1[0][0]) / 30
            LL012_ipsi['peak_asymmetry'][counter] = (b-a)/(b+a)
            
            print(LL012_ipsi['widths'][counter], LL012_ipsi['peak_asymmetry'][counter])
            counter = counter + 1
        else:
            pass
    else:
        pass
    
### LL012_contra
channels = LL012_contra['data'].template_info
LL012_contra['widths'] = np.zeros((len(LL012_contra['templates'])))
LL012_contra['peak_asymmetry'] = np.zeros((len(LL012_contra['templates'])))
for i in range(np.shape(channels)[0]):
    if channels[i,1] in LL012_contra['templates']:
        pass
    else:
        channels[i,:] = 0

#channels = channels[~np.all(channels == 0, axis=1)] # removes all rows with 0s

counter = 0
for t in range(np.shape(channels)[0]):    
    if counter < len(LL012_contra['templates']):
        if channels[t,1] == LL012_contra['templates'][counter]:
            waveforms = extract_waveforms(LL012_contra['binary_path'], unit_srate, 32, LL012_contra['spike_times'][counter]*unit_srate, 5, channels[t,0])
            
            for i in range(np.shape(waveforms)[1]):
                if np.argmin(waveforms[:,i], axis=0) < 70 or np.argmin(waveforms[:,i], axis=0) > 80:
                    waveforms[:,i] = np.nan
                    
            waveform_mean = np.nanmean(waveforms, axis=1)
            p1 = find_peaks(waveform_mean, height = np.mean(waveform_mean[0:60]+10), distance=15, width = 2)
            n1 = find_peaks(waveform_mean*-1, height = np.mean(waveform_mean[0:60] + 100), distance=40, width = 4)
            
            plt.plot(waveform_mean)
            plt.scatter(p1[0], p1[1]['peak_heights'], c = 'r', marker = 'x')
            plt.scatter(n1[0], n1[1]['peak_heights']*-1, c = 'r', marker = 'x')
            
            savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\Wide_narrow\test_output\LL012_contra'
            if os.path.isdir(savepath):
                #print("Saving Directory Exists")
                pass
            else:
                print("Created Saving Directory")
                os.makedirs(savepath)
                              
            #plt.savefig(savepath + '\\template' + str(LL012_contra['templates'][counter]) + '.png' , dpi=300)
            plt.close()
            
            a = p1[1]['peak_heights'][0]
            b = p1[1]['peak_heights'][1]
            
            LL012_contra['widths'][counter] = (p1[0][1] - p1[0][0]) / 30
            LL012_contra['peak_asymmetry'][counter] = (b-a)/(b+a)
            
            print(LL012_contra['widths'][counter], LL012_contra['peak_asymmetry'][counter])
            counter = counter + 1
        else:
            pass
    else:
        pass


### LL017
channels = LL017['data'].template_info
LL017['widths'] = np.zeros((len(LL017['templates'])))
LL017['peak_asymmetry'] = np.zeros((len(LL017['templates'])))
for i in range(np.shape(channels)[0]):
    if channels[i,1] in LL017['templates']:
        pass
    else:
        channels[i,:] = 0

#channels = channels[~np.all(channels == 0, axis=1)] # removes all rows with 0s

counter = 0
for t in range(np.shape(channels)[0]):    
    if counter < len(LL017['templates']):
        if channels[t,1] == LL017['templates'][counter]:
            waveforms = extract_waveforms(LL017['binary_path'], unit_srate, 32, LL017['spike_times'][counter]*unit_srate, 5, channels[t,0])
            
            for i in range(np.shape(waveforms)[1]):
                if np.argmin(waveforms[:,i], axis=0) < 70 or np.argmin(waveforms[:,i], axis=0) > 80:
                    waveforms[:,i] = np.nan
                    
            waveform_mean = np.nanmean(waveforms, axis=1)
            p1 = find_peaks(waveform_mean, height = np.mean(waveform_mean[0:60]+10), distance=15, width = 2)
            n1 = find_peaks(waveform_mean*-1, height = np.mean(waveform_mean[0:60] + 100), distance=40, width = 4)
            
            plt.plot(waveform_mean)
            plt.scatter(p1[0], p1[1]['peak_heights'], c = 'r', marker = 'x')
            plt.scatter(n1[0], n1[1]['peak_heights']*-1, c = 'r', marker = 'x')
            
            savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\Wide_narrow\test_output\LL017'
            if os.path.isdir(savepath):
                #print("Saving Directory Exists")
                pass
            else:
                print("Created Saving Directory")
                os.makedirs(savepath)
                              
            #plt.savefig(savepath + '\\template' + str(LL017['templates'][counter]) + '.png' , dpi=300)
            plt.close()
            
            a = p1[1]['peak_heights'][0]
            b = p1[1]['peak_heights'][1]
            
            LL017['widths'][counter] = (p1[0][1] - p1[0][0]) / 30
            LL017['peak_asymmetry'][counter] = (b-a)/(b+a)
            
            print(LL017['widths'][counter], LL017['peak_asymmetry'][counter])
            counter = counter + 1
        else:
            pass
    else:
        pass
    
    
### LL023
channels = LL023['data'].template_info
LL023['widths'] = np.zeros((len(LL023['templates'])))
LL023['peak_asymmetry'] = np.zeros((len(LL023['templates'])))
for i in range(np.shape(channels)[0]):
    if channels[i,1] in LL023['templates']:
        pass
    else:
        channels[i,:] = 0

#channels = channels[~np.all(channels == 0, axis=1)] # removes all rows with 0s

counter = 0
for t in range(np.shape(channels)[0]):    
    if counter < len(LL023['templates']):
        if channels[t,1] == LL023['templates'][counter]:
            waveforms = extract_waveforms(LL023['binary_path'], unit_srate, 32, LL023['spike_times'][counter]*unit_srate, 5, channels[t,0])
            
            for i in range(np.shape(waveforms)[1]):
                if np.argmin(waveforms[:,i], axis=0) < 70 or np.argmin(waveforms[:,i], axis=0) > 80:
                    waveforms[:,i] = np.nan
                    
            waveform_mean = np.nanmean(waveforms, axis=1)
            p1 = find_peaks(waveform_mean, height = np.mean(waveform_mean[0:60]+10), distance=15, width = 2)
            n1 = find_peaks(waveform_mean*-1, height = np.mean(waveform_mean[0:60] + 100), distance=40, width = 4)
            
            plt.plot(waveform_mean)
            plt.scatter(p1[0], p1[1]['peak_heights'], c = 'r', marker = 'x')
            plt.scatter(n1[0], n1[1]['peak_heights']*-1, c = 'r', marker = 'x')
            
            savepath = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\11_LC_neuron_characteristics\Wide_narrow\test_output\LL023'
            if os.path.isdir(savepath):
                #print("Saving Directory Exists")
                pass
            else:
                print("Created Saving Directory")
                os.makedirs(savepath)
                              
            #plt.savefig(savepath + '\\template' + str(LL023['templates'][counter]) + '.png' , dpi=300)
            plt.close()
            
            a = p1[1]['peak_heights'][0]
            b = p1[1]['peak_heights'][1]
            
            LL023['widths'][counter] = (p1[0][1] - p1[0][0]) / 30
            LL023['peak_asymmetry'][counter] = (b-a)/(b+a)
            
            print(LL023['widths'][counter], LL023['peak_asymmetry'][counter])
            counter = counter + 1
        else:
            pass
    else:
        pass
    
    
all_neurons_width = np.concatenate((
                                    LL001['widths'],
                                    LL004['widths'],
                                    LL005['widths'],
                                    LL007['widths'],
                                    LL009['widths'],
                                    LL011_ipsi['widths'],
                                    LL011_contra['widths'],
                                    LL012_ipsi['widths'],
                                    LL012_contra['widths'],
                                    LL017['widths'],
                                    LL023['widths'],                                    
                                    ))


all_neurons_peak_asymmetry = np.concatenate((
                                    LL001['peak_asymmetry'],
                                    LL004['peak_asymmetry'],
                                    LL005['peak_asymmetry'],
                                    LL007['peak_asymmetry'],
                                    LL009['peak_asymmetry'],
                                    LL011_ipsi['peak_asymmetry'],
                                    LL011_contra['peak_asymmetry'],
                                    LL012_ipsi['peak_asymmetry'],
                                    LL012_contra['peak_asymmetry'],
                                    LL017['peak_asymmetry'],
                                    LL023['peak_asymmetry'],                                    
                                    ))

all_templates = np.concatenate((
                                    LL001['templates'],
                                    LL004['templates'],
                                    LL005['templates'],
                                    LL007['templates'],
                                    LL009['templates'],
                                    LL011_ipsi['templates'],
                                    LL011_contra['templates'],
                                    LL012_ipsi['templates'],
                                    LL012_contra['templates'],
                                    LL017['templates'],
                                    LL023['templates'],                                    
                                    ))

plt.scatter(all_neurons_width, all_neurons_peak_asymmetry)

data = {}

data['templates'] = all_templates
data['widths'] = all_neurons_width
data['peak_asymmetry'] = all_neurons_peak_asymmetry

np.save(savepath_file, data)

