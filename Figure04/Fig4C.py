# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 10:53:38 2022

@author: User
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, r'E:\Manuscript_analysis_files')

from vies.parse.phy import phy_data
from itertools import combinations
import pandas as pd

savepath_distance_fig = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\output\Fig4C.png'

unit_srate = 30000

LL001 = {}
LL004 = {}
LL005 = {}
LL007 = {}
#LL009 = {}
#LL011_ipsi = {}
LL011_contra = {}
LL012_ipsi = {}
LL012_contra = {}
LL017 = {}
LL023 = {}


LL001['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL001_ranked.npy'
LL004['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL004_ranked.npy'
LL005['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL005_ranked.npy'
LL007['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL007_ranked.npy'
#LL009['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL009_ranked.npy'
#LL011_ipsi['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL011_ipsi_ranked.npy'
LL011_contra['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL011_contra_ranked.npy'
LL012_ipsi['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL012_ipsi_ranked.npy'
LL012_contra['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL012_contra_ranked.npy'
LL017['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL017_ranked.npy'
LL023['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad\LL023_ranked.npy'


LL001 = np.load(LL001['path'], allow_pickle='TRUE').item()
LL004 = np.load(LL004['path'], allow_pickle='TRUE').item()
LL005 = np.load(LL005['path'], allow_pickle='TRUE').item()
LL007 = np.load(LL007['path'], allow_pickle='TRUE').item()
#LL009 = np.load(LL009['path'], allow_pickle='TRUE').item()
#LL011_ipsi = np.load(LL011_ipsi['path'], allow_pickle='TRUE').item()
LL011_contra = np.load(LL011_contra['path'], allow_pickle='TRUE').item()
LL012_ipsi = np.load(LL012_ipsi['path'], allow_pickle='TRUE').item()
LL012_contra = np.load(LL012_contra['path'], allow_pickle='TRUE').item()
LL017 = np.load(LL017['path'], allow_pickle='TRUE').item()
LL023 = np.load(LL023['path'], allow_pickle='TRUE').item()

significant_parametric = np.concatenate((LL001['significant_parametric'],
                              LL004['significant_parametric'], 
                              LL005['significant_parametric'], 
                              LL007['significant_parametric'], 
                              LL011_contra['significant_parametric'], 
                              LL012_ipsi['significant_parametric'], 
                              LL012_contra['significant_parametric'], 
                              LL017['significant_parametric'], 
                              LL023['significant_parametric']))

latencies = np.concatenate((LL001['peak_latency'],
                              LL004['peak_latency'], 
                              LL005['peak_latency'], 
                              LL007['peak_latency'], 
                              LL011_contra['peak_latency'], 
                              LL012_ipsi['peak_latency'], 
                              LL012_contra['peak_latency'], 
                              LL017['peak_latency'], 
                              LL023['peak_latency']))

LL001['templates'] = [5, 6, 10, 48, 59, 91, 92, 101, 110, 114, 118, 125, 149, 164, 165, 176, 178, 180, 197]
LL004['templates'] = [18, 30, 40, 51, 52]
LL005['templates'] = [17, 20, 26, 36, 38, 40, 47, 53, 61, 63, 77, 84, 95, 100, 112, 120, 124, 128, 136, 146]
LL007['templates'] = [81, 107, 112, 135, 157, 159, 181]
LL011_contra['templates'] = [15, 33, 95, 104, 107, 140]
LL012_ipsi['templates'] = [122, 137]
LL012_contra['templates'] = [12, 65, 156, 166, 170, 200, 204, 218, 317, 318, 329, 348, 354]
LL017['templates'] = [276, 293, 300, 317]
LL023['templates'] = [4, 5, 9, 12, 33, 48, 49]

LL001['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
LL004['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL004-final.GUI'
LL005['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'
LL007['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'
#LL009['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL009-final.GUI'
#LL011_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_datai\LL011_ipsi-final.GUI'
LL011_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL011_contra-final.GUI'
LL012_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_ipsi-final.GUI'
LL012_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_contra-final.GUI'
LL017['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL017-final.GUI'
LL023['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL023-final.GUI'

def compute_elec_distances(phy_path, templates, unit_srate=30000):
    
    template_pairs = list(combinations(templates, 2))
    data = phy_data(phy_path, unit_srate)
    
    info = data.template_info
    
    indices = np.zeros((len(templates)))
    for i in range(len(templates)):
        indices[i] = np.where(info[:,1] == templates[i])[0]
        
    electrode_contacts = info
    
    coordinates = data.channel_coordinates
    
    distances = np.zeros((np.shape(template_pairs)[0]))
    for i in range(np.shape(template_pairs)[0]):
        temp1 = template_pairs[i][0]
        temp2 = template_pairs[i][1]

        elec1 = electrode_contacts[np.where(electrode_contacts[:,1] == temp1)[0], 0]
        elec2 = electrode_contacts[np.where(electrode_contacts[:,1] == temp2)[0], 0]

        
        x1 = coordinates[elec1, 0]
        x2 = coordinates[elec2, 0]
        y1 = coordinates[elec1, 1]
        y2 = coordinates[elec2, 1]

        x_dif = x1 - x2
    
        y_dif = y1 - y2
    
        if x_dif == 0 or y_dif == 0:
            distance = np.abs(np.max([x_dif, y_dif]))
        else:
            distance = np.sqrt(x_dif**2 + y_dif**2)
        
        distances[i] = distance
    
    return distances

LL001['distances'] = compute_elec_distances(LL001['phy_path'], LL001['templates'])
LL004['distances'] = compute_elec_distances(LL004['phy_path'], LL004['templates'])
LL005['distances'] = compute_elec_distances(LL005['phy_path'], LL005['templates'])
LL007['distances'] = compute_elec_distances(LL007['phy_path'], LL007['templates'])
LL011_contra['distances'] = compute_elec_distances(LL011_contra['phy_path'], LL011_contra['templates'])
LL012_ipsi['distances'] = compute_elec_distances(LL012_ipsi['phy_path'], LL012_ipsi['templates'])
LL012_contra['distances'] = compute_elec_distances(LL012_contra['phy_path'], LL012_contra['templates'])
LL017['distances'] = compute_elec_distances(LL017['phy_path'], LL017['templates'])
LL023['distances'] = compute_elec_distances(LL023['phy_path'], LL023['templates'])

distances = np.concatenate((LL001['distances'],
                              LL004['distances'], 
                              LL005['distances'], 
                              LL007['distances'], 
                              LL011_contra['distances'], 
                              LL012_ipsi['distances'], 
                              LL012_contra['distances'], 
                              LL017['distances'], 
                              LL023['distances']))

columns = ['Distances', 'Latencies', 'significant_parametric']
data = pd.DataFrame(data=np.array([distances, latencies, significant_parametric]).T, columns = columns)

levels = 75 # micrometers

n_levels = int(550/levels)

for level in range(n_levels):
    start = level * levels
    stop = start + levels
    for i in range(np.shape(data)[0]):
        if data['Distances'][i] >= start and data['Distances'][i] < stop:
            data['Distances'][i] = start + levels/2
    
data_sigificant = data[data['significant_parametric'] == 1].groupby(['Distances']).count()
data_all = data.groupby(['Distances']).count()

fraction_broad = data_sigificant/np.sum(data_all, axis = 0)



########### SHARP SHARP SHARP SHARP SHARP SHARP SHARP ##################
LL001 = {}
LL004 = {}
LL005 = {}
LL007 = {}
#LL009 = {}
#LL011_ipsi = {}
LL011_contra = {}
LL012_ipsi = {}
LL012_contra = {}
LL017 = {}
LL023 = {}


LL001['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL001_ranked.npy'
LL004['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL004_ranked.npy'
LL005['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL005_ranked.npy'
LL007['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL007_ranked.npy'
#LL009['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL009_ranked.npy'
#LL011_ipsi['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL011_ipsi_ranked.npy'
LL011_contra['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL011_contra_ranked.npy'
LL012_ipsi['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL012_ipsi_ranked.npy'
LL012_contra['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL012_contra_ranked.npy'
LL017['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL017_ranked.npy'
LL023['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp\LL023_ranked.npy'


LL001 = np.load(LL001['path'], allow_pickle='TRUE').item()
LL004 = np.load(LL004['path'], allow_pickle='TRUE').item()
LL005 = np.load(LL005['path'], allow_pickle='TRUE').item()
LL007 = np.load(LL007['path'], allow_pickle='TRUE').item()
#LL009 = np.load(LL009['path'], allow_pickle='TRUE').item()
#LL011_ipsi = np.load(LL011_ipsi['path'], allow_pickle='TRUE').item()
LL011_contra = np.load(LL011_contra['path'], allow_pickle='TRUE').item()
LL012_ipsi = np.load(LL012_ipsi['path'], allow_pickle='TRUE').item()
LL012_contra = np.load(LL012_contra['path'], allow_pickle='TRUE').item()
LL017 = np.load(LL017['path'], allow_pickle='TRUE').item()
LL023 = np.load(LL023['path'], allow_pickle='TRUE').item()

significant_parametric = np.concatenate((LL001['significant_parametric'],
                              LL004['significant_parametric'], 
                              LL005['significant_parametric'], 
                              LL007['significant_parametric'], 
                              LL011_contra['significant_parametric'], 
                              LL012_ipsi['significant_parametric'], 
                              LL012_contra['significant_parametric'], 
                              LL017['significant_parametric'], 
                              LL023['significant_parametric']))

latencies = np.concatenate((LL001['peak_latency'],
                              LL004['peak_latency'], 
                              LL005['peak_latency'], 
                              LL007['peak_latency'], 
                              LL011_contra['peak_latency'], 
                              LL012_ipsi['peak_latency'], 
                              LL012_contra['peak_latency'], 
                              LL017['peak_latency'], 
                              LL023['peak_latency']))

LL001['templates'] = [5, 6, 10, 48, 59, 91, 92, 101, 110, 114, 118, 125, 149, 164, 165, 176, 178, 180, 197]
LL004['templates'] = [18, 30, 40, 51, 52]
LL005['templates'] = [17, 20, 26, 36, 38, 40, 47, 53, 61, 63, 77, 84, 95, 100, 112, 120, 124, 128, 136, 146]
LL007['templates'] = [81, 107, 112, 135, 157, 159, 181]
LL011_contra['templates'] = [15, 33, 95, 104, 107, 140]
LL012_ipsi['templates'] = [122, 137]
LL012_contra['templates'] = [12, 65, 156, 166, 170, 200, 204, 218, 317, 318, 329, 348, 354]
LL017['templates'] = [276, 293, 300, 317]
LL023['templates'] = [4, 5, 9, 12, 33, 48, 49]

LL001['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
LL004['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL004-final.GUI'
LL005['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'
LL007['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'
#LL009['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL009-final.GUI'
#LL011_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL011_ipsi-final.GUI'
LL011_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL011_contra-final.GUI'
LL012_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_ipsi-final.GUI'
LL012_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_contra-final.GUI'
LL017['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL017-final.GUI'
LL023['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL023-final.GUI'

def compute_elec_distances(phy_path, templates, unit_srate=30000):
    
    template_pairs = list(combinations(templates, 2))
    data = phy_data(phy_path, unit_srate)
    
    info = data.template_info
    
    indices = np.zeros((len(templates)))
    for i in range(len(templates)):
        indices[i] = np.where(info[:,1] == templates[i])[0]
        
    electrode_contacts = info
    
    coordinates = data.channel_coordinates
    
    distances = np.zeros((np.shape(template_pairs)[0]))
    for i in range(np.shape(template_pairs)[0]):
        temp1 = template_pairs[i][0]
        temp2 = template_pairs[i][1]
        
        #print(temp1, temp2)
        elec1 = electrode_contacts[np.where(electrode_contacts[:,1] == temp1)[0], 0]
        elec2 = electrode_contacts[np.where(electrode_contacts[:,1] == temp2)[0], 0]
        
        #print(elec1, elec2)
        
        x1 = coordinates[elec1, 0]
        x2 = coordinates[elec2, 0]
        y1 = coordinates[elec1, 1]
        y2 = coordinates[elec2, 1]
        
        #print(x1, x2, y1, y2)
        
        x_dif = x1 - x2
    
        y_dif = y1 - y2
    
        if x_dif == 0 or y_dif == 0:
            distance = np.abs(np.max([x_dif, y_dif]))
        else:
            distance = np.sqrt(x_dif**2 + y_dif**2)
        
        distances[i] = distance
    
    return distances

LL001['distances'] = compute_elec_distances(LL001['phy_path'], LL001['templates'])
LL004['distances'] = compute_elec_distances(LL004['phy_path'], LL004['templates'])
LL005['distances'] = compute_elec_distances(LL005['phy_path'], LL005['templates'])
LL007['distances'] = compute_elec_distances(LL007['phy_path'], LL007['templates'])
LL011_contra['distances'] = compute_elec_distances(LL011_contra['phy_path'], LL011_contra['templates'])
LL012_ipsi['distances'] = compute_elec_distances(LL012_ipsi['phy_path'], LL012_ipsi['templates'])
LL012_contra['distances'] = compute_elec_distances(LL012_contra['phy_path'], LL012_contra['templates'])
LL017['distances'] = compute_elec_distances(LL017['phy_path'], LL017['templates'])
LL023['distances'] = compute_elec_distances(LL023['phy_path'], LL023['templates'])

distances = np.concatenate((LL001['distances'],
                              LL004['distances'], 
                              LL005['distances'], 
                              LL007['distances'], 
                              LL011_contra['distances'], 
                              LL012_ipsi['distances'], 
                              LL012_contra['distances'], 
                              LL017['distances'], 
                              LL023['distances']))


columns = ['Distances', 'Latencies', 'significant_parametric']
data = pd.DataFrame(data=np.array([distances, latencies, significant_parametric]).T, columns = columns)



levels = 75 # micrometers

n_levels = int(550/levels)

for level in range(n_levels):
    start = level * levels
    stop = start + levels
    for i in range(np.shape(data)[0]):
        if data['Distances'][i] >= start and data['Distances'][i] < stop:
            data['Distances'][i] = start + levels/2
    
data_sigificant = data[data['significant_parametric'] == 1].groupby(['Distances']).count()
data_all = data.groupby(['Distances']).count()

fraction_sharp = data_sigificant/np.sum(data_all, axis = 0)
data_significant_parametric_latencies = data[data['significant_parametric'] == 1].groupby(['Distances']).mean()





### BROAD PLOT EXAMPLE
x_broad = np.arange(1.3, 5.3, 1)
x_sharp = np.arange(1.7, 5.7, 1)

fig, ax = plt.subplots(1, figsize=(4,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
      
ax.bar(x_broad, fraction_broad.iloc[0:4,0]*100, width = 0.4, color='gray', label = 'Broad')
ax.bar(x_sharp, fraction_sharp.iloc[0:4,0]*100, width = 0.4, color='darkgray', label = 'Sharp')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.set_xticks([1.5, 2.5, 3.5, 4.5])
ax.set_xticklabels(['<75', '75-150', '150-225', '<225']) 
ax.set_xlabel('Electrode distance ($\mu$m)', fontsize=14, fontweight='bold')
ax.set_ylabel('% of all pairs\nrecorded at each distance', fontsize=14, fontweight='bold')
ax.legend(frameon=False, loc='upper right')
plt.tight_layout()       
plt.savefig(savepath_distance_fig, dpi=600)




