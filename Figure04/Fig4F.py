# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 10:53:38 2022

@author: User
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, r'E:\OneDrive - UGent\python_functions')

from vies.parse.phy import phy_data
from itertools import combinations

from scipy.spatial.distance import jaccard

from scipy.stats import fisher_exact

unit_srate = 30000

savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\output\Seizure_correlation.png'

path_interictal_data = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\excited_coupling_broad.npy'

LL001 = {}
LL004 = {}
LL005 = {}
LL007 = {}
#LL009 = {}
#LL011_ipsi = {}
LL011_contra = {}
LL012_ipsi = {}
LL012_contra = {}
#LL017 = {}
LL023 = {}


LL001['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad_sz\LL001_ranked.npy'
LL004['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad_sz\LL004_ranked.npy'
LL005['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad_sz\LL005_ranked.npy'
LL007['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad_sz\LL007_ranked.npy'
#LL009['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad_sz\LL009.npy'
#LL011_ipsi['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad_sz\LL011_ipsi.npy'
#LL011_contra['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad_sz\LL011_contra_ranked.npy'
LL012_ipsi['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad_sz\LL012_ipsi_ranked.npy'
LL012_contra['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad_sz\LL012_contra_ranked.npy'
#LL017['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad_sz\LL017_ranked.npy'
LL023['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\broad_sz\LL023_ranked.npy'


LL001 = np.load(LL001['path'], allow_pickle='TRUE').item()
LL004 = np.load(LL004['path'], allow_pickle='TRUE').item()
LL005 = np.load(LL005['path'], allow_pickle='TRUE').item()
LL007 = np.load(LL007['path'], allow_pickle='TRUE').item()
#LL009 = np.load(LL009['path'], allow_pickle='TRUE').item()
#LL011_ipsi = np.load(LL011_ipsi['path'], allow_pickle='TRUE').item()
#LL011_contra = np.load(LL011_contra['path'], allow_pickle='TRUE').item()
LL012_ipsi = np.load(LL012_ipsi['path'], allow_pickle='TRUE').item()
LL012_contra = np.load(LL012_contra['path'], allow_pickle='TRUE').item()
#LL017 = np.load(LL017['path'], allow_pickle='TRUE').item()
LL023 = np.load(LL023['path'], allow_pickle='TRUE').item()

significant_parametric = np.concatenate((LL001['significant_parametric'],
                              LL004['significant_parametric'], 
                              LL005['significant_parametric'], 
                              LL007['significant_parametric'], 
                              LL012_ipsi['significant_parametric'], 
                              LL012_contra['significant_parametric'], 
                              #LL017['significant_parametric'], 
                              LL023['significant_parametric']))

latencies = np.concatenate((LL001['peak_latency'],
                              LL004['peak_latency'], 
                              LL005['peak_latency'], 
                              LL007['peak_latency'], 
                              LL012_ipsi['peak_latency'], 
                              LL012_contra['peak_latency'], 
                              #LL017['peak_latency'], 
                              LL023['peak_latency']))

LL001['templates'] = [5, 6, 10, 48, 59, 91, 92, 101, 110, 114, 118, 125, 149, 164, 165, 176, 178, 180, 197]
LL004['templates'] = [18, 30, 40, 51, 52]
LL005['templates'] = [17, 20, 26, 36, 38, 40, 47, 53, 61, 63, 77, 84, 95, 100, 112, 120, 124, 128, 136, 146]
LL007['templates'] = [81, 107, 112, 135, 157, 159, 181]
#LL011_contra['templates'] = [15, 33, 95, 104, 107, 140]
LL012_ipsi['templates'] = [137]
LL012_contra['templates'] = [12, 65, 156, 166, 170, 200, 204, 218, 317, 318, 329, 348, 354]
#LL017['templates'] = [276, 293, 300, 317]
LL023['templates'] = [4, 5, 9, 12, 33, 48, 49]

LL001['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
LL004['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL004-final.GUI'
LL005['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'
LL007['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'
#LL009['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL009-final.GUI'
#LL011_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL011_ipsi-final.GUI'
#LL011_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL011-final.GUI'
LL012_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_ipsi-final.GUI'
LL012_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_contra-final.GUI'
#LL017['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL017-final.GUI'
LL023['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL023-final.GUI'

LL001['pairs'] = list(combinations(LL001['templates'], 2))
LL001_data = phy_data(LL001['phy_path'], unit_srate)

LL001['info'] = LL001_data.template_info

indices = np.zeros((len(LL001['templates'])))
for i in range(len(LL001['templates'])):
    indices[i] = np.where(LL001['info'][:,1] == LL001['templates'][i])[0]
    
electrode_contacts = LL001['info']

coordinates = LL001_data.channel_coordinates

distances = np.zeros((np.shape(LL001['pairs'])[0]))
for i in range(np.shape(LL001['pairs'])[0]):
    temp1 = LL001['pairs'][i][0]
    temp2 = LL001['pairs'][i][1]
    elec1 = electrode_contacts[np.where(electrode_contacts[:,1] == temp1)[0], 0]
    elec2 = electrode_contacts[np.where(electrode_contacts[:,1] == temp2)[0], 0]
    
    x1 = coordinates[elec1, 0]
    x2 = coordinates[elec2, 0]
    y1 = coordinates[elec1, 0]
    y2 = coordinates[elec2, 0]
    
    x_dif = x1 - x2

    y_dif = y1 - y2

    if x_dif == 0 or y_dif == 0:
        distance = np.abs(np.max([x_dif, y_dif]))
    else:
        distance = np.sqrt(x_dif**2 + y_dif**2)
    
    distances[i] = distance

interictal_coupling = np.load(path_interictal_data, allow_pickle='TRUE')[:,0].astype(dtype='float')

comparison = np.column_stack((interictal_coupling, significant_parametric))




positives_2_p1 = 0
positives_2_n1 = 0

negatives_2_p1 = 0
negatives_2_n1 = 0

for i in range(np.shape(comparison)[0]):
    if comparison[i,0] == 1 and comparison[i,1] == 1:
        positives_2_p1 = positives_2_p1 + 1
    elif comparison[i,0] == 0 and comparison[i,1] == 1:
        positives_2_n1 = positives_2_n1 + 1        
    elif comparison[i,0] == 1 and comparison[i,1] == 0:
        negatives_2_p1 = negatives_2_p1 + 1    
    elif comparison[i,0] == 0 and comparison[i,1] == 0:
        negatives_2_n1 = negatives_2_n1 + 1

table_broad = np.array([[positives_2_p1, negatives_2_p1],
                  [positives_2_n1, negatives_2_n1]])

fishers_broad = fisher_exact(table_broad)



unit_srate = 30000

path_interictal_data = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\excited_coupling_sharp.npy'

LL001 = {}
LL004 = {}
LL005 = {}
LL007 = {}
#LL009 = {}
#LL011_ipsi = {}
LL011_contra = {}
LL012_ipsi = {}
LL012_contra = {}
#LL017 = {}
LL023 = {}


LL001['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL001_ranked.npy'
LL004['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL004_ranked.npy'
LL005['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL005_ranked.npy'
LL007['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL007_ranked.npy'
#LL009['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL009.npy'
#LL011_ipsi['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL011_ipsi.npy'
#LL011_contra['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL011_contra_ranked.npy'
LL012_ipsi['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL012_ipsi_ranked.npy'
LL012_contra['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL012_contra_ranked.npy'
#LL017['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL017_ranked.npy'
LL023['path'] = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure04\input\sharp_sz\LL023_ranked.npy'


LL001 = np.load(LL001['path'], allow_pickle='TRUE').item()
LL004 = np.load(LL004['path'], allow_pickle='TRUE').item()
LL005 = np.load(LL005['path'], allow_pickle='TRUE').item()
LL007 = np.load(LL007['path'], allow_pickle='TRUE').item()
#LL009 = np.load(LL009['path'], allow_pickle='TRUE').item()
#LL011_ipsi = np.load(LL011_ipsi['path'], allow_pickle='TRUE').item()
#LL011_contra = np.load(LL011_contra['path'], allow_pickle='TRUE').item()
LL012_ipsi = np.load(LL012_ipsi['path'], allow_pickle='TRUE').item()
LL012_contra = np.load(LL012_contra['path'], allow_pickle='TRUE').item()
#LL017 = np.load(LL017['path'], allow_pickle='TRUE').item()
LL023 = np.load(LL023['path'], allow_pickle='TRUE').item()

significant_parametric = np.concatenate((LL001['significant_parametric'],
                              LL004['significant_parametric'], 
                              LL005['significant_parametric'], 
                              LL007['significant_parametric'], 
                              LL012_ipsi['significant_parametric'], 
                              LL012_contra['significant_parametric'], 
                              #LL017['significant_parametric'], 
                              LL023['significant_parametric']))

latencies = np.concatenate((LL001['peak_latency'],
                              LL004['peak_latency'], 
                              LL005['peak_latency'], 
                              LL007['peak_latency'], 
                              LL012_ipsi['peak_latency'], 
                              LL012_contra['peak_latency'], 
                              #LL017['peak_latency'], 
                              LL023['peak_latency']))

LL001['templates'] = [5, 6, 10, 48, 59, 91, 92, 101, 110, 114, 118, 125, 149, 164, 165, 176, 178, 180, 197]
LL004['templates'] = [18, 30, 40, 51, 52]
LL005['templates'] = [17, 20, 26, 36, 38, 40, 47, 53, 61, 63, 77, 84, 95, 100, 112, 120, 124, 128, 136, 146]
LL007['templates'] = [81, 107, 112, 135, 157, 159, 181]
#LL011_contra['templates'] = [15, 33, 95, 104, 107, 140]
LL012_ipsi['templates'] = [137]
LL012_contra['templates'] = [12, 65, 156, 166, 170, 200, 204, 218, 317, 318, 329, 348, 354]
#LL017['templates'] = [276, 293, 300, 317]
LL023['templates'] = [4, 5, 9, 12, 33, 48, 49]

LL001['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
LL004['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL004-final.GUI'
LL005['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL005-final.GUI'
LL007['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL007-final.GUI'
#LL009['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL009-final.GUI'
#LL011_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL011_ipsi-final.GUI'
#LL011_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL011-final.GUI'
LL012_ipsi['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_ipsi-final.GUI'
LL012_contra['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL012_contra-final.GUI'
#LL017['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL017-final.GUI'
LL023['phy_path'] = r'E:\Manuscript_analysis_files\data\phy_data\LL023-final.GUI'

LL001['pairs'] = list(combinations(LL001['templates'], 2))
LL001_data = phy_data(LL001['phy_path'], unit_srate)

LL001['info'] = LL001_data.template_info

indices = np.zeros((len(LL001['templates'])))
for i in range(len(LL001['templates'])):
    indices[i] = np.where(LL001['info'][:,1] == LL001['templates'][i])[0]
    
electrode_contacts = LL001['info']

coordinates = LL001_data.channel_coordinates

distances = np.zeros((np.shape(LL001['pairs'])[0]))
for i in range(np.shape(LL001['pairs'])[0]):
    temp1 = LL001['pairs'][i][0]
    temp2 = LL001['pairs'][i][1]
    elec1 = electrode_contacts[np.where(electrode_contacts[:,1] == temp1)[0], 0]
    elec2 = electrode_contacts[np.where(electrode_contacts[:,1] == temp2)[0], 0]
    
    x1 = coordinates[elec1, 0]
    x2 = coordinates[elec2, 0]
    y1 = coordinates[elec1, 0]
    y2 = coordinates[elec2, 0]
    
    x_dif = x1 - x2

    y_dif = y1 - y2

    if x_dif == 0 or y_dif == 0:
        distance = np.abs(np.max([x_dif, y_dif]))
    else:
        distance = np.sqrt(x_dif**2 + y_dif**2)
    
    distances[i] = distance

interictal_coupling = np.load(path_interictal_data, allow_pickle='TRUE')[:,0].astype(dtype='float')

comparison = np.column_stack((interictal_coupling, significant_parametric))


positives_2_p1 = 0
positives_2_n1 = 0

negatives_2_p1 = 0
negatives_2_n1 = 0

for i in range(np.shape(comparison)[0]):
    if comparison[i,0] == 1 and comparison[i,1] == 1:
        positives_2_p1 = positives_2_p1 + 1
    elif comparison[i,0] == 0 and comparison[i,1] == 1:
        positives_2_n1 = positives_2_n1 + 1        
    elif comparison[i,0] == 1 and comparison[i,1] == 0:
        negatives_2_p1 = negatives_2_p1 + 1    
    elif comparison[i,0] == 0 and comparison[i,1] == 0:
        negatives_2_n1 = negatives_2_n1 + 1

table_sharp = np.array([[positives_2_p1, negatives_2_p1],
                  [positives_2_n1, negatives_2_n1]])

fishers_sharp = fisher_exact(table_sharp)



### BROAD PLOT EXAMPLE
x_broad = [4, 5]
x_sharp = [1, 2]

fig, ax = plt.subplots(1, figsize=(5,5))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
      
ax.bar(x_sharp, table_sharp[0,:], width = 0.6, color='gray', label = 'Correlated')
ax.bar(x_sharp, table_sharp[1,:], width = 0.6, color='darkgray', label = 'Non-correlated', bottom = table_sharp[0,:])

ax.bar(x_broad, table_broad[0,:], width = 0.6, color='gray') 
ax.bar(x_broad, table_broad[1,:], width = 0.6, color='darkgray', bottom = table_broad[0,:])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

#ax.set_yticks([1.5, 2.5, 3.5, 4.5]) 
#ax.set_yticklabels(['', 0.2, 0.4, 0.6], fontsize=12, fontweight='bold') 
ax.set_xticks([1, 2, 4, 5])
ax.set_xticklabels(['Sz correlated', 'Sz non-\ncorrelated', 'Sz correlated', 'Sz non-\ncorrelated'], rotation=45) 

ax.set_xlabel('Sharp                           Broad', fontsize=14, fontweight='bold')
ax.set_ylabel('# pairs', fontsize=14, fontweight='bold')
ax.set_ylim((0, 65))
ax.legend(frameon=False, loc='upper right')
plt.tight_layout()       
plt.savefig(savepath, dpi=600)

