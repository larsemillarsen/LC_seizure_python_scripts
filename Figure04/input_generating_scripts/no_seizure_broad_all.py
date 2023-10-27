# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 10:53:38 2022

@author: User
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, r'C:\Users\User\OneDrive - UGent\python_functions')

from vies.parse.phy import phy_data
from itertools import combinations

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


LL001['path'] = r'C:\Users\User\OneDrive - UGent\88_LCseizureproject\results\Report_figures\09_Cross_correlation\output\LL001\broad_no_sz\LL001_ranked.npy'
LL004['path'] = r'C:\Users\User\OneDrive - UGent\88_LCseizureproject\results\Report_figures\09_Cross_correlation\output\LL004\broad_no_sz\LL004_ranked.npy'
LL005['path'] = r'C:\Users\User\OneDrive - UGent\88_LCseizureproject\results\Report_figures\09_Cross_correlation\output\LL005\broad_no_sz\LL005_ranked.npy'
LL007['path'] = r'C:\Users\User\OneDrive - UGent\88_LCseizureproject\results\Report_figures\09_Cross_correlation\output\LL007\broad_no_sz\LL007_ranked.npy'
#LL009['path'] = r'C:\Users\User\OneDrive - UGent\88_LCseizureproject\results\Report_figures\09_Cross_correlation\output\LL009\broad_no_sz\LL009.npy'
#LL011_ipsi['path'] = r'C:\Users\User\OneDrive - UGent\88_LCseizureproject\results\Report_figures\09_Cross_correlation\output\LL011_ipsi\broad_no_sz\LL011_ipsi.npy'
LL011_contra['path'] = r'C:\Users\User\OneDrive - UGent\88_LCseizureproject\results\Report_figures\09_Cross_correlation\output\LL011_contra\broad_no_sz\LL011_contra_ranked.npy'
LL012_ipsi['path'] = r'C:\Users\User\OneDrive - UGent\88_LCseizureproject\results\Report_figures\09_Cross_correlation\output\LL012_ipsi\broad_no_sz\LL012_ipsi_ranked.npy'
LL012_contra['path'] = r'C:\Users\User\OneDrive - UGent\88_LCseizureproject\results\Report_figures\09_Cross_correlation\output\LL012_contra\broad_no_sz\LL012_contra_ranked.npy'
LL017['path'] = r'C:\Users\User\OneDrive - UGent\88_LCseizureproject\results\Report_figures\09_Cross_correlation\output\LL017\broad_no_sz\LL017_ranked.npy'
LL023['path'] = r'C:\Users\User\OneDrive - UGent\88_LCseizureproject\results\Report_figures\09_Cross_correlation\output\LL023\broad_no_sz\LL023_ranked.npy'


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
LL012_ipsi['templates'] = [137]
LL012_contra['templates'] = [12, 65, 156, 166, 170, 200, 204, 218, 317, 318, 329, 348, 354]
LL017['templates'] = [276, 293, 300, 317]
LL023['templates'] = [4, 5, 9, 12, 33, 48, 49]

LL001['phy_path'] = r'E:\LCSeizureData\LCSeizure_OpenEphys\13_LL001\analysis_108_400_4000\LL001\LL001-final.GUI'
LL004['phy_path'] = r'E:\LCSeizureData\LCSeizure_OpenEphys\11_LL004\analysis108_400_4000\LL004\LL004-final.GUI'
LL005['phy_path'] = r'E:\LCSeizureData\LCSeizure_OpenEphys\14_LL005\analysis_108_400_4000\LL005\LL005-final.GUI'
LL007['phy_path'] = r'E:\LCSeizureData\LCSeizure_OpenEphys\16_LL007\analysis_108_400_4000\LL007\LL007-final.GUI'
#LL009['phy_path'] = r'E:\LCSeizureData\LCSeizure_OpenEphys\LL009\2_Ipsilateral\2020-10-19_18-20-39\Record Node 122\Analysis108\LL009\LL009-final.GUI'
#LL011_ipsi['phy_path'] = r'E:\LCSeizureData\LCSeizure_OpenEphys\LL011\2_ipsilateral\2020-10-26_14-00-16\Record Node 131\analysis108\LL011_ipsi\LL011_ipsi-final.GUI'
LL011_contra['phy_path'] = r'E:\LCSeizureData\LCSeizure_OpenEphys\LL011\3_contralateral\2020-10-26_16-45-38\Record Node 131\analysis108\LL011\LL011-final.GUI'
LL012_ipsi['phy_path'] = r'E:\LCSeizureData\LCSeizure_OpenEphys\LL012\1_ipsilateral\Record Node 138\analysis108\LL012_ipsi\LL012_ipsi-final.GUI'
LL012_contra['phy_path'] = r'E:\LCSeizureData\LCSeizure_OpenEphys\LL012\3_contralateral\analysis108\LL012_contra\LL012_contra-final.GUI'
LL017['phy_path'] = r'E:\LCSeizureData\LCSeizure_OpenEphys\LL017\Analysis108\LL017\LL017-final.GUI'
LL023['phy_path'] = r'E:\LCSeizureData\LCSeizure_OpenEphys\LL023\2020-12-18_13-43-19\Analysis108\LL023\LL023-final.GUI'

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


#LL001['data'] = phy_data(LL001['phy_path'], unit_srate)



'''
x_dif = coordinates[0,0] - coordinates[1,0]

y_dif = coordinates[0,1] - coordinates[1,1]

if x_dif == 0 or y_dif == 0:
    distance = np.abs(np.max([x_dif, y_dif]))
else:
    distance = np.sqrt(x_dif**2 + y_dif**2)




'''

