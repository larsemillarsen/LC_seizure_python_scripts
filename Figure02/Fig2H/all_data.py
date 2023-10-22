# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 12:15:00 2022

@author: llarsen
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sn

sys.path.insert(1, r'E:\OneDrive - UGent\python_functions')

directory = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\01_Fig_pinch_light\revision\non_LC_neurons_seizure\data'
save1 = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\01_Fig_pinch_light\revision\non_LC_neurons_seizure\sz1.jpg'
save2 = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\01_Fig_pinch_light\revision\non_LC_neurons_seizure\sz2.jpg'
save3 = r'E:\OneDrive - UGent\88_LCseizureproject\results\Report_figures\01_Fig_pinch_light\revision\non_LC_neurons_seizure\sz3.jpg'
files = os.listdir(directory)

no_neurons_sz1 = 0
no_neurons_sz2 = 0
no_neurons_sz3 = 0

for i in range(len(files)):
    file = directory + '\\' + files[i]
    data = np.load(file)
    if np.shape(data)[2] == 3:
        no_neurons_sz1 = no_neurons_sz1 + np.shape(data)[1]
        no_neurons_sz2 = no_neurons_sz2 + np.shape(data)[1]
        no_neurons_sz3 = no_neurons_sz3 + np.shape(data)[1]
    else:
        no_neurons_sz1 = no_neurons_sz1 + np.shape(data)[1]
        no_neurons_sz2 = no_neurons_sz2 + np.shape(data)[1]

data_sz1 = np.zeros((np.shape(data)[0], no_neurons_sz1))
data_sz2 = np.zeros((np.shape(data)[0], no_neurons_sz2))
data_sz3 = np.zeros((np.shape(data)[0], no_neurons_sz3))

neuron_counter_1 = 0
neuron_counter_2 = 0
for i in range(len(files)):
    file = directory + '\\' + files[i]
    data = np.load(file)
    
    if i == 0:
        if np.shape(data)[2] == 3:
            data_sz1 = data[:,:,0]
            data_sz2 = data[:,:,1]
            data_sz3 = data[:,:,2]
        else:
            data_sz1 = data[:,:,0]
            data_sz2 = data[:,:,1]            
    else:
        if np.shape(data)[2] == 3:
            data_sz1 = np.column_stack((data_sz1, data[:,:,0]))
            data_sz2 = np.column_stack((data_sz2, data[:,:,1]))
            data_sz3 = np.column_stack((data_sz3, data[:,:,2]))
        else:
            data_sz1 = np.column_stack((data_sz1, data[:,:,0]))
            data_sz2 = np.column_stack((data_sz2, data[:,:,1]))


# z-score

sz1_std = np.nanstd(data_sz1[0:30], axis=0)
sz1_mean = np.nanmean(data_sz1[0:30], axis=0)    
z_sz1 = (data_sz1 - sz1_mean) / sz1_std
e_sz1 = ((data_sz1 - sz1_mean) / sz1_mean) * 100
sig_sz1 = np.zeros((no_neurons_sz1)) + 2
sig_sz1[z_sz1[32,:]>1.96] = 1
sig_sz1[z_sz1[32,:]<-1.96] = 0

sz2_std = np.nanstd(data_sz2[0:30], axis=0)
sz2_mean = np.nanmean(data_sz2[0:30], axis=0) 
z_sz2 = (data_sz2 - sz2_mean) / sz2_std
e_sz2 = ((data_sz2 - sz2_mean) / sz2_mean) * 100
sig_sz2 = np.zeros((no_neurons_sz2)) + 2
sig_sz2[z_sz2[32,:]>1.96] = 1
sig_sz2[z_sz2[32,:]<-1.96] = 0

sz3_std = np.nanstd(data_sz3[0:30], axis=0)
sz3_mean = np.nanmean(data_sz3[0:30], axis=0) 
z_sz3 = (data_sz3 - sz3_mean) / sz3_std
e_sz3 = ((data_sz3 - sz3_mean) / sz3_mean) * 100
sig_sz3 = np.zeros((no_neurons_sz3)) + 2
sig_sz3[z_sz3[32,:]>1.96] = 1
sig_sz3[z_sz3[32,:]<-1.96] = 0


x1=np.ones(no_neurons_sz1)
x2=np.ones(no_neurons_sz3)

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(3,6))
plot = sn.stripplot(x=x1, y=e_sz1[32,:], hue=sig_sz1, palette=['b', 'r', 'g'], ax=ax, s=7)
#plot = sn.swarmplot(data=df2, x=x2, y='sz2_effect', hue='effect', palette='bright', ax=ax)
#plot = sn.swarmplot(data=df2, x=x3, y='sz3_effect', hue='effect', palette='bright', ax=ax)

ax.hlines(0, -0.15, 0.15, colors='k', linewidth=1, linestyles='dashed')
#ax[0].vlines(25, -0.2, 9.2, colors='r', linewidth=2, linestyles='dashed')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_linewidth(2)
ax.set_xlim(-0.15,0.15)
ax.set_yscale('symlog')
ax.set_yticks([-100, -10, 0, 10, 100, 1000])
ax.set_ylim([-120,2500])
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_yticklabels([-100, -10, 0, 10, 100, 1000])
#ax.title.set_text('Template ' + cluster_id)
#ax.set_xticks(np.arange(0, 40.1, step=10))
#ax.set_xlabel('Time (s)', fontsize=14)
ax.set_ylabel('Non-LC firing (% from baseline mean)', fontsize=14)
ax.legend([],[],frameon=False)
plt.tight_layout()
plt.savefig(os.path.abspath(save1), dpi=600)
#plt.close()
#sn.swarmplot(data=df2, x=x, y='sz1_effect', hue='effect', palette='bright', )


#plt.yscale('symlog')
#plt.ylim([-100, 1400])


plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(2,6))
#plot = sn.swarmplot(data=df2, x=x1, y='sz1_effect', hue='effect', palette='bright', ax=ax)
plot = sn.stripplot(x=x1, y=e_sz2[32,:], hue=sig_sz2, palette=['b', 'r', 'g'], ax=ax, s=7)
#plot = sn.swarmplot(data=df2, x=x3, y='sz3_effect', hue='effect', palette='bright', ax=ax)

ax.hlines(0, -0.15, 0.15, colors='k', linewidth=1, linestyles='dashed')
#ax[0].vlines(25, -0.2, 9.2, colors='r', linewidth=2, linestyles='dashed')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlim(-0.15,0.15)
ax.set_ylim([-120, 2500])
ax.set_yscale('symlog')
ax.set_yticks([])
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_yticklabels([],[])
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
#ax.tick_params(left=False, bottom=False)
#ax.title.set_text('Template ' + cluster_id)
#ax.set_xticks(np.arange(0, 40.1, step=10))
#ax.set_xlabel('Time (s)', fontsize=14)
ax.set_ylabel('')
ax.legend([],[],frameon=False)
plt.tight_layout()
plt.savefig(os.path.abspath(save2), dpi=600)
plt.close()
#sn.swarmplot(data=df2, x=x, y='sz1_effect', hue='effect', palette='bright', )


#plt.yscale('symlog')
#plt.ylim([-100, 1400])

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(2,6))
#plot = sn.swarmplot(data=df2, x=x1, y='sz1_effect', hue='effect', palette='bright', ax=ax)
plot = sn.stripplot(x=x2, y=e_sz3[32,:], hue=sig_sz3, palette=['b', 'r', 'g'], ax=ax, s=7)
#plot = sn.swarmplot(data=df2, x=x3, y='sz3_effect', hue='effect', palette='bright', ax=ax)

ax.hlines(0, -0.15, 0.15, colors='k', linewidth=1, linestyles='dashed')
#ax[0].vlines(25, -0.2, 9.2, colors='r', linewidth=2, linestyles='dashed')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlim(-0.15,0.15)
ax.set_ylim([-120,2500])
ax.set_yscale('symlog')
ax.set_yticks([])
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
#ax.tick_params(left=False, bottom=False)
#ax.title.set_text('Template ' + cluster_id)
#ax.set_xticks(np.arange(0, 40.1, step=10))
#ax.set_xlabel('Time (s)', fontsize=14)
ax.set_ylabel('')
ax.legend([],[],frameon=False)
#ax.legend(bbox_to_anchor=(1.04, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.abspath(save3), dpi=600)
#plt.close()
#sn.swarmplot(data=df2, x=x, y='sz1_effect', hue='effect', palette='bright', )


