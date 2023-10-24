# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 17:23:36 2021

@author: llarsen
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, r'E:\OneDrive - UGent\python_functions')

import pandas as pd
import seaborn as sn

path = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\Fig2H\output_data'
savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\output\fig2H_sz1.png'
savepath2 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\output\fig2H_sz2.png'
savepath3 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\output\fig2H_sz3.png'

files = os.listdir(path)

df = np.load(path + '\\' + files[0])
for i in range(len(files)):
    if i == 0:
        df = np.load(path + '\\' + files[0])
    else:
        df = np.row_stack((df, np.load(path + '\\' + files[i])))

df = df[:,:-1].astype('float64')
    
change = np.zeros((3,3))
effect_1 = []
effect_2 = []
effect_3 = []
x1=np.ones((np.shape(df)[0]))
x2=np.ones((np.shape(df)[0])) * 20
x3=np.ones((np.shape(df)[0])) * 40
for i in range(np.shape(df)[0]):
    effect_1.append('empty')
    effect_2.append('empty')
    effect_3.append('empty')

## SZ1

p_threshold = 0.05
counter_inhibited = 0
counter_excited = 0
counter_nochange = 0
for y in range(np.shape(df)[0]):
    if df[y,4] < 0 and df[y,7] < p_threshold:
        counter_inhibited = counter_inhibited + 1
        effect_1[y] = 'inhibited (p<0.05)'
    elif df[y,4] > 0 and df[y,7] < p_threshold:
        counter_excited = counter_excited + 1
        effect_1[y] = 'excited (p<0.05)'
    elif df[y,7] > p_threshold:
        counter_nochange = counter_nochange + 1
        effect_1[y] = 'no change (p>0.05)'

    change[0,0] = counter_inhibited
    change[1,0] = counter_excited
    change[2,0] = counter_nochange


##SZ2
counter_inhibited = 0
counter_excited = 0
counter_nochange = 0
for y in range(np.shape(df)[0]):
    if df[y,5] < 0 and df[y,8] < p_threshold:
        counter_inhibited = counter_inhibited + 1
        effect_2[y] = 'inhibited (p<0.05)'
    elif df[y,5] > 0 and df[y,8] < p_threshold:
        counter_excited = counter_excited + 1
        effect_2[y] = 'excited (p<0.05)'
    elif df[y,8] > p_threshold:
        counter_nochange = counter_nochange + 1
        effect_2[y] = 'no change (p>0.05)'

    change[0,1] = counter_inhibited
    change[1,1] = counter_excited
    change[2,1] = counter_nochange 
    
##SZ3
counter_inhibited = 0
counter_excited = 0
counter_nochange = 0
for y in range(np.shape(df)[0]):
    if df[y,6] < 0 and df[y,9] < p_threshold:
        counter_inhibited = counter_inhibited + 1
        effect_3[y] = 'inhibited (p<0.05)'
    elif df[y,6] > 0 and df[y,9] < p_threshold:
        counter_excited = counter_excited + 1
        effect_3[y] = 'excited (p<0.05)'
    elif df[y,9] > p_threshold:
        counter_nochange = counter_nochange + 1
        effect_3[y] = 'no change (p>0.05)'

    change[0,2] = counter_inhibited
    change[1,2] = counter_excited
    change[2,2] = counter_nochange
    


headers = ['template ID', 'effect_sz1', 'effect_sz2', 'effect_sz3', 'z_sz1', 'z_sz2','z_sz3', 'p_sz1', 'p_sz2', 'p_sz3']
df2 = pd.DataFrame(data=df, columns=headers)
df2['effect_1'] = effect_1
df2['effect_2'] = effect_2
df2['effect_3'] = effect_3
df2.iloc[:,1] = y=df2.iloc[:,1] - 100
df2.iloc[:,2] = y=df2.iloc[:,2] - 100
df2.iloc[:,3] = y=df2.iloc[:,3] - 100



### Plotting ####
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(3,6))
plot = sn.stripplot(data=df2, x=x1, y='effect_sz1', hue='effect_1', palette=['b', 'r', 'g'], ax=ax, s=7)
ax.hlines(0, -0.15, 0.15, colors='k', linewidth=1, linestyles='dashed')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_linewidth(2)
ax.set_xlim(-0.15,0.15)
ax.set_yscale('symlog')
ax.set_yticks([-100, -10, 0, 10, 100, 1000])
ax.set_ylim([-120,1500])
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_yticklabels([-100, -10, 0, 10, 100, 1000])
ax.set_ylabel('LC firing (% from baseline mean)', fontsize=14)
ax.legend([],[],frameon=False)
plt.tight_layout()
plt.savefig(os.path.abspath(savepath), dpi=600)
#plt.close()

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(2,6))
plot = sn.stripplot(data=df2, x=x1, y='effect_sz2', hue='effect_2', palette=['b', 'r', 'g'], ax=ax, s=7)
ax.hlines(0, -0.15, 0.15, colors='k', linewidth=1, linestyles='dashed')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlim(-0.15,0.15)
ax.set_ylim([-120, 1500])
ax.set_yscale('symlog')
ax.set_yticks([])
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_yticklabels([],[])
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
ax.set_ylabel('')
ax.legend([],[],frameon=False)
plt.tight_layout()
plt.savefig(os.path.abspath(savepath2), dpi=600)
#plt.close()

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(1, figsize=(2,6))
plot = sn.stripplot(data=df2, x=x1, y='effect_sz3', hue='effect_3', palette=['b', 'r', 'g'], ax=ax, s=7)
ax.hlines(0, -0.15, 0.15, colors='k', linewidth=1, linestyles='dashed')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlim(-0.15,0.15)
ax.set_ylim([-120,1500])
ax.set_yscale('symlog')
ax.set_yticks([])
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
ax.set_ylabel('')
ax.legend([],[],frameon=False)
plt.tight_layout()
plt.savefig(os.path.abspath(savepath3), dpi=600)
#plt.close()
