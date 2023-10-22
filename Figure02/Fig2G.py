
import sys
sys.path.insert(1, r'E:\Manuscript_analysis_files')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import scipy.stats

from vies.parse.spike2 import load_spikematfile_bins

savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\output\Fig2G.png'

## NEURON CLUSTER 1
lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp1\Neuron1_Seizure1.mat'
time = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[0]
exp1_seizure1_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp1_seizure1_neuron1[0:24])
base_std = np.std(exp1_seizure1_neuron1[0:24])
p_exp1_seizure1_neuron1 = ((exp1_seizure1_neuron1[26]-base_mean) / base_mean) * 100
exp1_seizure1_neuron1 = (exp1_seizure1_neuron1 - base_mean) / base_std
exp1_seizure1_neuron1 = exp1_seizure1_neuron1[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp1\Neuron1_Seizure2.mat'
exp1_seizure2_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp1_seizure2_neuron1[0:24])
base_std = np.std(exp1_seizure2_neuron1[0:24])
p_exp1_seizure2_neuron1 = ((exp1_seizure2_neuron1[26]-base_mean) / base_mean) * 100
exp1_seizure2_neuron1 = (exp1_seizure2_neuron1 - base_mean) / base_std
exp1_seizure2_neuron1 = exp1_seizure2_neuron1[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp1\Neuron1_Seizure3.mat'
exp1_seizure3_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp1_seizure3_neuron1[0:24])
base_std = np.std(exp1_seizure3_neuron1[0:24])
p_exp1_seizure3_neuron1 = ((exp1_seizure3_neuron1[26]-base_mean) / base_mean) * 100
exp1_seizure3_neuron1 = (exp1_seizure3_neuron1 - base_mean) / base_std
exp1_seizure3_neuron1 = exp1_seizure3_neuron1[26]

e_cluster1=np.array([p_exp1_seizure1_neuron1, p_exp1_seizure2_neuron1, p_exp1_seizure3_neuron1])
cluster1=np.array([exp1_seizure1_neuron1,exp1_seizure2_neuron1, exp1_seizure3_neuron1])


## NEURON CLUSTER 2
lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp1\Neuron2_Seizure2.mat'
exp1_seizure2_neuron2 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp1_seizure2_neuron2[0:24])
base_std = np.std(exp1_seizure2_neuron2[0:24])
p_exp1_seizure2_neuron2 = ((exp1_seizure2_neuron2[26]-base_mean) / base_mean) * 100
exp1_seizure2_neuron2 = (exp1_seizure2_neuron2 - base_mean) / base_std
exp1_seizure2_neuron2 = exp1_seizure2_neuron2[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp1\Neuron2_Seizure3.mat'
exp1_seizure3_neuron2 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp1_seizure3_neuron2[0:24])
base_std = np.std(exp1_seizure3_neuron2[0:24])
p_exp1_seizure3_neuron2 = ((exp1_seizure3_neuron2[26]-base_mean) / base_mean) * 100
exp1_seizure3_neuron2 = (exp1_seizure3_neuron2 - base_mean) / base_std
exp1_seizure3_neuron2 = exp1_seizure3_neuron2[26]

e_cluster2=np.array([float('nan'), p_exp1_seizure2_neuron2, p_exp1_seizure3_neuron2])
cluster2=np.array([float('nan'), exp1_seizure2_neuron2, exp1_seizure3_neuron2])

## NEURON CLUSTER 3
lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp6\Neuron1_Seizure1.mat'
exp6_seizure1_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp6_seizure1_neuron1[0:24])
base_std = np.std(exp6_seizure1_neuron1[0:24])
p_exp6_seizure1_neuron1 = ((exp6_seizure1_neuron1[26]-base_mean) / base_mean) * 100
exp6_seizure1_neuron1 = (exp6_seizure1_neuron1 - base_mean) / base_std
exp6_seizure1_neuron1 = exp6_seizure1_neuron1[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp6\Neuron1_Seizure2.mat'
exp6_seizure2_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp6_seizure2_neuron1[0:24])
base_std = np.std(exp6_seizure2_neuron1[0:24])
p_exp6_seizure2_neuron1 = ((exp6_seizure2_neuron1[26]-base_mean) / base_mean) * 100
exp6_seizure2_neuron1 = (exp6_seizure2_neuron1 - base_mean) / base_std
exp6_seizure2_neuron1 = exp6_seizure2_neuron1[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp6\Neuron1_Seizure3.mat'
exp6_seizure3_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp6_seizure3_neuron1[0:24])
base_std = np.std(exp6_seizure3_neuron1[0:24])
p_exp6_seizure3_neuron1 = ((exp6_seizure3_neuron1[26]-base_mean) / base_mean) * 100
exp6_seizure3_neuron1 = (exp6_seizure3_neuron1 - base_mean) / base_std
exp6_seizure3_neuron1 = exp6_seizure3_neuron1[26]

e_cluster3=np.array([p_exp6_seizure1_neuron1, p_exp6_seizure2_neuron1, p_exp6_seizure3_neuron1])
cluster3=np.array([exp6_seizure1_neuron1, exp6_seizure2_neuron1, exp6_seizure3_neuron1])

## NEURON CLUSTER 4
lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp7\Neuron1_Seizure1.mat'
exp7_seizure1_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp7_seizure1_neuron1[0:24])
base_std = np.std(exp7_seizure1_neuron1[0:24])
p_exp7_seizure1_neuron1 = ((exp7_seizure1_neuron1[26]-base_mean) / base_mean) * 100
exp7_seizure1_neuron1 = (exp7_seizure1_neuron1 - base_mean) / base_std
exp7_seizure1_neuron1 = exp7_seizure1_neuron1[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp7\Neuron1_Seizure2.mat'
exp7_seizure2_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp7_seizure2_neuron1[0:24])
base_std = np.std(exp7_seizure2_neuron1[0:24])
p_exp7_seizure2_neuron1 = ((exp7_seizure2_neuron1[26]-base_mean) / base_mean) * 100
exp7_seizure2_neuron1 = (exp7_seizure2_neuron1 - base_mean) / base_std
exp7_seizure2_neuron1 = exp7_seizure2_neuron1[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp7\Neuron1_Seizure3.mat'
exp7_seizure3_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp7_seizure3_neuron1[0:24])
base_std = np.std(exp7_seizure3_neuron1[0:24])
p_exp7_seizure3_neuron1 = ((exp7_seizure3_neuron1[26]-base_mean) / base_mean) * 100
exp7_seizure3_neuron1 = (exp7_seizure3_neuron1 - base_mean) / base_std
exp7_seizure3_neuron1 = exp7_seizure3_neuron1[26]

e_cluster4=np.array([p_exp7_seizure1_neuron1, p_exp7_seizure2_neuron1, p_exp7_seizure3_neuron1])
cluster4=np.array([exp7_seizure1_neuron1, exp7_seizure2_neuron1, exp7_seizure3_neuron1])

## NEURON CLUSTER 5
lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp8\Neuron1_Seizure1.mat'
exp8_seizure1_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp8_seizure1_neuron1[0:24])
base_std = np.std(exp8_seizure1_neuron1[0:24])
p_exp8_seizure1_neuron1 = ((exp8_seizure1_neuron1[26]-base_mean) / base_mean) * 100
exp8_seizure1_neuron1 = (exp8_seizure1_neuron1 - base_mean) / base_std
exp8_seizure1_neuron1 = exp8_seizure1_neuron1[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp8\Neuron1_Seizure2.mat'
exp8_seizure2_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp8_seizure2_neuron1[0:24])
base_std = np.std(exp8_seizure2_neuron1[0:24])
p_exp8_seizure2_neuron1 = ((exp8_seizure2_neuron1[26]-base_mean) / base_mean) * 100
exp8_seizure2_neuron1 = (exp8_seizure2_neuron1 - base_mean) / base_std
exp8_seizure2_neuron1 = exp8_seizure2_neuron1[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp8\Neuron1_Seizure3.mat'
exp8_seizure3_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp8_seizure3_neuron1[0:24])
base_std = np.std(exp8_seizure3_neuron1[0:24])
p_exp8_seizure3_neuron1 = ((exp8_seizure3_neuron1[26]-base_mean) / base_mean) * 100
exp8_seizure3_neuron1 = (exp8_seizure3_neuron1 - base_mean) / base_std
exp8_seizure3_neuron1 = exp8_seizure3_neuron1[26]

e_cluster5=np.array([p_exp8_seizure1_neuron1, p_exp8_seizure2_neuron1, p_exp8_seizure3_neuron1])
cluster5=np.array([exp8_seizure1_neuron1, exp8_seizure2_neuron1, exp8_seizure3_neuron1])

## NEURON CLUSTER 6
lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp8\Neuron2_Seizure1.mat'
exp8_seizure1_neuron2 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp8_seizure1_neuron2[0:24])
base_std = np.std(exp8_seizure1_neuron2[0:24])
p_exp8_seizure1_neuron2 = ((exp8_seizure1_neuron2[26]-base_mean) / base_mean) * 100
exp8_seizure1_neuron2 = (exp8_seizure1_neuron2 - base_mean) / base_std
exp8_seizure1_neuron2 = exp8_seizure1_neuron2[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp8\Neuron2_Seizure2.mat'
exp8_seizure2_neuron2 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp8_seizure2_neuron2[0:24])
base_std = np.std(exp8_seizure2_neuron2[0:24])
p_exp8_seizure2_neuron2 = ((exp8_seizure2_neuron2[26]-base_mean) / base_mean) * 100
exp8_seizure2_neuron2 = (exp8_seizure2_neuron2 - base_mean) / base_std
exp8_seizure2_neuron2 = exp8_seizure2_neuron2[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp8\Neuron2_Seizure3.mat'
exp8_seizure3_neuron2 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp8_seizure3_neuron2[0:24])
base_std = np.std(exp8_seizure3_neuron2[0:24])
p_exp8_seizure3_neuron2 = ((exp8_seizure3_neuron2[26]-base_mean) / base_mean) * 100
exp8_seizure3_neuron2 = (exp8_seizure3_neuron2 - base_mean) / base_std
exp8_seizure3_neuron2 = exp8_seizure3_neuron2[26]

e_cluster6=np.array([p_exp8_seizure1_neuron2, p_exp8_seizure2_neuron2, p_exp8_seizure3_neuron2])
cluster6=np.array([exp8_seizure1_neuron2, exp8_seizure2_neuron2, exp8_seizure3_neuron2])

## NEURON CLUSTER 7
lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp8\Neuron3_Seizure1.mat'
exp8_seizure1_neuron3 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp8_seizure1_neuron3[0:24])
base_std = np.std(exp8_seizure1_neuron3[0:24])
p_exp8_seizure1_neuron3 = ((exp8_seizure1_neuron3[26]-base_mean) / base_mean) * 100
exp8_seizure1_neuron3 = (exp8_seizure1_neuron3 - base_mean) / base_std
exp8_seizure1_neuron3 = exp8_seizure1_neuron3[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp8\Neuron3_Seizure2.mat'
exp8_seizure2_neuron3 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp8_seizure2_neuron3[0:24])
base_std = np.std(exp8_seizure2_neuron3[0:24])
p_exp8_seizure2_neuron3 = ((exp8_seizure2_neuron3[26]-base_mean) / base_mean) * 100
exp8_seizure2_neuron3 = (exp8_seizure2_neuron3 - base_mean) / base_std
exp8_seizure2_neuron3 = exp8_seizure2_neuron3[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp8\Neuron3_Seizure3.mat'
exp8_seizure3_neuron3 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp8_seizure3_neuron3[0:24])
base_std = np.std(exp8_seizure3_neuron3[0:24])
p_exp8_seizure3_neuron3 = ((exp8_seizure3_neuron3[26]-base_mean) / base_mean) * 100
exp8_seizure3_neuron3 = (exp8_seizure3_neuron3 - base_mean) / base_std
exp8_seizure3_neuron3 = exp8_seizure3_neuron3[26]

e_cluster7=np.array([p_exp8_seizure1_neuron3, p_exp8_seizure2_neuron3, p_exp8_seizure3_neuron3])
cluster7=np.array([exp8_seizure1_neuron3, exp8_seizure2_neuron3, exp8_seizure3_neuron3])

## NEURON CLUSTER 8
lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp9\Neuron1_Seizure1.mat'
exp9_seizure1_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp9_seizure1_neuron1[0:24])
base_std = np.std(exp9_seizure1_neuron1[0:24])
p_exp9_seizure1_neuron1 = ((exp9_seizure1_neuron1[26]-base_mean) / base_mean) * 100
exp9_seizure1_neuron1 = (exp9_seizure1_neuron1 - base_mean) / base_std
exp9_seizure1_neuron1 = exp9_seizure1_neuron1[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp9\Neuron1_Seizure2.mat'
exp9_seizure2_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp9_seizure2_neuron1[0:24])
base_std = np.std(exp9_seizure2_neuron1[0:24])
p_exp9_seizure2_neuron1 = ((exp9_seizure2_neuron1[26]-base_mean) / base_mean) * 100
exp9_seizure2_neuron1 = (exp9_seizure2_neuron1 - base_mean) / base_std
exp9_seizure2_neuron1 = exp9_seizure2_neuron1[26]

lc_data = r'E:\Manuscript_analysis_files\data\spike_data\Exp9\Neuron1_Seizure3.mat'
exp9_seizure3_neuron1 = load_spikematfile_bins(lc_data, 1, 0, 180, binsize = 5)[1]
base_mean = np.mean(exp9_seizure3_neuron1[0:24])
base_std = np.std(exp9_seizure3_neuron1[0:24])
p_exp9_seizure3_neuron1 = ((exp9_seizure3_neuron1[26]-base_mean) / base_mean) * 100
exp9_seizure3_neuron1 = (exp9_seizure3_neuron1 - base_mean) / base_std
exp9_seizure3_neuron1 = exp9_seizure3_neuron1[26]

e_cluster8=np.array([p_exp9_seizure1_neuron1, p_exp9_seizure2_neuron1, p_exp9_seizure3_neuron1])
cluster8=np.array([exp9_seizure1_neuron1, exp9_seizure2_neuron1, exp9_seizure3_neuron1])


### COLLECT
z_seizure_response=np.column_stack((cluster1, cluster2, cluster3, cluster4, cluster5, cluster6, cluster7, cluster8)) # z-scores of responses
e_seizure_response=np.column_stack((e_cluster1, e_cluster2, e_cluster3, e_cluster4, e_cluster5, e_cluster6, e_cluster7, e_cluster8)) # effect of responses in %

p_values = scipy.stats.norm.sf(abs(z_seizure_response))*2 # p-values for responses

index = ['Seizure 1', 'Seizure 2', 'Seizure 3']


columns = ['Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8']
df_neurons = pd.DataFrame(data=e_seizure_response, index=index, columns=columns)
new_df_neurons = pd.melt(df_neurons.reset_index(), id_vars='index',value_vars=columns)
new_df_neurons.rename(columns={'index': 'Seizure', 'variable': 'Neuron', 'value': 'Data'}, inplace=True)

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["figure.figsize"] = (12,12)
ax = sns.stripplot(data=new_df_neurons, x='Neuron', y='Data', size=20,clip_on=False, hue='Seizure')
sns.despine()
sns.set_context("poster")
sns.set_style("white")
ax.set_ylabel('Normalized firing frequency during seizure (%)', fontsize=28)
sns.set(font_scale = 3)
ax.set_ylim([-110, 500])
ax.set_xlim([-0.4,7.4])
ax.set_yticks([-100, 0, 100, 200, 300, 400, 500])
ax.hlines(0,-0.4, 7.4, colors='k', linewidth=4, linestyles='dashed', zorder=1)
ax.spines['left'].set_position(('data', -0.4))
ax.spines['bottom'].set_position(('data', -110))
ax.spines['left'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)
ax.legend(loc='upper left', ncol=1, prop={'size': 25}, markerscale=1.5, frameon=False)
ax.set_xticklabels(columns, rotation=45, fontsize=28)
ax.set_yticklabels([-100, 0, 100, 200, 300, 400, 500], fontsize=25)
plt.xlabel('')
plt.tight_layout()
plt.savefig(os.path.abspath(savepath), dpi=300)


