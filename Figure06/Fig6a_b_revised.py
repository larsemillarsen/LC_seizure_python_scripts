
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import os
from scipy.stats import pearsonr

import sys
sys.path.insert(1, r'E:\OneDrive - UGent\python_functions')

from vies.lfp.lfp_analysis import spectrogram
from vies.lfp.filter import butter_bandpass_filter, butter_bandstop_filter
from vies.parse.neuron import load_neuronfile
from vies.parse.phy import phy_data
from vies.spike.extract_frequency_bins import extract_frequency_bins

from scipy.stats import pearsonr, spearmanr
from vies.lfp.filter import movingaverage
from scipy.signal import correlate, correlation_lags
from f_cross_corr_function import eegpower_spikerate_crosscorr


## Params
lowcut_eeg=2
highcut_eeg=500
srate_eeg=30000
desired_eeg_srate=2000
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL001\Seizure1.mat'
#lc_data_inhibition = r'E:\Data\LCSeizureData\LCSeizureNIDAQ\01_Exp1\redone\try1\try3\matneurons\Neuron2_Seizure2.mat'
savepath = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure06\output\Fig6A.png'
savepath2 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure06\output\Fig6B.png'
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
seizure_start = 3782.72
spike_templates = [178, 164, 165]
unit_srate = 30000
LL001_spike_templates = [5, 6, 10, 48, 59, 91, 92, 101, 110, 114, 118, 125, 149, 164, 165, 176, 178, 180, 197]


## Load and preprocess data
time_eeg, hip_eeg = load_neuronfile(eeg_path,srate=30000,channel=1, gain=100, inputrange=20)
hip_eeg=butter_bandpass_filter(hip_eeg, lowcut_eeg, highcut_eeg, srate_eeg, order=2)
hip_eeg = signal.decimate(hip_eeg, int(srate_eeg/desired_eeg_srate))
time_eeg = np.arange(0,np.max(time_eeg),1/desired_eeg_srate)

spectrotime, frequencies, spectraldata = spectrogram(hip_eeg,desired_eeg_srate, windowlength=1, overlap=0.5, highfreq=1000)
#spectraldata = np.log(spectraldata)

new_time = np.zeros(int(np.ceil(len(spectrotime) / 2)))
unit_data = phy_data(phy_path, unit_srate)

seizure_start = 3782.72
unit_start = seizure_start - 120
unit_stop = seizure_start + 120

unit_bins = np.zeros((len(new_time), len(spike_templates)))
for i in range(len(spike_templates)):
    spikes = unit_data.extract_spike_trains(spike_templates[i])
    unit_bins[:,i] = extract_frequency_bins(spikes, 1, unit_start, unit_stop)[1]
    unit_bins[:,i] = movingaverage(unit_bins[:,i], 5)


original_spectraldata = spectraldata

start = int(130*2)
stop = int(180*2-1)
spectraldata_postictal = spectraldata[:, start:stop]
spectrotime_postictal = spectrotime[start:stop]
new_time_postictal = np.zeros(int(np.ceil(len(spectrotime_postictal) / 2)))
new_spectraldata_postictal = np.zeros((np.shape(spectraldata_postictal)[0], int(np.ceil(len(spectrotime_postictal) / 2))))

for i in range(len(new_time_postictal)):
    if i == 0:
        new_time_postictal[i] = spectrotime_postictal[i]
        new_spectraldata_postictal[:,i] = spectraldata_postictal[:,i]
    elif i == len(new_time_postictal) - 1:
        print('yes')
        new_time_postictal[i] = spectrotime_postictal[-1]
        new_spectraldata_postictal[:,i] = spectraldata_postictal[:,-1]
    else:
        start = int(i*2 - 1)
        stop = int(start+3)
        new_time_postictal[i] = np.mean(spectrotime_postictal[start:stop])
        new_spectraldata_postictal[:,i] = np.mean(spectraldata_postictal[:, start:stop], axis=1)

post_ictal_power = np.sqrt(np.sum(new_spectraldata_postictal, axis=0))        
new_spectraldata_postictal = np.sqrt(new_spectraldata_postictal)        

start = int(60*2)
stop = int(120*2-1)
spectraldata_preictal = spectraldata[:, start:stop]
spectrotime_preictal = spectrotime[start:stop]
new_time_preictal = np.zeros(int(np.ceil(len(spectrotime_preictal) / 2)))
new_spectraldata_preictal = np.zeros((np.shape(spectraldata_preictal)[0], int(np.ceil(len(spectrotime_preictal) / 2))))

for i in range(len(new_time_preictal)):
    if i == 0:
        new_time_preictal[i] = spectrotime_preictal[i]
        new_spectraldata_preictal[:,i] = spectraldata_preictal[:,i]
    elif i == len(new_time_preictal) - 1:
        print('yes')
        new_time_preictal[i] = spectrotime_preictal[-1]
        new_spectraldata_preictal[:,i] = spectraldata_preictal[:,-1]
    else:
        start = int(i*2 - 1)
        stop = int(start+3)
        new_time_preictal[i] = np.mean(spectrotime_preictal[start:stop])
        new_spectraldata_preictal[:,i] = np.mean(spectraldata_preictal[:, start:stop], axis=1)

pre_ictal_power = np.sqrt(np.sum(new_spectraldata_preictal, axis=0))

original_spectraldata = np.sqrt(original_spectraldata)

frequency_bins_in = unit_bins[:,0]
frequency_bins_ex = unit_bins[:,1]
frequency_bins_uncoupled = unit_bins[:,2]
time_bins = np.arange(0.5, 60*4, 1)

#plt.scatter(post_ictal_power, frequency_bins_ex[130:180])


#plt.scatter(post_ictal_power, frequency_bins_ex[130:180])
new_time = np.zeros(int(np.ceil(len(spectrotime) / 2)))
new_spectraldata = np.zeros((np.shape(spectraldata)[0], int(np.ceil(len(spectrotime) / 2))))

for i in range(len(new_time)):
    if i == 0:
        new_time[i] = spectrotime[i]
        new_spectraldata[:,i] = spectraldata[:,i]
        
    elif i == len(new_time) - 1:
        start = int(i*2 - 1)
        new_time[i] = np.mean(spectrotime[start:])
        new_spectraldata[:,i] = np.mean(spectraldata[:, start:], axis=1)
    else:
        start = int(i*2 - 1)
        stop = int(start+3)
        new_time[i] = np.mean(spectrotime[start:stop])
        new_spectraldata[:,i] = np.mean(spectraldata[:, start:stop], axis=1)

total_amplitude = np.squeeze(np.sum(new_spectraldata, axis=0))
#total_amplitude = movingaverage(total_amplitude, 5)
new_time = new_time + 1

total_amplitude[119:129] = 0



colors = np.arange(0, 500, 1)
### PLOTTING ###
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(5, figsize=(6,6))

ax[0].plot(time_eeg-120,hip_eeg,'k', label='Baseline (pre-stimulation)', linewidth=0.3)
ax[0].vlines(0, -15, 15, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[0].vlines(10, -15, 15, colors='r', linewidth=3, linestyles='dashed', zorder=3)
#ax[1].contourf(spectrotime-120, frequencies, original_spectraldata, levels = colors, cmap = 'jet', extend = 'both')
ax[1].bar(new_time-120,total_amplitude, color = 'k', width = 0.8)
ax[1].vlines(0, 0, np.max(total_amplitude), colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[1].vlines(10, 0, np.max(total_amplitude), colors='r', linewidth=3, linestyles='dashed', zorder=3)

#ax[1].vlines(0, 0, 1000, colors='r', linewidth=3, linestyles='dashed', zorder=3)
#ax[1].vlines(10, 0, 1000, colors='r', linewidth=3, linestyles='dashed', zorder=3)

ax[3].bar(time_bins-120,frequency_bins_ex, color = 'r', width = 0.8)
ax[3].vlines(0, 0, 20, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[3].vlines(10, 0, 20, colors='r', linewidth=3, linestyles='dashed', zorder=3)

ax[2].bar(time_bins-120,frequency_bins_in, color = 'b', width = 0.8)
ax[2].vlines(0, 0, 20, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[2].vlines(10, 0, 20, colors='r', linewidth=3, linestyles='dashed', zorder=3)

ax[4].bar(time_bins-120,frequency_bins_uncoupled, color = 'g', width = 0.8)
ax[4].vlines(0, 0, 20, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[4].vlines(10, 0, 20, colors='r', linewidth=3, linestyles='dashed', zorder=3)

ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_linewidth(2)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_linewidth(2)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[2].spines['left'].set_linewidth(2)
ax[3].spines['top'].set_visible(False)
ax[3].spines['right'].set_visible(False)
ax[3].spines['bottom'].set_visible(False)
ax[3].spines['left'].set_linewidth(2)
ax[4].spines['top'].set_visible(False)
ax[4].spines['right'].set_visible(False)
ax[4].spines['bottom'].set_linewidth(2)
ax[4].spines['left'].set_linewidth(2)
#ax[1].set_yscale('log')
ax[0].set_ylim([-15,15])
#ax[1].set_ylim([1,1000])
ax[2].set_ylim([0,5])
ax[3].set_ylim([0,30])
ax[4].set_ylim([0,5])
ax[0].set_xlim([-20,60])
ax[1].set_xlim([-20,60])
ax[2].set_xlim([-20,60])
ax[3].set_xlim([-20,60])
ax[4].set_xlim([-20,60])
ax[0].set_xticklabels([],[])
ax[0].tick_params(axis='x', which='both', bottom=False)
ax[1].set_xticklabels([],[])
ax[1].tick_params(axis='x', which='both', bottom=False)
ax[2].set_xticklabels([],[])
ax[2].tick_params(axis='x', which='both', bottom=False)
ax[3].set_xticklabels([],[])
ax[3].tick_params(axis='x', which='both', bottom=False)
ax[4].xaxis.set_tick_params(width=2)
ax[3].tick_params(axis="x", labelsize=14)
ax[0].set_yticks([-10, 10])
ax[0].set_yticklabels([-10, 10])
ax[1].set_yticks([0, 5])
ax[1].set_yticklabels([0, 5])
#ax[1].tick_params(axis='y', which='both', left=False)
ax[2].set_yticks([0, 5])
ax[2].set_yticklabels([0, 5])
ax[3].set_yticks([0, 20])
ax[3].set_yticklabels([0, 20])
ax[4].set_yticks([0, 5])
ax[4].set_yticklabels([0, 5])
ax[0].set_ylabel('mV')
ax[1].set_ylabel('a.u.')
ax[2].set_ylabel('Hz')
ax[3].set_ylabel('Hz')
ax[4].set_ylabel('Hz')
ax[4].set_xlabel('Time (s)', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.abspath(savepath), dpi=600)
#plt.close()


### CROSS CORRELATIONS

### LL001
eeg_path = r'E:\Manuscript_analysis_files\data\lfp_data\LL001\Seizure1.mat'
eeg_totalpowerate = 30000
eeg_stimstart = 120
phy_path = r'E:\Manuscript_analysis_files\data\phy_data\LL001-final.GUI'
LL001_spike_templates = [5, 6, 10, 48, 59, 91, 92, 101, 110, 114, 118, 125, 149, 164, 165, 176, 178, 180, 197]
seizure_start = 3782.72
eeg_downsample_rate = 2000
unit_totalpowerate = 30000

LL001_all_freqs, LL001_totalpower = eegpower_spikerate_crosscorr(eeg_path, eeg_totalpowerate,
                             eeg_downsample_rate,
                             eeg_stimstart,
                             phy_path,
                             unit_totalpowerate,
                             LL001_spike_templates,
                             seizure_start
                             )

a = LL001_totalpower[:,16]
b = LL001_totalpower[:,13]
c = LL001_totalpower[:,14] # 165

lags = correlation_lags(50, 50, mode='same')


a_max = int(np.where(np.max(np.abs(a)) == np.abs(a))[0])
b_max = int(np.where(np.max(np.abs(b)) == np.abs(b))[0])
c_max = int(np.where(np.max(np.abs(c)) == np.abs(c))[0])





#plt.plot(corrs_in)

#m_frequency_bins_ex, b_frequency_bins_ex = np.polyfit(post_ictal_power, frequency_bins_ex[130:], 1)
#corr_ex = pearsonr(post_ictal_power, frequency_bins_ex[130:])[0] ** 2

### correlations
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(3, figsize=(4,6))

ax[0].plot(lags, a, c = 'b', linewidth = 2)
ax[0].vlines(0,-1,1,color='k', linestyles='--', linewidth=1.5)
#ax[0].set_xscale('log')
#ax[0].text(50, 12, '$R^2$ = ' + str(np.around(corr_ex, 2)), fontsize=14)
#ax[0].plot(post_ictal_power, m_frequency_bins_ex*post_ictal_power + b_frequency_bins_ex)
ax[1].plot(lags, b, c = 'r', linewidth = 2)
ax[1].vlines(0,-1,1,color='k', linestyles='--', linewidth=1.5)
#ax[1].set_xscale('log')

ax[2].plot(lags, c, c = 'g', linewidth = 2)
ax[2].vlines(0,-1,1,color='k', linestyles='--', linewidth=1.5)

ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_linewidth(2)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_linewidth(2)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_linewidth(2)
ax[2].spines['left'].set_linewidth(2)

ax[0].set_xticklabels([],[])
ax[0].tick_params(axis='x', which='both', bottom=False)
ax[1].set_xticklabels([],[])
ax[1].tick_params(axis='x', which='both', bottom=False)
#ax[1].set_xticks([1, 10, 100, 1000])
#ax[1].set_xticklabels([1, 10, 100, 1000], fontsize = 12)

ax[0].set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
ax[0].set_yticklabels([-1.0, -0.5, 0, 0.5, 1.0], fontsize = 12)
ax[1].set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
ax[1].set_yticklabels([-1.0, -0.5, 0, 0.5, 1.0], fontsize = 12)
ax[2].set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
ax[2].set_yticklabels([-1.0, -0.5, 0, 0.5, 1.0], fontsize = 12)


ax[0].set_ylabel('xcorr',fontsize=16, fontweight='bold')
ax[2].set_xlabel('Lag (s)',fontsize=16, fontweight='bold')
ax[1].set_ylabel('xcorr',fontsize=16, fontweight='bold')
ax[2].set_ylabel('xcorr',fontsize=16, fontweight='bold')

ax[0].text(10, 0.6, 'max R = ' + str(np.around(a[a_max], decimals = 2)) + '\nlag = ' + str(lags[a_max]) + ' s', fontsize = 10)
ax[1].text(10, 0.6, 'max R = ' + str(np.around(b[b_max], decimals = 2)) + '\nlag = ' + str(lags[b_max]) + ' s', fontsize = 10)
ax[2].text(10, -0.8, 'max R = ' + str(np.around(c[c_max], decimals = 2)) + '\nlag = ' + str(lags[c_max]) + ' s', fontsize = 10)

ax[0].set_ylim(-1,1)
ax[1].set_ylim(-1,1)
ax[2].set_ylim(-1,1)
ax[0].set_xlim(-25,25)
ax[1].set_xlim(-25,25)
ax[2].set_xlim(-25,25)
plt.tight_layout()
plt.savefig(os.path.abspath(savepath2), dpi=600)





