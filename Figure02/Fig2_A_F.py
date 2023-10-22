
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import os
import neo

import sys
sys.path.insert(1, r'E:\Manuscript_analysis_files')

from vies.lfp.filter import butter_bandpass_filter 
from vies.parse.neuron import load_neuronfile
from vies.lfp.lfp_analysis import spectrogram
from vies.parse.spike2 import load_spikematfile, load_spikematfile_bins
from vies.spike.autocorrelogram import compute_autocorrelogram, spike_train_correlogram

## Params
lowcut_eeg=2
highcut_eeg=500
srate_eeg=30000
desired_eeg_srate=2000
hip_lfp_data = r'E:\Manuscript_analysis_files\data\lfp_data\Exp1\Seizure2.mat'
lc_data_excitation = r'E:\Manuscript_analysis_files\data\spike_data\Exp1\Neuron1_Seizure2.mat'
lc_data_inhibition = r'E:\Manuscript_analysis_files\data\spike_data\Exp1\Neuron2_Seizure2.mat'
savepath1 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\output\ExampleResponses_part1.png'
savepath2 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\output\ExampleResponses_part2.png'
savepath3 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\output\waveform_inh.png'
savepath4 = r'E:\Manuscript_analysis_files\LC_seizure_python_scripts\Figure02\output\waveform_ex.png'

## Load and preprocess data
time_eeg, hip_eeg = load_neuronfile(hip_lfp_data,srate=30000,channel=3, gain=200, inputrange=20)
hip_eeg=butter_bandpass_filter(hip_eeg, lowcut_eeg, highcut_eeg, srate_eeg, order=2)
hip_eeg = signal.decimate(hip_eeg, int(srate_eeg/desired_eeg_srate))
time_eeg = np.arange(0,np.max(time_eeg),1/desired_eeg_srate)

spectrotime, frequencies, spectraldata = spectrogram(hip_eeg,desired_eeg_srate, windowlength=1, overlap=0.9, highfreq=100)
spectraldata = np.log(spectraldata)

LC_spiketimes_ex, LC_waveformtimes_ex, LC_waveforms_ex = load_spikematfile(lc_data_excitation, 1, 30000, 1000)
time_bins_ex, frequency_bins_ex = load_spikematfile_bins(lc_data_excitation, 1, 0, 180, binsize=5)

LC_spiketimes_in, LC_waveformtimes_in, LC_waveforms_in = load_spikematfile(lc_data_inhibition, 1, 30000, 1000)
time_bins_in, frequency_bins_in = load_spikematfile_bins(lc_data_inhibition, 1, 0, 180, binsize=5)

seizurestart = 65.4+60
seizurestop = 98.2+60



### PLOTTING ### A C E
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(3, figsize=(9,7))
ax[0].plot(time_eeg-120,hip_eeg,'b', label='Baseline (pre-stimulation)', linewidth=0.3)
ax[0].vlines(0, -15, 15, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[0].vlines(10, -15, 15, colors='r', linewidth=3, linestyles='dashed', zorder=3)

for i in range(LC_waveforms_ex.shape[1]):
    ax[1].plot(LC_waveformtimes_ex[:,i]-120,LC_waveforms_ex[:,i],'g')   
ax[1].vlines(0, -0.2, 0.2, colors='r', linewidth=3, linestyles='dashed', zorder=3)

for i in range(LC_waveforms_in.shape[1]):
    ax[2].plot(LC_waveformtimes_in[:,i]-120,LC_waveforms_in[:,i],'m')
ax[2].vlines(0, -0.2, 0.2, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[2].vlines(10, -0.2, 0.2, colors='r', linewidth=3, linestyles='dashed', zorder=3)

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

ax[0].set_ylim([-10,10])

ax[1].set_ylim([-0.2,0.2])
ax[2].set_ylim([-0.2,0.2])
ax[0].set_xlim([-30,50])
ax[1].set_xlim([-30,50])
ax[2].set_xlim([-30,50])
ax[0].set_xticklabels([],[])
ax[0].tick_params(axis='x', which='both', bottom=False)
ax[0].tick_params(axis='x', which='both', bottom=False)
ax[1].set_xticklabels([],[])
ax[1].tick_params(axis='x', which='both', bottom=False)
ax[0].yaxis.set_tick_params(width=2)
ax[1].yaxis.set_tick_params(width=2)
ax[2].yaxis.set_tick_params(width=2)
ax[0].tick_params(axis="y", labelsize=14)
ax[1].tick_params(axis="y", labelsize=14)
ax[2].tick_params(axis="y", labelsize=14)
ax[2].tick_params(axis="x", labelsize=14)
ax[0].set_ylabel('Hip LFP\nAmplitude (mV)', fontsize=16, fontweight='bold')
ax[1].set_ylabel('LC Neuron 1\nAmplitude (mV)', fontsize=16, fontweight='bold')
ax[2].set_ylabel('LC Neuron 2\nAmplitude (mV)', fontsize=16, fontweight='bold')
ax[2].set_xlabel('Time (s)', fontsize=16, fontweight='bold')
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.tight_layout()
plt.savefig(os.path.abspath(savepath1), dpi=600)
#plt.close()

#plt.show()

### B D F
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(3, figsize=(9,7))

ax[0].contourf(spectrotime-120,frequencies,spectraldata,300,cmap='jet')
ax[0].vlines(0, 0, np.max(frequencies), colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[0].vlines(10, 0, np.max(frequencies), colors='r', linewidth=3, linestyles='dashed', zorder=3)

ax[1].bar(time_bins_ex-120,frequency_bins_ex, color = 'g', width = 4.5)
ax[1].vlines(0, 0, 20, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[1].vlines(10, 0, 20, colors='r', linewidth=3, linestyles='dashed', zorder=3)

ax[2].bar(time_bins_in-120,frequency_bins_in, color = 'm', width = 4.5)
ax[2].vlines(0, 0, 20, colors='r', linewidth=3, linestyles='dashed', zorder=3)
ax[2].vlines(10, 0, 20, colors='r', linewidth=3, linestyles='dashed', zorder=3)

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
ax[2].spines['left'].set_linewidth(2)
ax[2].spines['bottom'].set_linewidth(2)

ax[0].set_ylim([0,np.ceil(np.max(frequencies))])
ax[1].set_ylim([0,17])
ax[2].set_ylim([0,3])

ax[0].set_xlim([-30,50])
ax[1].set_xlim([-30,50])
ax[2].set_xlim([-30,50])

ax[0].set_xticklabels([],[])
ax[1].set_xticklabels([],[])
ax[1].tick_params(axis='x', which='both', bottom=False)
ax[2].xaxis.set_tick_params(width=2)
ax[0].yaxis.set_tick_params(width=2)

ax[1].yaxis.set_tick_params(width=2)
ax[2].yaxis.set_tick_params(width=2)

ax[0].tick_params(axis="y", labelsize=14)
ax[1].tick_params(axis="y", labelsize=14)
ax[2].tick_params(axis="y", labelsize=14)
ax[2].tick_params(axis="x", labelsize=14)

ax[0].set_ylabel('Hip LFP\nFrequency (Hz)', fontsize=16, fontweight='bold')
ax[1].set_ylabel('LC Neuron 1\nFrequency (Hz)', fontsize=16, fontweight='bold')
ax[2].set_ylabel('LC Neuron 2\nFrequency (Hz)', fontsize=16, fontweight='bold')
ax[2].set_xlabel('Time (s)', fontsize=16, fontweight='bold')

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.tight_layout()
plt.savefig(os.path.abspath(savepath2), dpi=600)
#plt.close()


waveform_in = np.mean(LC_waveforms_in, axis=1)
waveform_ex = np.mean(LC_waveforms_ex, axis=1)


srate = 30000
dt = 1/srate
t_waveform = (np.arange(0,126) / srate) * 1000 
t_waveform = t_waveform - np.mean(t_waveform)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

plot = ax.plot(t_waveform, waveform_in, linewidth = 8, color='m')    # first image on screen
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlim([np.min(t_waveform)/2, np.max(t_waveform)/2])
ax.set_ylim([np.min(waveform_in), np.max(waveform_in)])
ax.set_xlabel('ms', fontsize=16*3.5)
ax.set_yticks([]) 
ax.set_yticklabels([], fontsize=10, fontweight='bold') 
ax.set_xticks([-1,0,1])
ax.set_xticklabels([-1, 0, 1], fontsize=14*3.5, fontweight='bold') 

plt.tight_layout()
plt.savefig(savepath3, dpi=300)
#plt.close()


## inh insert
srate = 30000
dt = 1/srate
t_waveform = (np.arange(0,126) / srate) * 1000 
t_waveform = t_waveform - np.mean(t_waveform)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
#plt.xticks(fontsize=14)

plot = ax.plot(t_waveform, waveform_ex, linewidth = 8, color='g')    # first image on screen
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlim([np.min(t_waveform)/2, np.max(t_waveform)/2])
ax.set_ylim([np.min(waveform_ex) + 0.1*np.min(waveform_ex), np.max(waveform_ex) + 0.1*np.max(waveform_ex)])
ax.set_xlabel('ms', fontsize=16*3.5)
ax.set_yticks([]) 
ax.set_yticklabels([], fontsize=10, fontweight='bold') 
ax.set_xticks([-1,0,1])
ax.set_xticklabels([-1, 0, 1], fontsize=14*3.5, fontweight='bold') 

plt.tight_layout()
plt.savefig(savepath4, dpi=300)


