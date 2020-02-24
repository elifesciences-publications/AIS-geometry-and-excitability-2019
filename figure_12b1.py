

"""
AIS geometry and excitability, supplementary figure.

Extended AIS with myelin that starts beyond it.

"""
from __future__ import print_function
from brian2 import *
from shared import params_model_description, model_Na_Kv1, model_Na_Kv1_myelin, measure_current_threshold, measure_voltage_threshold
from joblib import Parallel, delayed
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.ticker as ticker
from scipy import stats

only_plotting = True # to plot the figure without running the simulations

# Parameters
defaultclock.dt = 0.005*ms

n=7
params = params_model_description
ra = 4.*params.Ri/(pi*params.axon_diam**2) # axial resistance per unit length
length = 30.*um # AIS length
GNa = 400.*nS # total Na conductance at the AIS
starts = linspace(0., 30., n)*um # AIS start positions

### SIMULATIONS in the biophysical model

if only_plotting: # loading data
    data = load('figure_14b1.npz')
    starts = data['arr_0']*um
    length = data['arr_1']*um
    thresholds = data['arr_2']*mV
    thresholds_myelin = data['arr_3']*mV


else: # running simulations
    
    # A function to measure of the threshold BIO model in CC, for different AIS start positions and Gna
    def BIO_model_threshold_in_CC(ais_start, ais_end, gna_tot):
        print ('AIS start:', ais_start, 'Gna:', gna_tot)
        params = params_model_description
        pulse_length = 50.*ms
        
        neuron = model_Na_Kv1(params=params, resting_vm =-75.*mV, Na_start = ais_start, Na_end = ais_end, gna_tot=gna_tot, density=False)
        i_rheo = measure_current_threshold(params, neuron, -75.*mV, ais_start, ais_end, pulse_length)
        print ('Rheobase:', i_rheo)
        
        neuron = model_Na_Kv1(params=params, resting_vm =-75.*mV, Na_start = ais_start, Na_end = ais_end, gna_tot=gna_tot, density=False)
        v_thres, _,_,_ = measure_voltage_threshold(params, neuron, -75.*mV, ais_start, ais_end, i_rheo, pulse_length) 
        print ('Voltage threshold:', v_thres)
        return v_thres 
    
    print ('BIO model')
        
    # run the simulations with joblib
    if __name__ == '__main__':
        results = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC)(start, start+length, GNa) for start in starts)
        thresholds = array(results)*1e3
    
    
    # A function to measure of the threshold BIO model in CC, for different AIS start positions and Gna
    def BIO_model_threshold_in_CC_myelin(ais_start, ais_end, gna_tot):
        print ('AIS start:', ais_start, 'Gna:', gna_tot)
        params = params_model_description
        pulse_length = 50.*ms
        
        neuron = model_Na_Kv1_myelin(params=params, resting_vm =-75.*mV, Na_start = ais_start, Na_end = ais_end, gna_tot=gna_tot, density=False)
        i_rheo = measure_current_threshold(params, neuron, -75.*mV, ais_start, ais_end, pulse_length)
        print ('Rheobase:', i_rheo)
        
        neuron = model_Na_Kv1_myelin(params=params, resting_vm =-75.*mV, Na_start = ais_start, Na_end = ais_end, gna_tot=gna_tot, density=False)
        v_thres, _,_,_ = measure_voltage_threshold(params, neuron, -75.*mV, ais_start, ais_end, i_rheo, pulse_length) 
        print ('Voltage threshold:', v_thres)
        return v_thres 
    
    print ('BIO model with myelin')
    
    # run the simulations with joblib
    if __name__ == '__main__':
        results_myelin = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC_myelin)(start, start+length, GNa) for start in starts)
        thresholds_myelin = array(results_myelin)*1e3
        
    # Save the data in an npz file
    savez('figure_14b1', starts/um, length/um, thresholds/mV, thresholds_myelin/mV)
    
x_mids = starts+length/2 # AIS middle position

# Slopes
slopes_myelin, _, _, _, _ = stats.linregress(-log(x_mids/um), thresholds_myelin)
slopes, _, _, _, _ = stats.linregress(-log(x_mids/um), thresholds)

print ('Slopes:')
print ('Biophysical model in current-clamp:')
print('Threshold vs delta:', slopes, 'mV')
print('Threshold vs delta myelin:', slopes_myelin, 'mV')

### Plots

fig = figure(2, figsize=(4,3))

ax3 = subplot(111)

ax3.semilogx(x_mids/um, thresholds, '--', color='k', label='original axon: k=%.2f mV' %slopes)
ax3.semilogx(x_mids/um, thresholds_myelin, '-', color='forestgreen', label='myelin: k=%.2f mV' %slopes_myelin)
ax3.set_xlabel('$x_{1/2}$ ($\mu$m)') 
ax3.set_ylabel('$V_s$ (mV)') 
ax3.set_ylim(-70,-50)
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax3.xaxis.set_minor_formatter(ScalarFormatter())
ax3.legend(frameon=False, fontsize=8)

tight_layout()

show()    






