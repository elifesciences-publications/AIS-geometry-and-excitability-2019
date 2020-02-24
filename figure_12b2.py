

"""
AIS geometry and excitability, supplementary figure.

Extended AIS.

"""
from __future__ import print_function
from brian2 import *
from shared import params_model_description, model_Na_Kv1, measure_current_threshold, measure_voltage_threshold
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
    data = load('figure_14b2.npz')
    starts = data['arr_0']*um
    length = data['arr_1']*um
    thresholds = data['arr_2']*mV
    thresholds_long = data['arr_3']*mV


else: # running simulations
    
    # A function to measure of the threshold BIO model in CC, for different AIS start positions and Gna
    def BIO_model_threshold_in_CC(ais_start, ais_end, gna_tot, morpho):
        print ('AIS start:', ais_start, 'Gna:', gna_tot)
        params = params_model_description
        pulse_length = 50.*ms
        
        neuron = model_Na_Kv1(params=params, resting_vm =-75.*mV, Na_start = ais_start, Na_end = ais_end, gna_tot=gna_tot, density=False, morpho=morpho)
        i_rheo = measure_current_threshold(params, neuron, -75.*mV, ais_start, ais_end, pulse_length)
        print ('Rheobase:', i_rheo)
        
        neuron = model_Na_Kv1(params=params, resting_vm =-75.*mV, Na_start = ais_start, Na_end = ais_end, gna_tot=gna_tot, density=False, morpho=morpho)
        v_thres, _,_,_ = measure_voltage_threshold(params, neuron, -75.*mV, ais_start, ais_end, i_rheo, pulse_length) 
        print ('Voltage threshold:', v_thres)
        return v_thres 
    
    print ('BIO model')
    
    dend_diam = 6.*um 
    dend_length = 1000.*um 
    axon_diam = 1. * um
    axon_length = 500. * um 
    soma_diameter = 30. * um
    morpho = Soma(diameter=soma_diameter)
    dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=500)
    axon = Cylinder(diameter=axon_diam, length=axon_length, n=500)
    morpho.dendrite = dendrite
    morpho.axon = axon
    
    # run the simulations with joblib
    if __name__ == '__main__':
        results = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC)(start, start+length, GNa, morpho) for start in starts)
        thresholds = array(results)*1e3
    
    print ('BIO model with longer axon')
    
    # Changing axon length
    dend_diam1 = 6.*um 
    dend_length1 = 1000.*um 
    axon_diam1 = 1. * um
    axon_length1 = 2000. * um 
    soma_diameter1 = 30. * um
    morpho1 = Soma(diameter=soma_diameter1)
    dendrite1 = Cylinder(diameter=dend_diam1, length=dend_length1, n=500)
    axon1 = Cylinder(diameter=axon_diam1, length=axon_length1, n=2000)
    morpho1.dendrite = dendrite1
    morpho1.axon = axon1
    
    # run the simulations with joblib
    if __name__ == '__main__':
        results_long = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC)(start, start+length, GNa, morpho1) for start in starts)
        thresholds_long = array(results_long)*1e3
        
    # Save the data in an npz file
    savez('figure_14b2', starts/um, length/um, thresholds/mV, thresholds_long/mV)
    
x_mids = starts+length/2 # AIS middle position

# Slopes
slopes_long, _, _, _, _ = stats.linregress(-log(ra*x_mids* GNa), thresholds_long)
slopes_500, _, _, _, _ = stats.linregress(-log(ra*x_mids* GNa), thresholds)

print ('Slopes:')
print ('Biophysical model in current-clamp:')
print('Threshold vs delta:', slopes_500, 'mV')
print('Threshold vs delta long axon:', slopes_long, 'mV')

### Plots

fig = figure(2, figsize=(4,3))

ax3 = subplot(111)

ax3.semilogx(x_mids/um, thresholds, '--', color='k', label='original axon: k=%.2f mV' %slopes_500)
ax3.semilogx(x_mids/um, thresholds_long, '-', color='forestgreen', label='2 mm  long axon: k=%.2f mV' %slopes_long)
ax3.set_xlabel('$x_{1/2}$ ($\mu$m)') 
ax3.set_ylabel('$V_s$ (mV)') 
ax3.set_ylim(-65,-50)
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax3.xaxis.set_minor_formatter(ScalarFormatter())
ax3.legend(frameon=False, fontsize=8)

tight_layout()

show()    






