

"""
AIS geometry and excitability, supplementary figure.

Extended AIS in a larger neuron.

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
GNa = 400.*nS # total Na conductance at the AIS


### SIMULATIONS in the biophysical model

if only_plotting: # loading data
    data = load('figure_14a.npz')
    starts = data['arr_0']*um
    starts_large = data['arr_1']*um
    length = data['arr_2']*um
    scaling = data['arr_3']
    thresholds = data['arr_4']
    thresholds_large = data['arr_5']
    
    length_large = sqrt(scaling) * length

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
    starts = linspace(0., 30., n)*um # AIS start positions
    length = 30.*um # AIS length
    if __name__ == '__main__':
        results = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC)(start, start+length, GNa, morpho) for start in starts)
        thresholds = array(results)
    
    print ('BIO model in a larger neuron')
    
    scaling = 3.
    axon_diam_large = scaling * axon_diam
    soma_diameter_large = (axon_diam_large/metre)**(3./4) * 1e6 * um
    dend_diam_large = dend_diam
    dend_length_large = 1000.*um 
    axon_length_large = 500. * um 
    morpho_large = Soma(diameter=soma_diameter_large)
    dendrite_large = Cylinder(diameter=dend_diam_large, length=dend_length_large, n=500)
    axon_large = Cylinder(diameter=axon_diam_large, length=axon_length_large, n=500)
    morpho_large.dendrite = dendrite_large
    morpho_large.axon = axon_large
    
    # run the simulations with joblib
    starts_large = sqrt(scaling) * linspace(0., 30., n)*um # AIS start positions scale with sqrt(axon_diam) 
    length_large = sqrt(scaling) * length # AIS length
    
    if __name__ == '__main__':
        results_large = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC)(start, start+length_large, scaling**(3./2)*GNa, morpho_large) for start in starts_large)
        thresholds_large = array(results_large)
        
    # # Save the data in an npz file
    savez('figure_14a', starts/um, starts_large/um, length/um, scaling, thresholds, thresholds_large)
    
x_mids = starts+length/2 # AIS middle position
x_mids_large = starts_large+length_large/2 # AIS middle position

# Slopes
slopes_large, _, _, _, _ = stats.linregress(-log(x_mids_large/um), thresholds_large)
slopes_500, _, _, _, _ = stats.linregress(-log(x_mids/um), thresholds)

print ('Slopes:')
print ('Biophysical model in current-clamp:')
print('Threshold vs delta:', slopes_500*1e3, 'mV')
print('Threshold vs delta large neuron:', slopes_large*1e3, 'mV')

### Plots

fig = figure(2, figsize=(4,3))

ax3 = subplot(111)

ax3.semilogx(x_mids/um, thresholds*1e3, '--', color='k', label='original axon: k=%.2f mV' %slopes_500)
ax3.semilogx(x_mids/um, thresholds*1e3 + 5./2*log(scaling), '--', color='k', label='original axon, predicted shift', alpha=0.5)
ax3.semilogx(x_mids_large/um, thresholds_large*1e3, '-', color='forestgreen', label='large neuron, k=%.2f mV' %slopes_large)
ax3.set_xlabel('$x_{1/2}$ ($\mu$m)') 
ax3.set_ylabel('$V_s$ (mV)') 
ax3.set_ylim(-65,-40)
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax3.xaxis.set_minor_formatter(ScalarFormatter())
ax3.legend(frameon=False, fontsize=8)

tight_layout()

show()   
 






