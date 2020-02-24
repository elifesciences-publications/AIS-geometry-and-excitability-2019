

"""
AIS geometry and excitability, supplementary figure.

Extended AIS with the axon originating from the soma.

"""
from __future__ import print_function
from brian2 import *
from shared import params_model_description, model_Na_Kv1, model_Na_Kv1_axodendritic, measure_current_threshold, measure_voltage_threshold
from joblib import Parallel, delayed
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.ticker as ticker
from scipy import stats

only_plotting = True # to plot the figure without running the simulations

# Parameters
defaultclock.dt = 0.005*ms

n=7
params = params_model_description#params_all
length = 30.*um # AIS length
GNa = 400.*nS # total Na conductance at the AIS
starts = linspace(0., 30., n)*um # AIS start positions

### SIMULATIONS in the biophysical model

if only_plotting: # loading data
    data = load('figure_14d.npz')
    starts = data['arr_0']*um
    length = data['arr_1']*um
    Ra = data['arr_2']*Mohm
    Ra_axo_dend = data['arr_3']*Mohm
    thresholds = data['arr_4']*mV
    thresholds_axo_dend = data['arr_5']*mV
    
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
    def BIO_model_threshold_in_CC_axodendritic(ais_start, ais_end, gna_tot):
        print ('AIS start:', ais_start, 'Gna:', gna_tot)
        params = params_model_description
        pulse_length = 50.*ms
        
        neuron = model_Na_Kv1_axodendritic(params=params, resting_vm =-75.*mV, Na_start = ais_start, Na_end = ais_end, gna_tot=gna_tot, density=False)
        i_rheo = measure_current_threshold(params, neuron, -75.*mV, ais_start, ais_end, pulse_length, axo_dend=True)
        print ('Rheobase:', i_rheo)
        
        neuron = model_Na_Kv1_axodendritic(params=params, resting_vm =-75.*mV, Na_start = ais_start, Na_end = ais_end, gna_tot=gna_tot, density=False)
        v_thres, _,_,_ = measure_voltage_threshold(params, neuron, -75.*mV, ais_start, ais_end, i_rheo, pulse_length, axo_dend=True) 
        print ('Voltage threshold:', v_thres)
        return v_thres
    
    print ('BIO model axo-dendritic neuron')
    
    # run the simulations with joblib
    if __name__ == '__main__':
        results_axo_dend = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC_axodendritic)(start, start+length,  GNa) for start in starts)
        thresholds_axo_dend = array(results_axo_dend)*1e3

# Axial resistance
# aaxo-somatic
x_mids = starts+length/2 # AIS middle position
Ra = 4*params.Ri*x_mids/(pi*params.axon_diam**2)

# axo-dendritic
stem_diam = 2.*um
stem_length = 7.*um
branch_diam = stem_diam/(2**(2./3))
Ra_axo_dend = 4*params.Ri*stem_length/(pi*stem_diam**2) + 4*params.Ri*x_mids/(pi*branch_diam**2)

# Save the data in an npz file
savez('figure_14d', starts/um, length/um, Ra/Mohm, Ra_axo_dend/Mohm, thresholds/mV, thresholds_axo_dend/mV)

# Slopes
slopes_axo_dend, _, _, _, _ = stats.linregress(-log(Ra_axo_dend/Mohm), thresholds_axo_dend)
slopes_500, _, _, _, _ = stats.linregress(-log(Ra/Mohm), thresholds)

print ('Slopes:')
print ('Biophysical model in current-clamp:')
print('Threshold vs delta:', slopes_500, 'mV')
print('Threshold vs delta axo_dend:', slopes_axo_dend, 'mV')

### Plots

fig = figure(2, figsize=(4,3))

ax3 = subplot(111)

ax3.semilogx(Ra/Mohm, thresholds, '--', color='k', label='original axon: k=%.2f mV' %slopes_500)
ax3.semilogx(Ra_axo_dend/Mohm, thresholds_axo_dend, '-', color='forestgreen', label='axo-dendritic, k=%.2f mV' %slopes_axo_dend)
ax3.set_xlabel('$R_a$ (M$\Omega$)') 
ax3.set_ylabel('$V_s$ (mV)') 
ax3.set_ylim(-65,-45)
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax3.xaxis.set_minor_formatter(ScalarFormatter())
ax3.legend(frameon=False, fontsize=8)

tight_layout()

show()   
 






