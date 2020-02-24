

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
params = params_model_description#params_all
ra = 4.*params.Ri/(pi*params.axon_diam**2) # axial resistance per unit length
length = 30.*um # AIS length
GNa = 400.*nS # total Na conductance at the AIS
starts = linspace(0., 30., n)*um # AIS start positions

### SIMULATIONS in the biophysical model

if only_plotting: # loading data
    data = load('figure_14c.npz')
    starts = data['arr_0']*um
    length =  data['arr_1']*um
    thresholds = data['arr_2']*mV
    thresholds_branched_rall = data['arr_3']*mV
    thresholds_branched_emp = data['arr_4']*mV
    
    branch_diam_rall = data['arr_6']*um #1.*um/(2**(2./3))
    branch_diam_emp = data['arr_5']*um #1.*um/(2**(5./2))

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
    
    print ('BIO model with branched axon')
    
    # Branching with Rall exponent 
    dend_diam2 = 6.*um 
    dend_length2 = 1000.*um 
    axon_diam2 = 1. * um
    axon_length2 = 100. * um 
    soma_diameter2 = 30. * um
    morpho2 = Soma(diameter=soma_diameter2)
    dendrite2 = Cylinder(diameter=dend_diam2, length=dend_length2, n=500)
    axon2 = Cylinder(diameter=axon_diam2, length=axon_length2, n=100)
    morpho2.dendrite = dendrite2
    morpho2.axon = axon2
    # axon branches
    branch_diam_rall = axon_diam2/(2**(2./3))
    morpho2.axon.branch1 = Cylinder(diameter=branch_diam_rall, length=1000*um, n=1000)
    morpho2.axon.branch2 = Cylinder(diameter=branch_diam_rall, length=1000*um, n=1000)
    
    # run the simulations with joblib
    if __name__ == '__main__':
        results_branched_rall = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC)(start, start+length,  GNa, morpho2) for start in starts)
        thresholds_branched_rall = array(results_branched_rall)*1e3
        
    # Branching with empirical exponent 
    dend_diam3 = 6.*um 
    dend_length3 = 1000.*um 
    axon_diam3 = 1. * um
    axon_length3 = 100. * um 
    soma_diameter3 = 30. * um
    morpho3 = Soma(diameter=soma_diameter3)
    dendrite3 = Cylinder(diameter=dend_diam3, length=dend_length3, n=500)
    axon3 = Cylinder(diameter=axon_diam3, length=axon_length3, n=100)
    morpho3.dendrite = dendrite3
    morpho3.axon = axon3
    # axon branches
    branch_diam_emp = axon_diam3/(2**(5./2))
    morpho3.axon.branch1 = Cylinder(diameter=branch_diam_emp, length=1000*um, n=1000)
    morpho3.axon.branch2 = Cylinder(diameter=branch_diam_emp, length=1000*um, n=1000)
    
    # run the simulations with joblib
    if __name__ == '__main__':
        results_branched_emp = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC)(start, start+length,  GNa, morpho3) for start in starts)
        thresholds_branched_emp = array(results_branched_emp)*1e3
    
    # Save the data in an npz file
    savez('figure_14c', starts/um, length/um, thresholds/mV, thresholds_branched_rall/mV, thresholds_branched_emp/mV, branch_diam_emp/um, branch_diam_rall/um)
    
x_mids = starts+length/2 # AIS middle position

# Slopes
slopes_branched_rall, _, _, _, _ = stats.linregress(-log(ra*x_mids* GNa), thresholds_branched_rall)
slopes_branched_emp, _, _, _, _ = stats.linregress(-log(ra*x_mids* GNa), thresholds_branched_emp)
slopes_500, _, _, _, _ = stats.linregress(-log(ra*x_mids* GNa), thresholds)

print ('Slopes:')
print ('Biophysical model in current-clamp:')
print('Threshold vs delta:', slopes_500, 'mV')
print('Threshold vs delta branched axon Rall:', slopes_branched_rall, 'mV')
print('Threshold vs delta branched axon empirical:', slopes_branched_emp, 'mV')

### Plots

fig = figure(2, figsize=(4,3))

ax3 = subplot(111)

ax3.semilogx(x_mids/um, thresholds, '--', color='k', label='original axon: k=%.2f mV' %slopes_500)
ax3.semilogx(x_mids/um, thresholds_branched_rall, '-', color='forestgreen', label='exp 3/2:$d_b=$%0.2f $\mu$m, k=%.2f mV' %(branch_diam_rall/um, slopes_branched_rall))
ax3.semilogx(x_mids/um, thresholds_branched_emp, '-', color='darkorange', label='exp 5/2:$d_b=$%0.2f $\mu$m, k=%.2f mV' %(branch_diam_emp/um,slopes_branched_emp))
ax3.set_xlabel('$x_{1/2}$ ($\mu$m)') 
ax3.set_ylabel('$V_s$ (mV)') 
ax3.set_ylim(-65,-50)
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax3.xaxis.set_minor_formatter(ScalarFormatter())
ax3.legend(frameon=False, fontsize=8)

tight_layout()

show()   
 






