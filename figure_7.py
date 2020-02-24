
"""
AIS geometry and excitability, section 3: "A spatially extended AIS" (figure 8)


We look at the threshold as a function of the AIS length and position in the biophysical model,
for a fixed Na conductance density at the AIS.

"""
from __future__ import print_function
from brian2 import *
from scipy import stats
from joblib import Parallel, delayed
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from shared.models import model_Na_Kv1, params_model_description
from shared.analysis import measure_current_threshold, measure_voltage_threshold

only_plotting = True # to plot the figure without running the simulations

# Parameters
m = 5
defaultclock.dt = 0.005*ms
params = params_model_description #params_all

### SIMULATIONS in the biophysical model

if only_plotting: # loading data
    data = load('figure_7.npz')
    starts = data['arr_0']*um # AIS start positions
    lengths = data['arr_1']*um # AIS length
    thresholds = data['arr_2']*mV
    
else: # running simulations
    starts = linspace(10., 30., m)*um # AIS start positions
    lengths = linspace (10., 40., m)*um # AIS length

    def BIO_model_in_CC_constant_density(ais_start, ais_end):
        defaultclock.dt = 0.005*ms
        resting_vm = -75.*mV
        pulse_length = 50.*ms
        print ('AIS start:', ais_start, 'AIS end:', ais_end)
            
        # current threshold
        neuron = model_Na_Kv1(params, resting_vm, ais_start, ais_end, density = True)
        i_rheo = measure_current_threshold(params, neuron, resting_vm, ais_start, ais_end, pulse_length = pulse_length)  
        print ('Rheobase:', i_rheo)
        
        # voltage threshold
        neuron = model_Na_Kv1(params, resting_vm, ais_start, ais_end, density=True)
        vs, va, _, _ = measure_voltage_threshold(params, neuron, resting_vm, ais_start, ais_end, i_rheo = i_rheo, pulse_length = pulse_length) 
        print ('Threshold:', vs)
    
        return vs
    
    if __name__ == '__main__':
        results = Parallel(n_jobs = 5)(delayed(BIO_model_in_CC_constant_density)(start, start+length) for start in starts for length in lengths)
    
    thresholds = array(results)
    thresholds = thresholds.reshape((m,m))*1e3
    
    # Save the data in an npz file
    savez('figure_8', starts/um, lengths/um, thresholds/mV)
        
### THEORY
ra = 4*params.Ri/(pi*params.axon_diam**2)
y0 = params.y0

### Threshold formula
def nu0_prediction(start, length):
    end = start + length
    return 2.*log(sqrt(2.)*y0*(length/end))-2.*log(cosh(y0*length/end))-2.*y0*(start/end)*tanh(y0*length/end)

def V0_prediction(start, length):
    vth = params.Va + params.Ka*nu0_prediction(start, length) - params.Ka*log(length*ra*params.gna_dens*length*pi*params.axon_diam*(params.ENa-params.Va)/params.Ka)
    return vth

thresholds_pred = zeros((m,m))
for i in range(m):
    for j in range(m):
        thresholds_pred[i,j] = V0_prediction(starts[i], lengths[j]) * 1e3

### Plots
        
# Custom colormap
top = get_cmap('Oranges_r', 128)
bottom = get_cmap('Blues', 128)
newcolors = vstack((top(linspace(0.5, 1, 128)),
                       bottom(linspace(0.5, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')
cmap = plt.get_cmap(newcmp)
colors = [cmap(i) for i in np.linspace(0, 1, m)]

# Figure
lengths_range_label = lengths/um
starts_range_label = starts/um

fig_pred = figure(3, figsize=(6,5))

# Panel A: BIO model: threshold vs AIS start position, for different AIS lengths
ax1 = fig_pred.add_subplot(221)
ax1.text(26.5, -42,'Biophysical model', fontsize=12)
for k in range(m):
    plot(starts/um, thresholds[:,k], color=colors[k]) 
ax1.set_ylabel('$V_s$ (mV)')
ax1.set_ylim(-75,-45)
ax1.text(15,-49,'$L$ = %i $\mu$m' %lengths_range_label[0], fontsize=10, color=colors[0])
ax1.text(22,-68,'$L$ = %i $\mu$m' %lengths_range_label[m-1], fontsize=10, color=colors[m-1])

ax1.text(1.5, -47.25,'A', fontsize=14, weight='bold')

# Panel B: BIO model: threshold vs AIS length, for different AIS start positions
ax2 = fig_pred.add_subplot(222)
for k in range(m):
    plot(lengths/um, thresholds[k,:], color =colors[k]) 
ax2.set_ylim(-75,-45)
ax2.text(15,-49,'$\Delta$ = %i $\mu$m' %starts_range_label[0], fontsize=10, color=colors[0])
ax2.text(27.5,-68,'$\Delta$ = %i $\mu$m' %starts_range_label[m-1], fontsize=10, color=colors[m-1])

ax2.text(-2.5, -47.25,'B', fontsize=14, weight='bold')

# Panel C: theory: threshold vs AIS start position, for different AIS lengths
ax3 = fig_pred.add_subplot(223)
ax3.text(32., -42,'Theory', fontsize=12)
for k in range(m):
    plot(starts/um, thresholds_pred[:,k], color=colors[k]) 
ax3.set_xlabel('$\Delta$ ($\mu$m)')
ax3.set_ylabel('$V_s$ (mV)')
ax3.set_ylim(-75,-45)

ax3.text(1.5, -47.25,'C', fontsize=14, weight='bold')

# Panel C: theory: threshold vs AIS length, for different AIS start positions
ax4 = fig_pred.add_subplot(224)
for k in range(m):
    plot(lengths/um, thresholds_pred[k,:], color =colors[k]) 
ax4.set_xlabel('$L$ ($\mu$m)')
ax4.set_ylim(-75,-45)

ax4.text(-2.5, -47.25,'D', fontsize=14, weight='bold')

subplots_adjust(hspace=0.35, wspace=0.4)

show()









