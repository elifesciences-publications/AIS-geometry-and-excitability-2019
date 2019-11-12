

"""
AIS geometry and excitability, section 3: "Non-sodium axonal currents" (Figure 11)

Extended AIS with Kv7 channels in the distal half of the AIS.

"""

from brian2 import *
from shared import params_all, model_Na_Kv1_Kv7, measure_current_threshold, measure_voltage_threshold
from joblib import Parallel, delayed
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.ticker as ticker

only_plotting = True # to plot the figure without running the simulations

# Parameters
defaultclock.dt = 0.005*ms

n=7
params = params_all
ra = 4.*params.Ri/(pi*params.axon_diam**2) # axial resistance per unit length
length = 30.*um # AIS length
GNa = 500.*nS # total Na conductance at the AIS

### SIMULATIONS in the biophysical model

if only_plotting: # loading data
    data = load('figure_11.npz')
    starts = data['arr_0']*um
    gms = data['arr_1']*nS
    thresholds_soma = data['arr_2']*mV
    thresholds_ais = data['arr_3']*mV

    x_mids = starts+length/2 # AIS middle position
else: # running simulations
    # A function to measure the threshold as a function of AIS position, length and Kv7 conductance at the AIS
    def BIO_model_in_CC_Kv7(resting_vm, ais_start, ais_end, gna_tot, gm_tot):
        params = params_Kv7
        defaultclock.dt = 0.005*ms
        pulse_length = 50.*ms
        
        # current threshold
        neuron = model_Na_Kv1_Kv7(params, resting_vm, Na_start = ais_start, Na_end = ais_end, density=False, gna_tot=gna_tot, gm_tot = gm_tot)
        i_rheo = measure_current_threshold(params=params, neuron=neuron, resting_vm=resting_vm, ais_start=ais_start, ais_end=ais_end,pulse_length = pulse_length)  
            
        # voltage threshold
        neuron = model_Na_Kv1_Kv7(params, resting_vm, Na_start = ais_start, Na_end = ais_end, density=False, gna_tot=gna_tot, gm_tot = gm_tot)
        vs, va, _, _ = measure_voltage_threshold(params=params, neuron=neuron, resting_vm=resting_vm, ais_start=ais_start, ais_end=ais_end,\
                                                 i_rheo = i_rheo, pulse_length = pulse_length)
        
        return i_rheo, vs, va

    starts = linspace(0., 30., n)*um # AIS start positions
    gm_densities = linspace(0, 300., n)*(siemens/meter**2) # Kv7 conductance density at the AIS
    gms = gm_densities * (pi*(length/2)*params.axon_diam) # corresponding total Kv7 conductance at the AIS
    
    # Measuring the threshold
    if __name__ == '__main__':
        traces = Parallel(n_jobs = 5)(delayed(BIO_model_in_CC_Kv7)(-75.*mV, start, start+length, GNa, gm) for start in starts for gm in gms)
    
    thresholds_soma = zeros((n,n)) # voltage threshold at the soma
    thresholds_ais = zeros((n,n)) # voltage threshold at the AIS
    
    # Re-arranging the array
    for i in range(n):
        for j in range(n):
            thresholds_soma[i,j] = traces[n*i+j][1]*1e3
            thresholds_ais[i,j] = traces[n*i+j][2]*1e3
    
    # Save the data to a npz file
    savez('figure_11', starts/um, gms/nS, thresholds_soma/mV, thresholds_ais/mV)
    
### Plots
    
# Custom colormap
top = get_cmap('Oranges_r', 128)
bottom = get_cmap('Blues', 128)
newcolors = vstack((top(linspace(0.5, 1, 128)),
                       bottom(linspace(0.5, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')
cmap = plt.get_cmap(newcmp)
colors = [cmap(i) for i in np.linspace(0, 1, n/2+1)]

gms_label = gms/(pi*(length/2)*params.axon_diam) / (siemens/meter**2) 

# Figure
fig_kv7 = figure(2, figsize=(6,3))

# Panel A: somatic threshold vs AIS middle position
ax1 = subplot(121)
for i in range(int(n/2)+1):
    semilogx(x_mids/um, thresholds_soma[:,2*i], color=colors[i], label='%.0f S/$m^2$' %gms_label[2*i])
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.xaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
ax1.set_ylim(-75, -40)
ax1.set_xlabel('$x_{1/2}$ ($\mu$m)')
ax1.set_ylabel('$V_s$ (mV)')

# legends
lines = ax1.get_lines()
legend1 = legend([lines[i] for i in [0,1]], ['0 S/$m^2$', '100 S/$m^2$', \
                 '150 S/$m^2$'], loc='lower left', frameon=False, fontsize=8)
legend2 = legend([lines[i] for i in [2,3]], ['200 S/$m^2$', '300 S/$m^2$'], loc='lower right', frameon=False, fontsize=8)
ax1.add_artist(legend1)
ax1.add_artist(legend2)
 
ax1.text(9.5, -40,'A', fontsize=14, weight='bold')

# Panel B: threshold at the end of the AIS vs AIS middle position
ax2 = subplot(122)
for i in range(int(n/2)+1):
    semilogx(x_mids/um, thresholds_ais[:,2*i], color=colors[i])
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.xaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
#ax2.set_yticks([-75,-65,-55,-45,-35])
#ax2.set_yticklabels(['-75','-65','-55','-45','-35'])
ax2.set_ylim(-75, -40)
ax2.set_xlabel('$x_{1/2}$ ($\mu$m)')
ax2.set_ylabel('$V_a$ (mV)')
ax2.legend(frameon=False, fontsize=8)

ax2.text(9.5, -40,'B', fontsize=14, weight='bold')

subplots_adjust(top=0.88, bottom=0.29, wspace=0.4)

show()






