
"""
AIS geometry and excitability, section 3: "Non-sodium axonal currents" (Figure 12)

Effect of a hyperpolarizing current on the resting membrane potential and the voltage threshold
in a point AIS. 

"""
from __future__ import print_function
from brian2 import *
from shared import params_all, model_Na_Kv1, measure_current_threshold, measure_voltage_threshold, calculate_resting_state
from joblib import Parallel, delayed
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.ticker as ticker

only_plotting = True # to plot the figure without running the simulations

# Parameters
defaultclock.dt = 0.005*ms
n = 5
m = 5
params = params_all
ra = 4.*params.Ri/(pi*params.axon_diam**2) # axial resistance per unit length
starts = linspace(10., 30., n)*um # AIS start positions
GNa = 500.*nS # total NA conductance in the AIS

### SIMULATIONS in the biophysical model

if only_plotting: # loading data
    data = load('figure_12.npz')
    starts = data['arr_0']*um
    inj_currents= data['arr_1']*nA
    thresholds_soma = data['arr_2']*1e-3
    thresholds_ais = data['arr_3']*1e-3
    v_rest_soma = data['arr_4']*1e-3
    v_rest_ais = data['arr_5']*1e-3
    delta_v_rest = data['arr_6']*1e-3
else: # running the simulations
    # A function to measure the threshold at the soma and at the AIS when a hyperpolarizing current is injected at the AIS
    def BIO_model_in_CC_Kv(resting_vm, ais_start, ais_end, i_inj, gna_tot):
        params = params_all
        defaultclock.dt = 0.005*ms
        pulse_length = 50.*ms

        # current threshold
        neuron = model_Na_Kv1(params, resting_vm, Na_start = ais_start, Na_end = ais_end, density=False, gna_tot=gna_tot)
        i_rheo = measure_current_threshold(params=params, neuron=neuron, resting_vm=resting_vm, ais_start=ais_start, \
                                           ais_end=ais_end, pulse_length=pulse_length, i_inj=i_inj)  
            
        # voltage threshold
        neuron = model_Na_Kv1(params, resting_vm, Na_start = ais_start, Na_end = ais_end, gna_tot = gna_tot, density=False)
        vs, va, _, _ = measure_voltage_threshold(params=params, neuron=neuron, resting_vm=resting_vm, ais_start=ais_start, \
                                           ais_end=ais_end, i_rheo = i_rheo, pulse_length=pulse_length, i_inj=i_inj)
        
        neuron = model_Na_Kv1(params, resting_vm, Na_start = ais_start, Na_end = ais_end, gna_tot=gna_tot,  density=False)
        t, v_soma, _, v_ais = calculate_resting_state(neuron, ais_start, ais_end, i_inj)
    
        return i_rheo, vs, va, t, v_soma, v_ais
    
    inj_currents = linspace(0., -400, m)*pA # injected hyperpolarizing current in the AIS
    
    # Measuring the thresholds
    if __name__ == '__main__':
        traces = Parallel(n_jobs = 5)(delayed(BIO_model_in_CC_Kv)(-75.*mV, start, start, i_inj, GNa) for start in starts for i_inj in inj_currents)
    
    # Re-arranging the array
    thresholds_soma = zeros((n,m))
    thresholds_ais = zeros((n,m))
    thresholds_ais_pulse = zeros((n,m))
    v_rest_soma = zeros((n,m))
    v_rest_ais = zeros((n,m))
    delta_v_rest = zeros((n,m))
    
    for i in range(n):
        for j in range(m):
            thresholds_soma[i,j] = traces[m*i+j][1]*1e3
            thresholds_ais[i,j] = traces[m*i+j][2]*1e3
            v_rest_soma[i,j] = mean(traces[m*i+j][4][int(10.*ms/defaultclock.dt):int(19.*ms/defaultclock.dt)])
            v_rest_ais[i,j] = mean(traces[m*i+j][5][int(10.*ms/defaultclock.dt):int(19.*ms/defaultclock.dt)])
            delta_v_rest[i,j] = v_rest_ais[i,j] - v_rest_soma[i,j]
    
    # Save the data in a npz file
    savez('figure_12', starts/um, inj_currents/nA, thresholds_soma/mV, thresholds_ais/mV, v_rest_soma/mV, v_rest_ais/mV, delta_v_rest/mV)
        
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

starts_label = starts/um
inj_curr_label = inj_currents/pA

fig_hyp = figure(2, figsize=(6,5))

# Panel A: resting membrane potential at the soma and AIS vs injected current
ax1 = subplot(221)
ax1.plot(abs(inj_currents)/pA, v_rest_soma[3,:], 'k', label = 'soma') 
ax1.plot(abs(inj_currents)/pA, v_rest_ais[3,:], 'k-.', label = 'AIS')  
ax1.set_ylim(-110, -70)
ax1.set_xlabel('| I |  (pA)')
ax1.set_ylabel('$V$ (mV)')
ax1.legend(frameon=False, fontsize=8)

ax1.text(-175, -70,'A', fontsize=14, weight='bold')

print ('panel A: AIS starts at %0.1f' %starts_label[3], 'um')

# Panel B: threshold at the soma and teh AIS vs injected current
ax4 = subplot(222)
ax4.plot(abs(inj_currents)/pA, thresholds_soma[3,:], 'k') 
ax4.plot(abs(inj_currents)/pA, thresholds_ais[3,:], 'k-.')  
ax4.set_yticks([-75,-65,-55,-45])
ax4.set_yticklabels(['-75','-65','-55','-45'])
ax4.set_ylim(-75, -45)
ax4.set_xlabel('| I |  (pA)')
ax4.set_ylabel('Threshold (mV)')

ax4.text(-160, -45,'B', fontsize=14, weight='bold')

print ('panel B: AIS starts at %0.1f' %starts_label[3], 'um')

# Panel C: difference in threshold vs difference in resting membrane potential
ax3 = subplot(223)
ax3.plot(linspace(delta_v_rest[-1,-1], delta_v_rest[0,0], 20), linspace(delta_v_rest[-1,-1], delta_v_rest[0,0], 20) + params.Ka/mV, '--', color='k', label='theory')
for i in range(int(m/2)+1):
    ax3.plot(delta_v_rest[:,2*i], thresholds_ais[:,2*i] - thresholds_soma[:,2*i], color=colors[2*i])
ax3.set_xlabel('$\Delta V$ (mV)')
ax3.set_ylabel('$\Delta$ Threshold (mV)')
ax3.set_ylim(-12.5, 7.5)
ax3.set_xlim(-17.5, 2.5)
ax3.legend(frameon=False, fontsize=8)

ax3.text(-24.5, 7.5,'C', fontsize=14, weight='bold')

# Panel D: threshold at the soma and AIS vs AIS position for different injected currents
ax2 = subplot(224)
for i in range(int(m/2)+1):
    semilogx(starts/um, thresholds_ais[:,2*i], '-.', color=colors[2*i])
    semilogx(starts/um, thresholds_soma[:,2*i], color=colors[2*i], label='%.0f pA' %inj_curr_label[2*i])
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.xaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
ax2.set_yticks([-75,-65,-55,-45])
ax2.set_yticklabels(['-75','-65','-55','-45'])
ax2.set_ylim(-75, -45)
#ax2.set_xlim(1,41)
ax2.set_xlabel('$\Delta$ ($\mu$m)')
ax2.set_ylabel('Threshold (mV)')
ax2.legend(frameon=False, fontsize=8)

ax2.text(6.5, -45,'D', fontsize=14, weight='bold')

tight_layout()

show()


