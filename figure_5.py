
'''

Presentation of a simple biophysical model of spike initiation in the AIS and how to use it.

'''

from brian2 import *
from shared.models import params_model_description, model_Na_Kv1
from shared.analysis.trace_analysis import *
import matplotlib.image as mpimg
import matplotlib.patches as patches

defaultclock.dt = 0.005*ms

### AIS geometry and parameters
start =5.*um # AIS start position
end = 35.*um # AIS end position
params = params_model_description

### Model
# The total conductances at the AIS is fixed:
#neuron = model_Na_Kv1(params, -75.*mV, start, end, density=False, gna_tot=350.*nS, gk_tot=150.*nS, morpho=params.morpho)

# The conductance densities at the AIS is fixed:
neuron = model_Na_Kv1(params, -75.*mV, start, end, density=True, morpho=params.morpho)

M = StateMonitor(neuron, ('v','I_VC', 'INa', 'IK'), record = 0) # recording at the soma
M_AIS = StateMonitor(neuron, ('v', 'm', 'h', 'n', 'INa', 'IK'), record = neuron.morphology.axon[end-1*um]) # recording at AIS end

### Simulation
# Protocol to elicit an AP with fixed initial voltage
neuron.V_VC[0] = -75.*mV
neuron.VC_on[0] = 1
run(20*ms)
neuron.VC_on[0] = 0
neuron.I_CC[0] = 1.5*nA
run(2*ms)
neuron.VC_on[0] = 0
neuron.I_CC[0] = 0*amp 
run(20*ms)

# Somatic AP features
print ('SOMA')

i_spike = spike_onsets(M.v[0], criterion = 30*volt/second * defaultclock.dt, v_peak=-20 * mV)

print ("Spike at time ", M.t[i_spike])
print ("Spike onset:", M.v[0][i_spike])
print ("Peak:", M.v[0][spike_peaks(M.v[0], v_peak=-20 * mV)[0]])
print ("Onset rapidness:", max_phase_slope(M.v[0], v_peak=0*mV)[0]/defaultclock.dt)
print ("Onset rapidness:", max_phase_slope(M.v[0], v_peak=-20*mV)[0]/defaultclock.dt)

rise_time, spike_width, _, _, _, _ = spike_duration(M.v[0], onsets = i_spike, full=True)
print ("Rise time:", rise_time*defaultclock.dt)
print ("Spike width:", spike_width*defaultclock.dt)
print ("Max dv/dt:", max_dv_dt(M.v[0], v_peak=-20 * mV)/defaultclock.dt)

half_height = (M.v[0][spike_peaks(M.v[0], v_peak=-20 * mV)[0]] + abs(M.v[0][i_spike]))/2 + M.v[0][i_spike]
print ('Half height:', half_height)
half_height_idx = where(M.v[0] > half_height )
half_width = (half_height_idx[0][-1] - half_height_idx[0][0])*defaultclock.dt
print ('Half width:', half_width)

# AIS AP features
print ('AIS')

i_spike_ais = spike_onsets(M_AIS.v[0], criterion = 50*volt/second * defaultclock.dt, v_peak=-20 * mV)

print ("Spike at time ",M_AIS.t[i_spike_ais[0]])
print ("Spike onset:",M_AIS.v[0][i_spike_ais[0]])
print ("Peak:", M_AIS.v[0][spike_peaks(M_AIS.v[0], v_peak=-20 * mV)[0]])
print ("Onset rapidness:", max_phase_slope(M_AIS.v[0], v_peak=0*mV)[0]/defaultclock.dt)
print ("Onset rapidness:", max_phase_slope(M_AIS.v[0], v_peak=-20*mV)[0]/defaultclock.dt)

rise_time_ais, spike_width_ais, _, _, _, _ = spike_duration(M_AIS.v[0], onsets = i_spike_ais, full=True)
print ("Rise time:", rise_time_ais*defaultclock.dt)
print( "Spike width:", spike_width_ais*defaultclock.dt)
print ("Max dv/dt:", max_dv_dt(M_AIS.v[0], v_peak=-20 * mV)/defaultclock.dt)

half_height = (M_AIS.v[0][spike_peaks(M_AIS.v[0], v_peak=-20 * mV)[0]] + abs(M_AIS.v[0][i_spike]))/2 + M_AIS.v[0][i_spike]
print ('Half height:', half_height)
half_height_idx = where(M_AIS.v[0] > half_height )
half_width = (half_height_idx[0][-1] - half_height_idx[0][0])*defaultclock.dt
print ('Half width:', half_width)

# Plotting
f1 = figure(1, figsize=(6, 10))

# Morphology
ax0 = subplot(411)
xlim(-12.5,12.5)
ylim(-4.5,4.5)
ax0.axis('off')
ax0.add_patch(patches.Circle((3.5,-1.7), 1.2, color='black')) # soma
ax0.add_patch(patches.Rectangle((-12.5, -2), 16, 0.6, color='black')) # dendrite
ax0.add_patch(patches.Rectangle((4.5, -1.7), 7, 0.1, color='black')) # axon
ax0.add_patch(patches.Rectangle((5, -1.7), 0.5, 0.1, color='red')) # AIS
ax0.text(-12.5, -3,'1000 $\mu$m') # dendrite length
ax0.text(-14.5, -1.95,'6 $\mu$m') # dendrite diam
ax0.text(9, -3,'500 $\mu$m') # axon length
ax0.text(11.75, -1.9,'1 $\mu$m') # axon diam
ax0.text(2.25, 0.2,'30 $\mu$m') # soma diam
ax0.text(5, -3,'AIS', color='red') # AIS label
ax0.text(5, -4,'L = 30 $\mu$m', color='red') 
ax0.text(5, -5,'$\Delta$ = 5 $\mu$m', color='red') # soma diam

# Panel A: equilibrium functions of the gating variables
ax1 = subplot(424)
ax1.plot(M_AIS.t/ms, M_AIS.m[0],  color='forestgreen', label='m')
ax1.plot(M_AIS.t/ms, M_AIS.h[0],  color='darkblue', label='h')
ax1.plot(M_AIS.t/ms, M_AIS.n[0]**8, 'darkorange', label='n$^8$')
ax1.set_xlim(20,23)
ax1.set_ylim(-0.05,1.05)
ax1.set_xlabel('Time (ms)')
ax1.set_xticks([20,21,22,23])
ax1.legend(frameon=False, fontsize=8)
ax1.text(19,1.05,'B', fontsize=14, weight='bold')

# Panel B: time course of the gating variables at AIS end
ax2 = subplot(423)
v = linspace(params.EL-25.*mV,20*mV,105)
minf = 1/(1+exp(-(v-params.Va)/params.Ka))
hinf = 1/(1+exp((v-params.Vh)/params.Kh))
ninf = 1/(1+exp(-(v-params.Vn)/params.Kn))
ninf8 = (1/(1+exp(-(v-params.Vn)/(params.Kn))))**8
ax2.plot(v/mV, minf, color='forestgreen', label='$m_\infty$')
ax2.plot(v/mV, hinf, color='darkblue', label='$h_\infty$')
ax2.plot(v/mV, ninf8,'darkorange',  label='$n_\infty^8$')
ax2.set_xlabel('V (mV)')#, fontsize=14)
ax2.set_xlim(-100,20)
ax2.set_ylim(-0.05,1.05)
ax2.legend(frameon=False, fontsize=8)
ax2.text(-140,1.05,'A', fontsize=14, weight='bold')

# Panel C: AP at AIS end and soma
ax3 = subplot(425)
ax3.plot(M.t/ms, M.v[0]/mV, 'k', label='soma') 
ax3.plot(M_AIS.t/ms, M_AIS.v[0]/mV, 'r', label='AIS end')
ax3.set_ylabel('V (mV)')
ax3.set_xlim(20,23)
ax3.set_ylim(-75, 50)
ax3.set_xlabel('Time (ms)')
ax3.set_xticks([20,21,22,23])
ax3.legend(frameon=False, fontsize=8)
ax3.text(19,50,'C', fontsize=14, weight='bold')

# Panel D: phase plot of the AP at AIS end and soma
ax4 = subplot(426)
t = M.t<30*ms
ax4.plot(M.v[0][t]/mV, diff(M.v[0])[t[:-1]]/defaultclock.dt, 'k')
ax4.plot(M_AIS.v[0][t]/mV, diff(M_AIS.v[0])[t[:-1]]/defaultclock.dt, 'r')
ax4.set_xlim(-75, 50)
ax4.set_ylim(-200,2500)
ax4.set_ylabel('dV/dt (V/s)')
ax4.set_xlabel('V (mV)')
ax4.text(-115, 2500,'D', fontsize=14, weight='bold')

# Panel E: Na and Kv1 currents at AIS
ax5 = subplot(427)
ax5.plot(M_AIS.t/ms, M_AIS.INa[0], 'r', label='Na')
ax5.plot(M_AIS.t/ms, -M_AIS.IK[0], 'r--', label='-K')
ax5.set_xlim(20,23)
ax5.set_ylim(-5, 100)
ax5.set_ylabel('Current (pA/$\mu$ m$^2$)')
ax5.set_xlabel('Time (ms)')
ax5.set_xticks([20,21,22,23])
ax5.legend(frameon=False, fontsize=8)
ax5.text(19, 100,'E', fontsize=14, weight='bold')

# Panel F: Na and Kv1 currents at soma
ax6 = subplot(428)
ax6.plot(M.t/ms, M.INa[0], 'k', label='Na')
ax6.plot(M.t/ms, -M.IK[0], 'k--', label='-K')
ax6.set_xlim(20,23)
ax6.set_ylim(-0.5, 10)
ax6.set_ylabel('Current (pA/$\mu$ m$^2$)')
ax6.set_xlabel('Time (ms)')
ax6.set_xticks([20,21,22,23])
ax6.legend(frameon=False, fontsize=8)
ax6.text(19, 10,'F', fontsize=14, weight='bold')

tight_layout()

show()
