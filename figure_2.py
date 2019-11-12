'''
Measurement of passive responses to current pulses in a ball-and-stick model.

Input resistances are calculated.

The neuron is a sphere + a long cylindrical axon. The diameter of the sphere is varied.

'''
from __future__ import print_function
from brian2 import *
from joblib import Parallel, delayed
import matplotlib.patches as patches
from shared.resistive_coupling.resistive_coupling_simulation import *

fig = figure(1, (8,4))

### SMALL NEURON
soma_diameter_small = axon_diameter

# Top panel: drawing of the axon 
ax1 = subplot(321)
title('Small soma')
xlim(-10,200)
ylim(-20,20)
ax1.axis('off')
ax1.add_patch(patches.Rectangle((-10,-5),210, 10, color='darkgrey'))
ax1.annotate('', xy=(20, -5), xytext=(20, -20),
            arrowprops=dict(facecolor='black', shrink=0.05))
ax1.annotate('', xy=(100, -5), xytext=(100, -20),
            arrowprops=dict(facecolor='red', shrink=0.05))

# Panel A: voltage across the axon at the end of the pulse
ax2 = subplot(323)

t, x, _ , _, v20 = pulse_simulation(soma_diameter_small, 20*um, current = 10*pA)
t, x, _, _, v100 = pulse_simulation(soma_diameter_small, 100*um, current = 10*pA)

ax2.plot(x/um, v20/mV,'k')
ax2.plot(x/um, v100/mV,'r')
ax2.set_xlim(-10,200)
ax2.set_ylim(-75,-60)
ax2.set_ylabel('Va (mV)')
ax2.text(-55,-60,'A', fontsize=14, weight='bold')

# Panel B: Input resistance vs. distance
ax3 = subplot(325)
distances = linspace(10,200,10)*um

# Simulation
results = Parallel(n_jobs = 4)(delayed(input_resistance)(soma_diameter_small,d) for d in distances)
_, _, Ra200ms, _ = zip(*results)

ax3.plot(distances/um,Ra200ms/Mohm,"k", label='R')
ax3.set_ylim(0,1000)
ax3.set_xlim(-10,200)
ax3.text(-55, 1000,'B', fontsize=14, weight='bold')

# Theory
space_constant = sqrt(axon_diameter/(gL*4*Ri))
ra = 4*Ri/(pi*axon_diameter**2)
#accurate_prediction = ra*space_constant/(1+tanh(distances/space_constant))
approximation = ra*(space_constant-distances)
ax3.plot(distances/um, approximation/Mohm,"k--", label='theory')
ax3.legend(frameon=False, fontsize=8)
ax3.set_xlabel('Distance ($\mu$m)')
ax3.set_ylabel('R (M$\Omega$)')

### LARGE NEURON
soma_diameter_large = 100*um
current = 100*pA

# Top panel: drawing of the axon
ax4 = subplot(322)
title('Large soma')
xlim(-10,200)
ylim(-20,20)
ax4.axis('off')
ax4.add_patch(patches.Rectangle((0, -5), 200, 10, color='darkgrey'))
ax4.add_patch(patches.Circle((-50, 0), 50, color='darkgrey'))
ax4.annotate('', xy=(20, -5), xytext=(20, -20),
            arrowprops=dict(facecolor='black', shrink=0.05))
ax4.annotate('', xy=(100, -5), xytext=(100, -20),
            arrowprops=dict(facecolor='red', shrink=0.05))

# Panel C: voltage across the axon at the end of the pulse
ax5 = subplot(324)

t, x, vs20 , va20, v20 = pulse_simulation(soma_diameter_large, 20*um, current = current)
t, x, vs100, va100, v100 = pulse_simulation(soma_diameter_large, 100*um, current = current)

ax5.plot(x/um, v20/mV,'k')
ax5.plot(x/um, v100/mV,'r')
ax5.set_xlim(-10,200)
ax5.set_ylim(-75,-60)
ax5.text(-50, -60,'C', fontsize=14, weight='bold')

# Panel D: Input resistance vs. distance
ax6 = subplot(326)
distances = linspace(10,200,10)*um

# Simulation
results = Parallel(n_jobs = 4)(delayed(input_resistance)(soma_diameter_large,d) for d in distances)
Ra300us, Rs300us, Ra200ms, _ = zip(*results)

ax6.plot(distances/um, Ra200ms/Mohm,"k", label='R')

# Theory
space_constant = sqrt(axon_diameter/(gL*4*Ri))
ra = 4*Ri/(pi*axon_diameter**2)
rs = 1/(pi*soma_diameter_large**2*gL) # soma resistance
rdistal = ra*space_constant # distal resistance
approximation = 1./(1/(ra*distances + rs) + 1./rdistal)
approximation2 = 1./(1/(ra*distances))# + 1/rdistal)
ax6.plot(distances/um, approximation/Mohm,"k--", label='theory (soma)')
ax6.plot(distances/um, approximation2/Mohm,"r--", label='theory (killed end)')
ax6.set_xlabel('Distance ($\mu$m)')
ax6.set_ylim(0,300)
ax6.set_xlim(-10,200)
ax6.text(-50, 300,'D', fontsize=14, weight='bold')

# legends
lines = ax6.get_lines()
legend1 = legend([lines[i] for i in [0]], ['R'], loc='lower right', frameon=False, fontsize=8)
legend2 = legend([lines[i] for i in [1,2]], [ 'theory (soma)', 'theory (killed end)'], loc='upper left', frameon=False, fontsize=8)
ax6.add_artist(legend1)
ax6.add_artist(legend2)

tight_layout()

show()
