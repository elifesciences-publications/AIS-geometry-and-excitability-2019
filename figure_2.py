
'''

Pulse responses at the axon and soma, in a ball-and-stick model in the large soma regime and in a L5 pyramidal neuron (Hu and Bean 2018).

'''

from __future__ import print_function
from pylab import *
from joblib import Parallel, delayed
from shared.resistive_coupling.hu_bean_shared import *
from shared.resistive_coupling.resistive_coupling_simulation import *

### DATA from Hu and Bean 2018

## Example traces (injection point at 75 um)

dt = 0.02 # in ms

n_before = 50*5 # number of time steps before the pulse (= 50 kHz * duration in ms)
n_after = 50*100 # number of time steps after the pulse

distance,t,Vs,Va,Is,Ia = load_file('data/hu_and_bean_2018/100420 1-1 75um.abf')


# Find places where there's an axonal but not somatic current pulse
i = 1*((abs(Ia)>0.008) & (abs(Is)<0.008))
start = where(diff(i) == 1)[0] +1
end = where(diff(i) == -1)[0]

# Look for the smallest positive current pulse
min_current = 1e10
for i,j in zip(start,end):
    if j-i>24000: # long pulse
        current = median(Ia[i:j]) # average current
        if (current>0) & (current<min_current):
            vs = Vs[i-n_before:j+n_after]
            va = Va[i-n_before:j+n_after]
            min_current = current
            
current_example = min_current

print ("EXAMPLE CELL:")
print ("Smallest positive current pulse at the axonal bleb:", current, 'nA')

# Resting potential calculated just before the pulse (not used)
vs_rest = median(vs[:n_before])
va_rest = median(va[:n_before])

## Rin at the start of the pulse

distances = []
Ra_short, Rs_short = [], []

dt = 0.02 # in ms

n_before = 50*5 # number of time steps before the pulse (= 50 kHz * duration in ms)

print ('INPUT RESISTANCE vs DISTANCE')
for distance,t,Vs,Va,Is,Ia in load_all_files():
    print ("Distance from soma to axonal bleb: ",distance,"um")
    # Find places where there's an axonal but not somatic current pulse
    i = 1*((abs(Ia)>0.008) & (abs(Is)<0.008))
    start = where(diff(i) == 1)[0] +1
    end = where(diff(i) == -1)[0]

    vss = []
    vaa = []
    current = []
    for i,j in zip(start,end):
        if j-i>24000: # long pulse
            current.append(median(Ia[i:j]))
            vss.append(Vs[i-n_before:j])
            vaa.append(Va[i-n_before:j])

    current = array(current)

    try:
        # Largest positive current
        ind = where(current>0)[0]

        i=ind[argmin(current[ind])]
        current_pos = current[i]
        vs_pos = vss[i]
        va_pos = vaa[i]

        print("Current pulse injected at the bleb:", current_pos, 'nA')

    except ValueError:
        print ("No current pulse found")
        continue

    n = len(vs_pos)

    # Rest
    vs_rest = median(vs_pos[:n_before])
    va_rest = median(va_pos[:n_before])
    # Resistances
    R_input_short = (median(va_pos[n_before+10:n_before+20])-va_rest)/current_pos
    Rs_input_short = (median(vs_pos[n_before+10:n_before+20])-vs_rest)/current_pos
    R_axial = R_input_short - Rs_input_short
    print ("Input resistance at axonal bleb:", R_input_short)
    print ("Input resistance at soma:", Rs_input_short) 
    print ("Axial resistance:", R_axial)
    distances.append(distance)
    Ra_short.append(R_input_short)
    Rs_short.append(Rs_input_short)

# Print sorted list
#for d,Rp, Rn in sorted(zip(array(distances),Ra_short,Rs_short)):
    #print (d,'um:',Rp, Rn)

print ("Number of recordings used in the analysis:", len(Ra_short))

### SIMULATION in the large neuron regime
soma_diameter = 100*um
current = current_example * 1e3 * pA # the same current pulse as in the example
t_sim, x, vs75 , va75, v75 = pulse_simulation(soma_diameter, 75*um, current = current)

### Figure

fig_data = figure(1, figsize=(8,4))

### SIMULATIONS

# Panel A: traces for a pulse at 75 um
ax1 = subplot(321)
ax1.plot(t_sim/ms, va75/mV,'r',label='axon')
ax1.plot(t_sim/ms, vs75/mV,'k',label='soma')
ax1.set_title('Model')
ax1.set_ylabel('V (mV)')
ax1.legend(frameon=False, fontsize=8)
ax1.set_ylim(-76, -66)
ax1.set_xlim(0,300)
ax1.text(-60, -65,'A', fontsize=14, weight='bold')

# Panel B: the voltage gradient across the axon
ax2 = subplot(323)
ax2.plot(t_sim/ms,(va75-vs75)/mV,'forestgreen')
ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('$\Delta$ V (mV)')
ax2.set_ylim(-0.5,5)
ax2.set_xlim(0,300)
ax2.text(-60, 5.55,'B', fontsize=14, weight='bold')

# Panel C: input (or axial) resistance at axon and soma vs distance, at 300 us
distances_sim = linspace(10,200,10)*um
results = Parallel(n_jobs = 4)(delayed(input_resistance)(soma_diameter,d) for d in distances_sim)

Ra300us, Rs300us, Ra200ms, _ = zip(*results)
ax3 = subplot(325)
Ra300us, Rs300us = array(Ra300us), array(Rs300us)
ax3.plot(distances_sim/um, Ra300us/Mohm, 'r', label='axon')
ax3.plot(distances_sim/um, Rs300us/Mohm, 'k', label='soma')
ax3.set_xlabel('Distance ($\mu$m)')
ax3.set_ylabel('R (M$\Omega$)')
ax3.set_xlim(0,210)
ax3.set_ylim(-5,120)
ax3.text(-42, 132,'C', fontsize=14, weight='bold')

### DATA
# Panel D: traces for a pulse at 75 um
ax4 = subplot(322)
t=arange(len(va))*dt
ax4.plot(t,va, 'r', label='axon')
ax4.plot(t,vs, 'k', label='soma')
ax4.set_title('L5 pyramidal cell')
ax4.set_ylim(-61,-51)
ax4.set_xlim(-25,600)
ax4.text(-110, -50,'D', fontsize=14, weight='bold')

# Panel E: the voltage gradient across the axon
ax5 = subplot(324)
ax5.plot(t,va-vs, 'forestgreen')
ax5.set_xlabel('Time (ms)')
ax5.set_ylim(-0.5,3.5)
ax5.set_xlim(-25,600)
ax5.text(-110, 3.9,'E', fontsize=14, weight='bold')

# Panel C: input (or axial) resistance at axon and soma vs distance, at 300 us
ax6 = subplot(326)
ax6.plot(distances, Ra_short, 'r.',label='axon')
ax6.plot(distances, Rs_short, 'k.',label='soma')
ax6.set_xlabel('Distance ($\mu$m)')
ax6.set_xlim(0,210)
ax6.set_ylim(-5,60)
ax6.text(-28, 66.5,'F', fontsize=14, weight='bold')

tight_layout()
subplots_adjust(hspace=0.7)

show()

