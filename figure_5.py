

"""

We compare the rheobase with the somatic voltage threshold and 
how they vary with the leak conductance at the soma, the total Na conductance at the AIS and the amount of current injected at the AIS.

"""
from __future__ import print_function
from brian2 import *
from joblib import Parallel, delayed
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from shared.analysis import measure_current_threshold_no_VC, measure_voltage_threshold_no_VC, measure_input_resistance
from shared.models import params_model_description, model_Na_Kv1

only_plotting = True # to plot the figure without running the simulations

# Parameters
n=4 
defaultclock.dt = 0.005*ms
start = 5.*um # AIS start position
end = 35.*um # AIS end position
resting_vm = -75.*mV # resting membrane potential
params = params_model_description # model parameters

# A function to measure current and voltage thresholds
def measure_different_thresholds(gna_ais, i_inj, gL_soma):
    '''
    A function to measure the charge, current and voltage thresholds at the soma. 
    
    Variables: the somatic leak conductance, the total Na conductance at the AIS and the amount of hyperpolarizing current injected at the AIS.
    '''

    # current threshold
    neuron = model_Na_Kv1(params, resting_vm, start, end, density=False, gna_tot=gna_ais, gk_tot=params.Gk, gL_soma = gL_soma)
    i_rheo = measure_current_threshold_no_VC(params, neuron, resting_vm, start, end, pulse_length = 50.*ms, i_inj = i_inj, plot_v = False)  
    print ('Current threshold:', i_rheo)
    
    # voltage threshold
    neuron = model_Na_Kv1(params, resting_vm, start, end, density=False, gna_tot=gna_ais, gk_tot=params.Gk, gL_soma = gL_soma)
    vs, va, _, _, v0, v0_ais= measure_voltage_threshold_no_VC(params, neuron, resting_vm,start, end, i_rheo = i_rheo, pulse_length = 50.*ms, i_inj = i_inj) 
    print ('Threshold:', vs)
    
    # input resistance
    neuron = model_Na_Kv1(params, resting_vm, start, end, density=False, gna_tot=gna_ais, gk_tot=params.Gk)
    rs,_,_ = measure_input_resistance(params, neuron, resting_vm)
    print (rs)
    
    return i_rheo, vs, va, rs, v0, v0_ais

if only_plotting: # loading data to plot figure
    data = load('figure_5.npz')
    
    gL_somas = data['arr_0']
    gnas = data['arr_1']*nS
    hyp_curr = data['arr_2']*nA
    current_thresholds_gl = data['arr_3']*nA
    current_thresholds_gna = data['arr_4']*nA
    current_thresholds_hyp = data['arr_5']*nA
    v_thresholds_gl = data['arr_6']*mV
    v_thresholds_gna = data['arr_7']*mV
    v_thresholds_hyp = data['arr_8']*mV
    v_thresholds_ais_hyp = data['arr_9']*mV
    v0_soma_gl = data['arr_10']*mV
    v0_ais_gl = data['arr_11']*mV
    v0_soma_gna = data['arr_12']*mV
    v0_ais_gna = data['arr_13']*mV
    v0_soma_hyp = data['arr_14']*mV
    v0_ais_hyp = data['arr_15']*mV
    
else: # running simulations
    # Variation of somatic leak conductance density
    print ("changing gL")
    gL_somas = linspace(1./40000., 1./1000., n) * (1./(ohm * cm**2))
    
    # run simulations with joblib
    if __name__ == '__main__':
        results_gl = Parallel(n_jobs = 4)(delayed(measure_different_thresholds)( 350.*nS, 0, gLs) for gLs in gL_somas)
    
    v_thresholds_gl = [results_gl[i][1] for i in range(n)] # voltage threshold at soma
    current_thresholds_gl = [results_gl[i][0] for i in range(n)] # current threshold
    rin_soma_gl = [results_gl[i][3] for i in range(n)] # input resistance at the soma
    v0_soma_gl = [results_gl[i][4] for i in range(n)] # resting membrane potential at the soma
    v0_ais_gl = [results_gl[i][5] for i in range(n)] # resting membrane potential at the AIS
    
    # Variation of Na total conductance at AIS
    print ("changing Na conductance")
    gnas = linspace(200.,500.,n)*nS
    
    # run simulations with joblib
    if __name__ == '__main__':
        results_gna = Parallel(n_jobs = 4)(delayed(measure_different_thresholds)( gna, 0, params.gL) for gna in gnas)
    
    v_thresholds_gna = [results_gna[i][1] for i in range(n)]# voltage threshold at soma
    current_thresholds_gna = [results_gna[i][0] for i in range(n)]# current threshold
    rin_soma_gna = [results_gna[i][3] for i in range(n)]# input resistance at the soma
    v0_soma_gna = [results_gna[i][4] for i in range(n)]# resting membrane potential at the soma
    v0_ais_gna = [results_gna[i][5] for i in range(n)]# resting membrane potential at the AIS
    
    # Variation of hyperpolarizing current injected at AIS end
    print ("changing AIS current")
    hyp_curr = linspace(-0.25, 0,n)*nA
    
    # run simulations with joblib
    if __name__ == '__main__':
        results_hyp = Parallel(n_jobs = 4)(delayed(measure_different_thresholds)( 450.*nS, i_hyp, params.gL) for i_hyp in hyp_curr)
    
    v_thresholds_hyp = [results_hyp[i][1] for i in range(n)]# voltage threshold at soma
    current_thresholds_hyp = [results_hyp[i][0] for i in range(n)]# current threshold
    v_thresholds_ais_hyp = [results_hyp[i][2] for i in range(n)]# voltage threshold at the AIS
    rin_soma_hyp = [results_hyp[i][3] for i in range(n)]# input resistance at the soma
    v0_soma_hyp = [results_hyp[i][4] for i in range(n)]# resting membrane potential at the soma
    v0_ais_hyp = [results_hyp[i][5] for i in range(n)]# resting membrane potential at the AIS
    
    # Save the data in an npz file
    savez('figure_6', gL_somas/(siemens/meter**2), gnas/nS, hyp_curr/nA, current_thresholds_gl/nA, current_thresholds_gna/nA, current_thresholds_hyp/nA,\
                  v_thresholds_gl/mV, v_thresholds_gna/mV, v_thresholds_hyp/mV,\
                  v_thresholds_ais_hyp/mV, v0_soma_gl/mV, v0_ais_gl/mV, v0_soma_gna/mV, v0_ais_gna/mV, v0_soma_hyp/mV, v0_ais_hyp/mV)

### Figure 

fcurrent = figure('Voltage threshold vs rheobase', figsize=(9.5,3))  

# Panel A: voltage and current threshold vs G
i1_range = max(current_thresholds_gna)/nA * (-75.*mV - (-35.*mV))/(resting_vm - max(v_thresholds_gna))
G_min = 200.
G_max = 500.
ax1 = fcurrent.add_subplot(131)
ax2 = ax1.twinx()

ax1.plot(gnas/nS, v_thresholds_gna/mV, 'cornflowerblue')
ax1.plot(gnas/nS, v0_soma_gna/mV, 'k--')
ax1.tick_params(axis='y', colors='cornflowerblue')
ax1.set_ylim(-75, -35)
ax1.set_xlim(G_min,G_max)
ax1.set_ylabel('$V_{s}$ (mV)', color='cornflowerblue')
ax1.set_xlabel('$G$ (nS)')
ax1.text(G_min+(G_max-G_min)/15, -72,'$V_0$', fontsize=10)

ax2.plot(gnas/nS, current_thresholds_gna/nA, color= 'darkblue')
ax2.set_ylim(0, i1_range)
ax2.set_ylabel('$I_{r}$ (nA)', color='darkblue') 
ax2.tick_params(colors='darkblue')

ax1.text((G_max-G_min)/3, -35,'A', fontsize=14, weight='bold')

# Panel B: voltage and current threshold vs gL
i2_range = min(current_thresholds_gl)/nA * (-75.*mV - (-35.*mV))/(resting_vm - max(v_thresholds_gl))
gl_min = 0.25
gl_max = 10.
ax3 = fcurrent.add_subplot(132)
ax4 = ax3.twinx()

ax3.plot(gL_somas, v_thresholds_gl/mV, 'cornflowerblue')
ax3.plot(gL_somas, v0_soma_gl/mV, 'k--')
ax3.set_ylim(-75, -35)
ax3.set_ylabel('$V_{s}$ (mV)', color='cornflowerblue')
ax3.set_xlabel('$g_L$ (S/$m^2$)')
ax3.tick_params(axis='y', colors='cornflowerblue')
ax3.text(gl_min+(gl_max-gl_min)/15, -72,'$V_0$', fontsize=10)

ax4.plot(gL_somas, current_thresholds_gl/nA, 'darkblue')
ax4.set_ylabel('$I_{r}$ (nA)', color='darkblue') 
ax4.tick_params(colors='darkblue')
ax4.set_xlim(gl_min,gl_max)
ax4.set_ylim(0, i2_range)

ax3.text(-(gl_max-gl_min)/3, -35,'B', fontsize=14, weight='bold')

# Panel C: somatic and axonal voltage threshold and current threshold vs Ih
i3_range = min(current_thresholds_hyp)/nA * (-75.*mV - (-35.*mV))/(-75.*mV - min(v_thresholds_hyp))
ih_max = 0
ih_min = -0.25
ax5 = fcurrent.add_subplot(133)
ax6 = ax5.twinx()

ax5.plot(hyp_curr/nA, v_thresholds_hyp/mV, 'cornflowerblue', label='soma')
ax5.plot(hyp_curr/nA, v_thresholds_ais_hyp/mV, '-.', color='cornflowerblue', label='AIS')
ax5.set_ylim(-75, -35)
ax5.set_xlim(ih_min, ih_max)
ax5.set_ylabel('Threshold (mV)', color='cornflowerblue')
ax5.set_xlabel('I (nA)')
ax5.tick_params(axis='y', colors='cornflowerblue')
ax5.legend(frameon=False, fontsize=8, loc='lower left')

ax6.plot(hyp_curr/nA, current_thresholds_hyp/nA, 'darkblue') 
ax6.set_ylim(0,  i3_range)
ax6.set_ylabel('$I_{r}$ (nA)', color='darkblue') 
ax6.tick_params(colors='darkblue')

ax5.text(ih_min + (ih_min-ih_max)/3, -35,'C', fontsize=14, weight='bold')

subplots_adjust(top=0.88, bottom=0.29, left=0.075,right=0.93,hspace=0.2,wspace=0.65)

show()


