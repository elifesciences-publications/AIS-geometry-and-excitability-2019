
"""

AIS geometry and excitability, section 3: "A point AIS" (figure 7)

The threshold measured:
    in voltage-clamp in a simple spike initiation model (SI model), 
    in the biophysical model in current-clamp (BIO model)
    with the theory

We show that the threshold goes as log(x gna).

"""

from brian2 import *
from shared.models import model_spike_initiation,  model_Na_Kv1, params_all
from shared.analysis import measure_threshold_in_vc, measure_current_threshold, measure_voltage_threshold
from shared.models.params_all import *
from scipy import stats
from joblib import Parallel, delayed
from mpl_toolkits.axes_grid1 import SubplotDivider, LocatableAxes, Size
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

only_plotting = True # to plot the figure without running the simulations

# Parameters
defaultclock.dt = 0.005*ms
n=5
ra = 4.*Ri/(pi*axon_diam**2) # axial resistance per unit length
starts = linspace(10., 30., n)*um # AIS starts
gnas = linspace (200, 600, n)*nS # total Na conductances

### SIMULATIONS

if only_plotting: # loading data
    data = load('figure_7.npz')
    starts = data['arr_0']*um
    gnas = data['arr_1']*nS
    thresholds_si = data['arr_2']*mV
    thresholds_bio = data['arr_3']*mV

else: # running simulations
    # A function to measure of the threshold in the SI model in VC, for different AIS start positions and Gna
    print ('SI model')
    def SI_model_threshold_in_VC(ais_start, ais_end, gna_tot):
        print (ais_start, gna_tot)
        params = params_all
        defaultclock.dt = 0.005*ms
        
        neuron = model_spike_initiation(params, ais_start, ais_end, gna_tot)
        _, v_thres = measure_threshold_in_vc(params, neuron, ais_start, ais_start, hyperpol_current = False, plot_v = True)
        print ('threshold:', v_thres)
        return v_thres
    
    # run the simulations with joblib
    if __name__ == '__main__':
        results_si = Parallel(n_jobs = 5)(delayed(SI_model_threshold_in_VC)(start, start, gna) for start in starts for gna in gnas)
        #print results_si
        thresholds_si = array(results_si)
        thresholds_si = thresholds_si.reshape((n,n))*1e3
    
    # A function to measure of the threshold BIO model in CC, for different AIS start positions and Gna
    print ('BIO model')
    def BIO_model_threshold_in_CC(ais_start, ais_end, gna_tot):
        print (ais_start, gna_tot)
        params = params_all
        pulse_length = 50.*ms
        
        neuron = model_Na_Kv1(params=params, resting_vm =-75.*mV, Na_start = ais_start, Na_end = ais_end, gna_tot=gna_tot, density=False)
        i_rheo = measure_current_threshold(params, neuron, -75.*mV, ais_start, ais_end, pulse_length)
        print ('current threshold:', i_rheo)
        
        neuron = model_Na_Kv1(params=params, resting_vm =-75.*mV, Na_start = ais_start, Na_end = ais_end, gna_tot=gna_tot, density=False)
        v_thres, _,_,_ = measure_voltage_threshold(params, neuron, -75.*mV, ais_start, ais_end, i_rheo, pulse_length) 
    
        return v_thres 
    
    # run the simulations with joblib
    if __name__ == '__main__':
        results_bio = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC)(start, start, gna) for start in starts for gna in gnas)
        thresholds_bio = array(results_bio)
        thresholds_bio = thresholds_bio.reshape((n,n))*1e3
    
    # Save the data in an npz file
    savez('figure_7', starts/um, gnas/nS, thresholds_si/mV, thresholds_bio/mV)
    
### THEORY
    
thresholds_pred  = zeros((n,n)) # withotu the contribution of inactivation
thresholds_pred_inact  = zeros((n,n)) # with the contribution of inactivation

for i in range(n):
    for j in range(n):
        v_thres_pred = Va - Ka - Ka * log(starts[i] * ra * gnas[j] * (ENa-Va)/Ka)
        v_thres_pred_inact = Va - Ka - Ka * log(starts[i] * ra * gnas[j] * (ENa-Va)/Ka) + Ka*log(1.+exp((-75.*mV-Vh)/Kh))
        thresholds_pred[i,j] = v_thres_pred/mV
        thresholds_pred_inact[i,j] = v_thres_pred_inact/mV

# Slopes
slopes_x_si = zeros(n)
slopes_gna_si = zeros(n)
slopes_x_bio = zeros(n)
slopes_gna_bio = zeros(n)
for k in range(n):
    sl_x_si, _, _, _, _ = stats.linregress(-log(ra*starts*gnas[k]), thresholds_si[:,k])
    sl_gna_si, _, _, _, _ = stats.linregress(-log(ra*starts[k]*gnas), thresholds_si[k,:])
    slopes_x_si[k] = sl_x_si
    slopes_gna_si[k] = sl_gna_si
    sl_x_bio, _, _, _, _ = stats.linregress(-log(ra*starts*gnas[k]), thresholds_bio[:,k])
    sl_gna_bio, _, _, _, _ = stats.linregress(-log(ra*starts[k]*gnas), thresholds_bio[k,:])
    slopes_x_bio[k] = sl_x_bio
    slopes_gna_bio[k] = sl_gna_bio

print ('Slopes:')
print ('SI model VC:', 'delta:', slopes_x_si, 'G:', slopes_gna_si)
print ('BIO model CC:', 'delta:', slopes_x_bio, 'G:', slopes_gna_bio)

### Plots

# Custom colormap
top = get_cmap('Oranges_r', 128)
bottom = get_cmap('Blues', 128)
newcolors = vstack((top(linspace(0.5, 1, 128)),
                       bottom(linspace(0.5, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')

# Figure
cmap = get_cmap(newcmp)
colors = [cmap(i) for i in np.linspace(0, 1, n)]
gnas_label = gnas/nS
starts_label = starts/um

fig = figure(2, figsize=(6,8))

# Panel A: theory: threshold vs AIS start position
ax1 = subplot(321)
ax1.text(33, -47,'Theory', fontsize=12)

for k in range(n):
    semilogx(starts/um, thresholds_pred[:,k],  color=colors[k], label='gna = %i nS' %gnas_label[k])
ax1.set_ylabel('$V_s$ (mV)') #, fontsize=14)
ax1.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
ax1.set_ylim(-75,-50)
ax1.text(10,-57, '$G$ = %i nS' %gnas_label[0], fontsize=10, color=colors[0]) 
ax1.text(17.5,-74,'$G$ = %i nS'%gnas_label[n-1], fontsize=10, color=colors[n-1])

ax1.text(6., -52,'A', fontsize=14, weight='bold')

# Panel B: theory: threshold vs total Na conductance in the AIS

ax12 = subplot(322, sharey = ax1)
for k in range(n):
    semilogx(gnas/nS, thresholds_pred[k,:], color=colors[k], label='x = %i um' %starts_label[k])
ax12.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    left=False,
    labelbottom=False,
    labelleft=False) # labels along the bottom edge are off
ax12.set_ylim(-75,-50)
ax12.text(200,-57, '$\Delta$ = %i $\mu$m' %starts_label[0], fontsize=10, color=colors[0])
ax12.text(350,-74,'$\Delta$ = %i $\mu$m' %starts_label[n-1], fontsize=10, color=colors[n-1])
ax12.text(130, -52,'B', fontsize=14, weight='bold')

# Panel C: SI model: threshold vs AIS start position

ax2 = subplot(323)
ax2.text(28., -47,'Simple model', fontsize=12)

for k in range(n):
    semilogx(starts/um, thresholds_si[:,k], color=colors[k], label='k=%.2f mV' %slopes_x_si[k])
ax2.set_ylabel('$V_s$ (mV)') 
ax2.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
ax2.set_ylim(-75,-50)

ax2.text(6., -52,'C', fontsize=14, weight='bold')

# Panel D: SI model: threshold vs total Na conductance in the AIS

ax22 = subplot(324, sharey = ax2)
#for k in range(int(n/2)+1):
for k in range(n):
    semilogx(gnas/nS, thresholds_si[k,:], color =colors[k], label='k=%.2f mV' %slopes_gna_si[k])
ax22.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    left=False,
    labelbottom=False,
    labelleft=False) # labels along the bottom edge are off
ax22.set_ylim(-75,-50)

ax22.text(130, -52,'D', fontsize=14, weight='bold')

# Panel E: BIO model: threshold vs AIS start position

ax3 = subplot(325)
ax3.text(24.5   , -47,'Biophysical model', fontsize=12)

for k in range(n):
    semilogx(starts/um, thresholds_bio[:,k], color=colors[k], label='k=%.2f mV' %slopes_x_bio[k])
ax3.set_xlabel('$\Delta$ ($\mu$m)') 
ax3.set_ylabel('$V_s$ (mV)') 
ax3.set_ylim(-75,-50)
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax3.xaxis.set_minor_formatter(ScalarFormatter())

ax3.text(6., -52,'E', fontsize=14, weight='bold')

# Panel F: BIO model: threshold vs total Na conductance in the AIS

ax32 = subplot(326, sharey = ax3)
for k in range(n):
    semilogx(gnas/nS, thresholds_bio[k,:], color =colors[k], label='k=%.2f mV' %slopes_gna_bio[k])
ax32.set_xlabel('$G$ (nS)') #, fontsize=14)
ax32.set_ylim(-75,-50)
ax32.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax32.xaxis.set_minor_formatter(ScalarFormatter())
ax32.text(130, -52,'F', fontsize=14, weight='bold')

subplots_adjust(hspace=0.35, wspace=0.4)

show()    






