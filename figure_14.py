

"""
AIS geometry and excitability, "Impact of the distal axon".

A script to plot the figure 14.

Extended AIS in a larger neuron.

"""
from __future__ import print_function
from brian2 import *
from shared import params_model_description, model_Na_Kv1, measure_current_threshold, measure_voltage_threshold
from joblib import Parallel, delayed
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.ticker as ticker
from scipy import stats

# Parameters
defaultclock.dt = 0.005*ms
params = params_model_description
GNa = 400.*nS # total Na conductance at the AIS

### Figure 

# Custom colormap
top = get_cmap('Oranges_r', 128)
bottom = get_cmap('Blues', 128)
newcolors = vstack((top(linspace(0.5, 1, 128)),
                       bottom(linspace(0.4, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')
cmap = plt.get_cmap(newcmp)
colors = [cmap(i) for i in np.linspace(0, 1, 4)]

fig = figure('Distal axon', figsize=(6,5))

### PANEL A: large neuron

# loading data
data = load('figure_14a.npz')
starts = data['arr_0']*um
starts_large = data['arr_1']*um
length = data['arr_2']*um
scaling = data['arr_3']
thresholds = data['arr_4']*mV
thresholds_large = data['arr_5']*mV

length_large = sqrt(scaling) * length

x_mids = starts+length/2 # AIS middle position original neuron
x_mids_large = starts_large+length_large/2 # AIS middle position large neuron

# Slopes
slopes, _, _, _, _ = stats.linregress(-log(x_mids/um), thresholds)
slopes_large, _, _, _, _ = stats.linregress(-log(x_mids_large/um), thresholds_large)

print ('Panel A: large neuron')
print ('Slopes:')
print('Threshold vs delta:', slopes, 'mV')
print('Threshold vs delta large neuron:', slopes_large, 'mV')

ax1 = fig.add_subplot(221)

ax1.semilogx(x_mids/um, thresholds, color=colors[3], label='original')#'darkblue')
ax1.semilogx(x_mids_large/um, thresholds_large, color=colors[2], label='large neuron')
ax1.semilogx(x_mids/um, thresholds + 5*log(scaling)*volt, '--', color=colors[3], label='expected shift')
ax1.set_ylabel('$V_s$ (mV)') 
ax1.set_xlabel('$x_{1/2}$ ($\mu$m)' )
ax1.set_ylim(-75,-40)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.xaxis.set_minor_formatter(ScalarFormatter())
ax1.legend(frameon=False, fontsize=8, loc='upper right')

ax1.text(7.5, -40,'A', fontsize=14, weight='bold')

### PANEL B: myelin

# loading data
data = load('figure_14b1.npz')
starts = data['arr_0']*um
length = data['arr_1']*um
thresholds = data['arr_2']*mV
thresholds_myelin = data['arr_3']*mV

x_mids = starts+length/2 # AIS middle position original neuron

# Slopes
slopes, _, _, _, _ = stats.linregress(-log(x_mids/um), thresholds)
slopes_myelin, _, _, _, _ = stats.linregress(-log(x_mids/um), thresholds_myelin)

print ('Panel B: myelin')
print ('Slopes:')
print('Threshold vs delta:', slopes, 'mV')
print('Threshold vs delta myelin:', slopes_myelin, 'mV')

ax2 = fig.add_subplot(222)

ax2.semilogx(x_mids/um, thresholds, color=colors[3], label='original')
ax2.semilogx(x_mids/um, thresholds_myelin,color=colors[2], label='myelin')
ax2.set_ylim(-75,-40)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.xaxis.set_minor_formatter(ScalarFormatter())

### PANEL B: long axon

# loading data
data = load('figure_14b2.npz')
starts = data['arr_0']*um
length = data['arr_1']*um
thresholds = data['arr_2']*mV
thresholds_long = data['arr_3']*mV

x_mids = starts+length/2 # AIS middle position original neuron

# Slopes
slopes, _, _, _, _ = stats.linregress(-log(x_mids/um), thresholds)
slopes_long, _, _, _, _ = stats.linregress(-log(x_mids/um), thresholds_long)

print ('Panel C: long axon')
print ('Slopes:')
print('Threshold vs delta:', slopes, 'mV')
print('Threshold vs delta myelin:', slopes_long, 'mV')

#ax3 = fig.add_subplot(323)

ax2.semilogx(x_mids/um, thresholds_long, '--', color=colors[0], label='long axon')
#ax2.semilogx(x_mids/um, thresholds, 'darkblue')
ax2.set_ylim(-75,-40)
ax2.set_xlabel('$x_{1/2}$ ($\mu$m)' )
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.xaxis.set_minor_formatter(ScalarFormatter())
ax2.legend(frameon=False, fontsize=8)

ax2.text(10, -40,'B', fontsize=14, weight='bold')


### PANEL C: branching after the AIS

# loading data
data = load('figure_14c.npz')
starts = data['arr_0']*um
length = data['arr_1']*um
thresholds = data['arr_2']*mV
thresholds_rall = data['arr_3']*mV
thresholds_emp = data['arr_4']*mV

x_mids = starts+length/2 # AIS middle position original neuron

# Slopes
slopes, _, _, _, _ = stats.linregress(-log(x_mids/um), thresholds)
slopes_rall, _, _, _, _ = stats.linregress(-log(x_mids/um), thresholds_rall)
slopes_emp, _, _, _, _ = stats.linregress(-log(x_mids/um), thresholds_emp)

print ('Panel D: long axon')
print ('Slopes:')
print('Threshold vs delta:', slopes, 'mV')
print('Threshold vs delta Rall:', slopes_rall, 'mV')
print('Threshold vs delta emp:', slopes_emp, 'mV')

ax4 = fig.add_subplot(223)

ax4.semilogx(x_mids/um, thresholds,color=colors[3], label='original')
ax4.semilogx(x_mids/um, thresholds_emp, color=colors[2], label='branched, $p$ = 5/2')
ax4.semilogx(x_mids/um, thresholds_rall, '--', color=colors[0], label='branched, $p$ = 3/2')
ax4.set_ylim(-75,-40)
ax4.set_ylabel('$V_s$ (mV)') 
ax4.set_xlabel('$x_{1/2}$ ($\mu$m)' )
ax4.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax4.xaxis.set_minor_formatter(ScalarFormatter())
ax4.legend(frameon=False, fontsize=8)

ax4.text(9.5, -40,'C', fontsize=14, weight='bold')

### PANEL D: axo-dendritic neuron

# loading data
data = load('figure_14d.npz')
starts = data['arr_0']*um
length = data['arr_1']*um
axial_resistance = data['arr_2']*Mohm
axial_resistance_axo_dend = data['arr_3']*Mohm
thresholds = data['arr_4']*mV
thresholds_axo_dend = data['arr_5']*mV

# Slopes
slopes, _, _, _, _ = stats.linregress(-log(axial_resistance/Mohm), thresholds)
slopes_axo_dend, _, _, _, _ = stats.linregress(-log(axial_resistance_axo_dend/Mohm), thresholds_axo_dend)

print ('Panel E: long axon')
print ('Slopes:')
print('Threshold vs delta:', slopes, 'mV')
print('Threshold vs delta axo-dend:', slopes_axo_dend, 'mV')

ax5 = fig.add_subplot(224)

ax5.semilogx(axial_resistance/Mohm, thresholds,color=colors[3], label='original')
ax5.semilogx(axial_resistance_axo_dend/Mohm, thresholds_axo_dend, color=colors[2], label='axo-dendritic neuron')
ax5.set_ylim(-75,-40)
ax5.set_xlabel('$R_a$ (M$\Omega$)') 
ax5.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax5.xaxis.set_minor_formatter(ScalarFormatter())
ax5.legend(frameon=False, fontsize=8)

ax5.text(8.75, -40,'D', fontsize=14, weight='bold')


tight_layout()

show()   
 






