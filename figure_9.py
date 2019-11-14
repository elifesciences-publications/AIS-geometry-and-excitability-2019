
"""

AIS geometry and excitability, section 3: "A spatially extended AIS" (figure 9)

Effect of the middle position of the AIS on the somatic threshold for spike initiation. 
The total Na conductance is constant in the AIS and the middle AIS position and length are varied.

"""
from __future__ import print_function
from brian2 import *
from joblib import Parallel, delayed
from scipy.optimize import fsolve
from mpl_toolkits.axes_grid1 import SubplotDivider, LocatableAxes, Size
from matplotlib.colors import ListedColormap
from shared.models import model_Na_Kv1, params_model_description
from shared.analysis import measure_current_threshold, measure_voltage_threshold

only_plotting = True # to plot the figure without running the simulations

# Parameters
defaultclock.dt = 0.005*ms
params = params_model_description #params_all
ra = 4*params.Ri/(pi*params.axon_diam**2)
n= 5
m= 9

### Divers functions
# Theory: threshold formulas
# Point AIS
def Threshold_point_AIS(x):
    vth = params.Va - params.Ka - params.Ka*log(x*ra*params.Gna*(params.ENa-params.Va)/params.Ka)
    return vth 
# Extended AIS that starts at the soma
def Threshold_AIS_starts_at_soma(L):
    vth = params.Va - 0.12959*params.Ka - params.Ka*log(L*ra*params.Gna*(params.ENa-params.Va)/params.Ka)
    return vth

### SIMULATIONS

if only_plotting: # loading data
    data = load('figure_9.npz')
    starts = data['arr_0']*um
    lengths = data['arr_1']*um
    x_mids_vals = data['arr_2']
    thresholds_pred = data['arr_3']
    thresholds_bio = data['arr_4']

else: # runnning the simulations
    
    # AIS start position and length
    starts = linspace(0., 20., n)*um # starts and length need to have the same interval between values
    lengths = linspace(0., 40., m)*um
    
    # AIS middle positions for each combination for start position and length
    x_mids = zeros((n,m))
    for i in range(n):
        for j in range(m):
            x_mids[i,j] = round(starts[i]/um +0.5*lengths[j]/um , 2) 

    # An array that contains the unique AIS middle positions
    x_mids_vals = unique(x_mids)

    # Neuron morphology: we need to more spatial discretization in the axon
    dend_diam = 6.*um 
    dend_length = 1000.*um 
    axon_diam = 1. * um
    axon_length = 500. * um
    dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=500)
    axon = Cylinder(diameter=axon_diam, length=axon_length, n=1000)
    morpho = Soma(30.*um)
    morpho.dendrite = dendrite
    morpho.axon = axon
    
    # A function to measure the trheshold in the biophysical model, for various AIS middle position and length
    def BIO_model_threshold_in_CC(x_star, length, gna_tot):
        defaultclock.dt = 0.005*ms
        resting_vm = -75.*mV
        pulse_length = 50.*ms
        
        # AIS start
        start = round(x_star, 2)*um - length/2
        # AIS end
        end = start + length
        
        # Removing impossible geometries: the length can not be bigger than two times the middle position
        if length/2 <= round(x_star, 2)*um:
            print ('x1/2:', round(x_star, 2)*um, 'AIS length:', length, 'AIS start:', start, 'AIS end:', end)
            
            # current threshold
            neuron = model_Na_Kv1(params, resting_vm, Na_start = start, Na_end = end, density=False, gna_tot=gna_tot, morpho=morpho)
            i_rheo = measure_current_threshold(params, neuron, resting_vm=resting_vm, ais_start=start, ais_end=end, pulse_length = pulse_length)  
            print ('Rheobase:', i_rheo)
            
            # voltage threshold
            neuron = model_Na_Kv1(params, resting_vm, Na_start = start, Na_end = end, density=False, gna_tot=gna_tot, morpho=morpho)
            vs, va, _, _ = measure_voltage_threshold(params, neuron, resting_vm=resting_vm, ais_start=start, ais_end=end, i_rheo = i_rheo, pulse_length = pulse_length) 
            print ('Threshold:', vs)
            res = vs
        else: # bottom right triangle
            res = nan
    
        return res
    
    if __name__ == '__main__':
        results = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC)(x_mid, length, params.Gna) for x_mid in x_mids_vals for length in lengths)
    
    thresholds_bio = array(results)
    thresholds_bio = thresholds_bio.reshape((len(x_mids_vals),m)) * 1e3
        
    ### THEORY
    
    # Extended AIS that starts away from the soma - exact solution
    # Equation to be solved
    func = lambda y, delta : (1+delta)*y*tanh(y) + delta*y**2 * (1-(tanh(y))**2) - 1
    nu0_func = lambda y, delta: 2*log(sqrt(2)*y) - 2*log(cosh(y)) - 2*delta*y*tanh(y) 
    V0_func = lambda y, delta, L: params.Va + params.Ka*nu0_func(y, delta) - params.Ka*log(L*ra*params.Gna*(params.ENa-params.Va)/params.Ka)
    
    # We solve V0_func to calculate the threshold, for each combination of(x*, L) provided that x* > L/2
    y_initial_guess = 1.2 # it corresponds to the AIS that starts at the soma
    thresholds_pred = zeros((len(x_mids_vals), m))
    
    for i in range(len(x_mids_vals)):
        for j in range(m):
            if j <= i: # unpossible geometries
                if j == 0: # for a point AIS, we need to use the corresponding threshold formula 
                    start = x_mids_vals[i]*um - lengths[j]/2
                    thresholds_pred[i,j] = Threshold_point_AIS(start)/mV
                else:
                    start = x_mids_vals[i]*um - lengths[j]/2
                    #print i, j, x_mids_vals[i]*um, lengths[j], start
                    delta = start/lengths[j]
                    y_solution = fsolve(func, y_initial_guess, args=delta)  
                    v0 = V0_func(y_solution, delta, lengths[j])
                    thresholds_pred[i,j] = v0/mV
            else:
                thresholds_pred[i,j] = nan
    
    #Save the data
    savez('figure_9', starts/um, lengths/um, x_mids_vals, thresholds_pred, thresholds_bio)
    
# Bifurcation condition: we search for the AIS position below which there is no bifuraction, for a fixed Gna
print ('THEORY')
x_pt = linspace(0., 50., 1000)*um

for i in range(len(x_pt)):
    RaGNa = ra*x_pt[i] * params.Gna
    if RaGNa > 1:
        print ('Minimal x1/2 with bifurcation:', x_pt[i])
        print ('Maximal x1/2 without bifurcation:', x_pt[i-1])
        x_mid_min = x_pt[i]
        break
            
### Figure paper

# Custom colormap
top = get_cmap('Oranges_r', 128)
bottom = get_cmap('Blues', 128)
newcolors = vstack((top(linspace(0, 1, 128)),
                   bottom(linspace(0, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')        

fig = figure('fig paper', figsize=(6,6))

# Panel A: equivalent point AIS
l = 501
x = linspace(0.,40., l)*um
x_min = where(x>=x_mid_min)[0][0]

thresholds_pt = zeros(l)
thresholds_ext = zeros(l)
thresholds_ext_mid = zeros(l)
for i in range(l):
    thresholds_pt[i] = Threshold_point_AIS(x[i])
    thresholds_ext[i] = Threshold_AIS_starts_at_soma(x[i]) 
    thresholds_ext_mid[i] = Threshold_point_AIS(x[i]/2)
print ('Threshold difference between point AIS at L and extended AIS that starts at the soma and ends at distance L:', (thresholds_ext[i] - thresholds_pt[i])*1e3, 'mV')
print ('Threshold difference between extended AIS that starts at the soma and ends at distance L and its equivalent point AIS in L/2:', (thresholds_ext[i] - thresholds_ext_mid[i])*1e3, 'mV')

ax1 = subplot2grid((3, 2), (0, 0), colspan=2)
ax1.plot(x[2*x_min:]/um, thresholds_pt[2*x_min:]/mV, 'k', label='Point AIS in L')
ax1.plot(x/um, thresholds_ext/mV, '--k', label='AIS from the soma to L')
ax1.set_xticks([0,10,20,30,40])
ax1.set_xticklabels(['0','10','20','30','40'])
ax1.set_ylabel('$V_s$ (mV)')
ax1.set_xlabel('$L$ ($\mu$m)')
ax1.set_xlim(0, 40)
ax1.set_ylim(-75,-35)

ax2 = ax1.twiny()
ax2.plot(x[2*x_min:]/2/um, thresholds_ext_mid[2*x_min:]/mV, '-', color='forestgreen')
ax2.set_xticks([0,5,10,15,20])
ax2.set_xticklabels(['0','5','10','15','20'])
ax2.set_xlabel('$x_{1/2}$ ($\mu$m)', color='forestgreen')
ax2.set_xlim(0,20)
ax2.tick_params(colors='forestgreen')

ax1.text(-5.25, -35,'A', fontsize=14, weight='bold')

# Colorplots
data_pred = thresholds_pred
data_bio = thresholds_bio

n_pred = len(data_pred)
levels_pred = arange(int(nanmin(data_pred[n_pred-1]))-1,int(nanmax(data_pred[1]))+1, 2)
n_bio = len(data_bio)
levels_bio = arange(int(nanmin(data_bio[n_pred-1]))-1,int(nanmax(data_bio[1]))+1, 2)

cmap_pred= plt.cm.get_cmap(newcmp, levels_pred+1)
cmap_bio= plt.cm.get_cmap(newcmp, levels_bio+1)

X, Y = meshgrid(lengths/um, x_mids_vals)
Z1 = data_pred
Z2 = data_bio

# Panel B: theory: threshold vs AIS middle position and length
    
ax3 = subplot2grid((3, 2), (1, 0), rowspan=2)
ax3.contourf(X, Y, Z1, levels_pred, cmap = cmap_pred)
cont3 = ax3.contour(X,Y,Z1, levels_pred, linewidths=1,  colors='k', linestyles='solid')
ax3.clabel(cont3, inline=1, fmt='%1.0f')
ax3.set_yticks([0,10,20,30,40])
ax3.set_yticklabels(['0','10','20','30','40'])
ax3.set_xlabel('$L$ ($\mu$m)')
ax3.set_ylabel('$x_{1/2}$ ($\mu$m)')
ax3.set_title('Theory')
ax3.set_aspect('equal')

# Region without bifurcation
xx = linspace(0.,10, 100)
y1 = x_mid_min/um * ones(len(xx))
y2 = xx/2
xx_max = where(y2>y1)[0][0]
ax3.fill_between(xx[:xx_max], y1[:xx_max], y2[:xx_max], color='lightgray')


ax3.text(-12, 40,'B', fontsize=14, weight='bold')

# Panel C: BIO model: threshold vs AIS middle position and length
ax4 = subplot2grid((3, 2), (1, 1), rowspan=2)
ax4.contourf(X, Y, Z2, levels_bio, cmap = cmap_bio)
cont4 = ax4.contour(X,Y,Z2, levels_bio, linewidths=1, colors='k', linestyles='solid')
ax4.clabel(cont4, inline=1, fmt='%1.0f')
ax4.set_yticks([0,10,20,30,40])
ax4.set_yticklabels(['0','10','20','30','40'])
ax4.set_xlabel('$L$ ($\mu$m)')
ax4.set_title('Biophysical model')
ax4.set_aspect('equal')

# Region without bifurcation
ax4.fill_between(xx[:xx_max], y1[:xx_max], y2[:xx_max], color='w')
ax4.fill_between(xx[:xx_max], y1[:xx_max], y2[:xx_max], color='lightgray')


ax4.text(-8, 40,'C', fontsize=14, weight='bold')

tight_layout()

show()




