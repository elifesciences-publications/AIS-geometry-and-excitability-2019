
"""
AIS geometry and excitability, section 3: "A spatially extended AIS" (figure 10)

A figure to summarize the logarithmic dependence of the threshold on the 
Na conductance density, AIS length and middle position.

"""
from __future__ import print_function
from brian2 import *
from shared import params_all, model_Na_Kv1, measure_current_threshold, measure_voltage_threshold
from scipy import stats
from joblib import Parallel, delayed
from mpl_toolkits.axes_grid1 import SubplotDivider, LocatableAxes, Size
from matplotlib.colors import ListedColormap

only_plotting = False # to plot the figure without running the simulations

# Parameters
n = 4
m = 3
p = 4
defaultclock.dt = 0.005*ms
params = params_all
ra = 4.*params.Ri/(pi*params.axon_diam**2)

### SIMULATIONS

if only_plotting: # load data
    data = load('figure_10.npz')
    ais_geom = data['arr_0']*um
    gna_densities = data['arr_1']* (siemens/meter**2)
    thresholds_gna = data['arr_2']
    starts = data['arr_3']*um
    lengths = data['arr_4']*um
    threshold_low = data['arr_5']
    threshold_high = data['arr_6']
    
    # AIS middle position from AIS start and length
    x_mids = zeros((m,p))
    for i in range(m):
        for j in range(p):
            x_mids[i,j] = starts[i]+0.5*lengths[j]
    x_mids_vals = unique(x_mids)*1e6
else: # running simulations
    
    # A function to measure threshold vs Na density for different AIS middle position and length
    def BIO_model_in_CC_threshold_density(ais_start, ais_end, gna_dens):
        defaultclock.dt = 0.005*ms
        params = params_all
        resting_vm = -75.*mV
        pulse_length = 50.*ms
        print (ais_start, ais_end)
            
        # current threshold
        neuron = model_Na_Kv1(params, resting_vm, ais_start, ais_end, density = True, gna_density = gna_dens)
        i_rheo = measure_current_threshold(params, neuron, resting_vm, ais_start, ais_end, pulse_length = pulse_length)  
        
        # voltage threshold
        neuron = model_Na_Kv1(params, resting_vm, ais_start, ais_end, density=True, gna_density = gna_dens)
        vs, va, _, _ = measure_voltage_threshold(params, neuron, resting_vm, ais_start, ais_end, i_rheo = i_rheo, pulse_length = pulse_length) 
        
        return vs
    
    print ('measuring thresholds for varying Na densities' )
    gna_densities = linspace(3000., 6000., n) * (siemens/meter**2)
    ais_geom = [[10.,20.], [0.,40.], [20.,20.],[10.,40.]]*um 
    
    # Measuring the threshold 
    if __name__ == '__main__':
        thresholds_gna = Parallel(n_jobs = 4)(delayed(BIO_model_in_CC_threshold_density)(start, start+length, gna_dens) \
                                 for (start, length) in ais_geom for gna_dens in gna_densities)
    
    thresholds_gna = array(thresholds_gna)
    thresholds_gna = thresholds_gna.reshape((2,2,n))*1e3
    
    # A function to measure threshold vs length and middle position for different Na density 
    def BIO_model_threshold_in_CC(x_star, length, gna_dens):
        defaultclock.dt = 0.005*ms
        params = params_all 
        resting_vm = -75.*mV
        pulse_length = 50.*ms
        
        # AIS start
        start = round(x_star, 2)*um - length/2
        # AIS end
        end = start + length
        
        # Removing impossible geometries: the length can not be bigger than two times the middle position
        if length/2 <= round(x_star, 2)*um:
            
            # current threshold
            neuron = model_Na_Kv1(params, resting_vm, Na_start = start, Na_end = end, density=True, gna_density = gna_dens)
            i_rheo = measure_current_threshold(params, neuron, resting_vm=resting_vm, ais_start=start, ais_end=end, pulse_length = pulse_length)  
            
            # voltage threshold
            neuron = model_Na_Kv1(params, resting_vm, Na_start = start, Na_end = end, density=True, gna_density = gna_dens)
            vs, va, _, _ = measure_voltage_threshold(params, neuron, resting_vm=resting_vm, ais_start=start, ais_end=end, i_rheo = i_rheo, pulse_length = pulse_length) 
            print ('threshold:', vs)
            res = vs
        else: # bottom right triangle
            res = nan
    
        return res
    
    starts = linspace(0., 20., m)*um
    lengths = linspace(10., 40., p)*um
    
    # Middle position from AIS start and end
    x_mids = zeros((m,p))
    for i in range(m):
        for j in range(p):
            x_mids[i,j] = starts[i]+0.5*lengths[j]
    x_mids_vals = unique(x_mids)*1e6
    
    # Measuring the threshold with high Na density
    print ('measuring threshold for varying AIS length and middle position with high Na density')
    if __name__ == '__main__':
        threshold_high = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC)(x_mid, length, 5000. * (siemens/meter**2)) for x_mid in x_mids_vals for length in lengths)
    
    threshold_high = array(threshold_high)
    threshold_high = threshold_high.reshape((len(x_mids_vals),p))*1e3
    
    # Measuring the threshold with low Na density
    print ('measuring threshold for varying AIS length and middle position with low Na density')
    if __name__ == '__main__':
        threshold_low = Parallel(n_jobs = 5)(delayed(BIO_model_threshold_in_CC)(x_mid, length, 3000. * (siemens/meter**2)) for x_mid in x_mids_vals for length in lengths)
    
    threshold_low = array(threshold_low)
    threshold_low = threshold_low.reshape((len(x_mids_vals),p))*1e3
    
    #Save the data
    savez('figure_10', ais_geom/um, gna_densities, thresholds_gna, starts/um, lengths/um, \
                              threshold_low, threshold_high)

### Re-arranging the arrays 

# Threshold vs middle position
# short AIS, low density: 
thresholds_mid_shortlow = threshold_low[:, 1] # L=20
# short AIS, high density
thresholds_mid_shorthigh = threshold_high[:, 1] # L=20
# long AIS, low density
thresholds_mid_longlow = threshold_low[:, 3] # L=30
# long AIS, high density
thresholds_mid_longhigh = threshold_high[:, 3] # L=30

# Threshold vs length
# proximal AIS, high density 
thresholds_len_proxhigh = threshold_high[3, :]  
# distal AIS, high density
thresholds_len_disthigh = threshold_high[5,  :]  
# distal AIS, low density
thresholds_len_distlow = threshold_low[5,  :] 
# proximal AIS, low density 
thresholds_len_proxlow = threshold_low[3, :]

###Plots

# Custom colormap
top = get_cmap('Oranges_r', 128)
bottom = get_cmap('Blues', 128)
newcolors = vstack((top(linspace(0.5, 1, 128)),
                       bottom(linspace(0.4, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')

# Figure
cmap = plt.get_cmap(newcmp)
colors = [cmap(i) for i in np.linspace(0, 1, 4)]

fig = figure(3, figsize=(9,3))

# Panel A: threshold vs Na conductance density for varying AIS start position and length
ax1 = subplot(131)
slope_gna_proxshort, _, _, _, _ = stats.linregress(log(gna_densities/ (siemens/meter**2)), thresholds_gna[0,0,:]) 
slope_gna_distshort, _, _, _, _ = stats.linregress(log(gna_densities/ (siemens/meter**2)), thresholds_gna[1,0,:]) 
slope_gna_proxlong, _, _, _, _ = stats.linregress(log(gna_densities/ (siemens/meter**2)), thresholds_gna[0,1,:]) 
slope_gna_distlong, _, _, _, _ = stats.linregress(log(gna_densities/ (siemens/meter**2)), thresholds_gna[1,1,:]) 
print ('Panel A')
print ('light blue:', 'x*=', ais_geom[0,0] + ais_geom[0,1]/2, 'L=',  ais_geom[0,1], 'ka=', slope_gna_proxshort )
print ('light orange:', 'x*=', ais_geom[1,0] + ais_geom[1,1]/2, 'L=',  ais_geom[1,1], 'ka=', slope_gna_distshort )
print ('dark blue: ', 'x*=', ais_geom[2,0] + ais_geom[2,1]/2, 'L=',  ais_geom[2,1], 'ka=', slope_gna_proxlong)
print ('dark orange: ',  'x*=', ais_geom[3,0] + ais_geom[3,1]/2, 'L=',  ais_geom[3,1],'ka=', slope_gna_distlong)

semilogx(gna_densities/ (siemens/meter**2), thresholds_gna[0,0,:], color=colors[2]) 
semilogx(gna_densities/ (siemens/meter**2), thresholds_gna[1,0,:], color=colors[1]) 
semilogx(gna_densities/ (siemens/meter**2), thresholds_gna[0,1,:], color=colors[3]) 
semilogx(gna_densities/ (siemens/meter**2), thresholds_gna[1,1,:], color=colors[0]) 
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.xaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
ax1.set_yticks([-70,-60,-50,-40])
ax1.set_yticklabels(['-70','-60','-50','-40'])
xlabel('$g$ (S/$m^2$)')
ylabel('$V_s$ (mV)')
ylim(-70,-40)

ax1.text(2300, -40,'A', fontsize=14, weight='bold')

# Panel B: threshold vs AIS middle position for varying AIS length and Na conductance density
ax2 = subplot(132)
slope_mid_shortlow, _, _, _, _ = stats.linregress(log(x_mids_vals[1:]), thresholds_mid_shortlow[1:]) 
slope_mid_shorthigh, _, _, _, _ = stats.linregress(log(x_mids_vals[1:]), thresholds_mid_shorthigh[1:])
slope_mid_longlow, _, _, _, _ = stats.linregress(log(x_mids_vals[3:]), thresholds_mid_longlow[3:])  
slope_mid_longhigh, _, _, _, _ = stats.linregress(log(x_mids_vals[3:]), thresholds_mid_longhigh[3:]) 
print ('Panel B')
print ('light blue: gna = 3000', 'L=', lengths[1], 'ka=', slope_mid_shortlow)
print ('light orange: gna = 3000', 'L=', lengths[3], 'ka=', slope_mid_longlow)
print ('dark blue: gna = 5000', 'L=', lengths[1],'ka=', slope_mid_shorthigh)
print ('dark orange: gna = 5000', 'L=', lengths[3],'ka=', slope_mid_longhigh)

semilogx(x_mids_vals, thresholds_mid_shortlow, color=colors[2])
semilogx(x_mids_vals, thresholds_mid_longlow, color=colors[1])
semilogx(x_mids_vals, thresholds_mid_shorthigh, color=colors[3])
semilogx(x_mids_vals, thresholds_mid_longhigh, color=colors[0])
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.xaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
ax2.set_yticks([-70,-60,-50,-40])
ax2.set_yticklabels(['-70','-60','-50','-40'])
xlabel('$x_{1/2}$ ($\mu$m)', )
ylim(-70,-40)

ax2.text(6.25, -40,'B', fontsize=14, weight='bold')

# Panel C: threshold vs AIS length for varying AIS start position and Na conductance density

ax3 = subplot(133)
slope_len_distlow, _, _, _, _ = stats.linregress(log(lengths/um), thresholds_len_distlow) 
slope_len_proxhigh, _, _, _, _ = stats.linregress(log(lengths/um), thresholds_len_proxhigh)
slope_len_disthigh, _, _, _, _ = stats.linregress(log(lengths/um), thresholds_len_disthigh)  
slope_len_proxlow, _, _, _, _ = stats.linregress(log(lengths/um), thresholds_len_proxlow)
print ('Panel C')
print ('light blue: gna = 3000', 'x=', x_mids_vals[3], 'ka=', slope_len_proxlow)
print ('light orange: gna = 5000', 'x=', x_mids_vals[3],'ka=', slope_len_proxhigh)
print ('dark blue: gna = 3000', 'x=', x_mids_vals[5], 'ka=', slope_len_distlow)
print ('dark orange: gna = 5000', 'x=', x_mids_vals[5], 'ka=', slope_len_disthigh)

semilogx(lengths/um, thresholds_len_proxhigh, color=colors[1])
semilogx(lengths/um, thresholds_len_proxlow, color=colors[2])
semilogx(lengths/um, thresholds_len_distlow, color=colors[3])
semilogx(lengths/um, thresholds_len_disthigh, color=colors[0])
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax3.xaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
ax3.set_yticks([-70,-60,-50,-40])
ax3.set_yticklabels(['-70','-60','-50','-40'])
xlabel('$L$ ($\mu$m)')
ylim(-70,-40)

ax3.text(6.25, -40,'C', fontsize=14, weight='bold')

subplots_adjust(top=0.88, bottom=0.29, wspace=0.3)

show()

print ('Distal shift of a 40 um long AIS with high g (x* from 25 to 30::', thresholds_mid_longhigh[5] - thresholds_mid_longhigh[4])












