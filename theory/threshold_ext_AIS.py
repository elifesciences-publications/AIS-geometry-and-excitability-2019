"""
A script to solve the equation for the threshold in an extended AIS that starts away from soma.
"""

from brian2 import *
from scipy.optimize import fsolve
from shared.models.params_all import *
from mpl_toolkits.axes_grid1 import SubplotDivider, LocatableAxes, Size

ra = 4*Ri/(pi*axon_diam**2)

# Equation to be solved for y 
func = lambda y, delta : (1+delta)*y*tanh(y) + delta*y**2 * (1-(tanh(y))**2) - 1

# The threshold is a function of y
nu0_func = lambda y, delta: 2*log(sqrt(2)*y) - 2*log(cosh(y)) - 2*delta*y*tanh(y) 
nuL_func = lambda y: 2*log(sqrt(2)*y) 
V0_func = lambda y, delta, L: Va + Ka*nu0_func(y, delta) - Ka*log(L*ra*Gna*(ENa-Va)/Ka)
VL_func = lambda y, L: Va + Ka*nuL_func(y) - Ka*log(L*ra*Gna*(ENa-Va)/Ka)
corrective_term = lambda y, delta: nu0_func(y, delta) + log(delta+0.5) + 1

# Analytical approximation for the threshold (to compare)
def nu0_prediction(start, length):
    end = start + length
    return 2.*log(sqrt(2.)*y0*(length/end))-2.*log(cosh(y0*length/end))-2.*y0*(start/end)*tanh(y0*length/end)

def V0_prediction(start, length):
    vth = Va + Ka*nu0_prediction(start, length) - Ka*log(length*ra*Gna*(ENa-Va)/Ka)
    return vth
    
def corrective_term_prediction(start, length):
    return nu0_prediction(start, length) + log(start/length + 0.5) + 1.
    
def V0_point_ais_prediction(start):
    vth = Va - Ka - Ka*log(start*ra*Gna*(ENa-Va)/Ka)
    return vth
    
# We solve for y to calculate the threshold
y_initial_guess = 1.2

n = 51
ais_starts = linspace(0., 20., n)*um
ais_lengths = linspace(.005, 40., n)*um

thresholds_num = zeros((n,n))
thresholds_num_ais = zeros((n,n))
thresholds_pred = zeros((n,n))
thresholds_pt_ais = zeros((n,n))

y_solutions = zeros((n,n))
y_approxs = zeros((n,n))

nu0s_num = zeros((n,n))
nuLs_num = zeros((n,n))
nu0s_pred = zeros((n,n))

corrections_num = zeros((n,n))
corrections_pred = zeros((n,n))

x_stars = zeros((n,n))

for i in range(n):
    start = ais_starts[i]
    for j in range(n):
        length = ais_lengths[j]
        delta = start/length
        
        y_solution = fsolve(func, y_initial_guess, args=delta)
        y_approx = y0*(1./(1.+delta))        
        y_solutions[i,j] =  y_solution
        y_approxs[i,j] = y_approx
        
        nu0s_num[i,j] = nu0_func(y_solution, delta)
        nuLs_num[i,j] = nuL_func(y_solution)
        nu0s_pred[i,j] = nu0_prediction(start, length)
        
        v0 = V0_func(y_solution, delta, length)
        vL = VL_func(y_solution, length)
        thresholds_num[i,j] = v0/mV
        thresholds_num_ais[i,j] = vL/mV
        thresholds_pred[i,j] = V0_prediction(start, length)/mV
        
        corr_num = corrective_term(y_solution, delta)
        corrections_num[i,j] = corr_num
        corrections_pred[i,j] = corrective_term_prediction(start, length)

lens = ais_lengths/um
starts = ais_starts/um

# Plots

figure(1, figsize=(11,9))

subplot(331)
plot(ais_starts/um, y_solutions[:,0], 'g', label='L$\simeq$0')
plot(ais_starts/um, y_approxs[:,0], '--g', label='approx')
plot(ais_starts/um, y_solutions[:,int(n/2)], 'b', label='L=%.1f' %lens[int(n/2)])
plot(ais_starts/um, y_approxs[:,int(n/2)], '--b')
plot(ais_starts/um, y_solutions[:,n-1], 'r', label='L=%.1f' %lens[n-1])
plot(ais_starts/um, y_approxs[:,n-1], '--r')
ylabel('$y^* (\Delta, L)$')
xlabel('$\Delta$ (um)')
legend(frameon=False)

subplot(332)
plot(ais_starts/um, nu0s_num[:,0], 'g', label='L$\simeq$0')
#plot(ais_starts/um, nu0s_pred[:,0], '--g', label='approx')
plot(ais_starts/um, nu0s_num[:,int(n/2)], 'b', label='L=%.1f' %lens[int(n/2)])
#plot(ais_starts/um, nu0s_pred[:,int(n/2)], '--b')
plot(ais_starts/um, nu0s_num[:,n-1], 'r', label='L=%.1f' %lens[n-1])
#plot(ais_starts/um, nu0s_pred[:,n-1], '--r')
ylabel('$ v_0(\Delta, L)$')
xlabel('$\Delta$ (um)')
ylim(-10,2)
legend(frameon=False)

subplot(333)
plot(ais_starts/um, nuLs_num[:,0], 'g', label='L$\simeq$0')
#plot(ais_starts/um, nu0s_pred[:,0], '--g', label='approx')
plot(ais_starts/um, nuLs_num[:,int(n/2)], 'b', label='L=%.1f' %lens[int(n/2)])
#plot(ais_starts/um, nu0s_pred[:,int(n/2)], '--b')
plot(ais_starts/um, nuLs_num[:,n-1], 'r', label='L=%.1f' %lens[n-1])
#plot(ais_starts/um, nu0s_pred[:,n-1], '--r')
ylabel('$ v_L(\Delta, L)$')
xlabel('$\Delta$ (um)')
ylim(-10,2)
legend(frameon=False)                          

subplot(334)
plot(ais_starts/um, thresholds_num[:,0], 'g', label='L$\simeq$ 0')
#plot(ais_starts/um, thresholds_pred[:,0], '--g', labe$l='approximation')
#plot(ais_starts/um, thresholds_pt_ais[:,0], '.g', label='equivalent point AIS in $\Delta$+L/2')
plot(ais_starts/um, thresholds_num[:,int(n/2)], 'b', label='L=%.2f' %lens[int(n/2)])
#plot(ais_starts/um, thresholds_pred[:,int(n/2)], '--b')
#plot(ais_starts/um, thresholds_pt_ais[:,int(n/2)], '.b')
plot(ais_starts/um, thresholds_num[:,n-1], 'r', label='L=%.2f' %lens[n-1])
#plot(ais_starts/um, thresholds_pred[:,n-1], '--r')
#plot(ais_starts/um, V0_point_ais_prediction(ais_starts+20*um)/mV, '--k', label='prediction point AIS')
#plot(ais_starts/um, thresholds_pt_ais[:,n-1], '.r')
ylabel('Threshold (mV)')
xlabel('$\Delta$ (um)')
ylim(-70,-45)
legend(frameon=False)

subplot(335)
plot(ais_starts/um, thresholds_num_ais[:,0], 'g', label='L$\simeq$ 0')
#plot(ais_starts/um, thresholds_pred[:,0], '--g', labe$l='approximation')
#plot(ais_starts/um, thresholds_pt_ais[:,0], '.g', label='equivalent point AIS in $\Delta$+L/2')
plot(ais_starts/um, thresholds_num_ais[:,int(n/2)], 'b', label='L=%.2f' %lens[int(n/2)])
#plot(ais_starts/um, thresholds_pred[:,int(n/2)], '--b')
#plot(ais_starts/um, thresholds_pt_ais[:,int(n/2)], '.b')
plot(ais_starts/um, thresholds_num_ais[:,n-1], 'r', label='L=%.2f' %lens[n-1])
#plot(ais_starts/um, thresholds_pred[:,n-1], '--r')
#plot(ais_starts/um, V0_point_ais_prediction(ais_starts+20*um)/mV, '--k', label='prediction point AIS')
#plot(ais_starts/um, thresholds_pt_ais[:,n-1], '.r')
ylabel('Threshold AIS (mV)')
xlabel('$\Delta$ (um)')
ylim(-70,-45)
legend(frameon=False)

subplot(336)
plot(ais_starts/um, thresholds_num_ais[:,0]-thresholds_num[:,0], 'g', label='L$\simeq$ 0')
#plot(ais_starts/um, thresholds_pred[:,0], '--g', labe$l='approximation')
#plot(ais_starts/um, thresholds_pt_ais[:,0], '.g', label='equivalent point AIS in $\Delta$+L/2')
plot(ais_starts/um, thresholds_num_ais[:,int(n/2)]-thresholds_num[:,int(n/2)], 'b', label='L=%.2f' %lens[int(n/2)])
#plot(ais_starts/um, thresholds_pred[:,int(n/2)], '--b')
#plot(ais_starts/um, thresholds_pt_ais[:,int(n/2)], '.b')
plot(ais_starts/um, thresholds_num_ais[:,n-1]-thresholds_num[:,n-1], 'r', label='L=%.2f' %lens[n-1])
#plot(ais_starts/um, thresholds_pred[:,n-1], '--r')
#plot(ais_starts/um, V0_point_ais_prediction(ais_starts+20*um)/mV, '--k', label='prediction point AIS')
#plot(ais_starts/um, thresholds_pt_ais[:,n-1], '.r')
ylabel('$\Delta$V (mV)')
xlabel('$\Delta$ (um)')
ylim(4.5,6)
legend(frameon=False)

subplot(337)
plot(ais_lengths/um, thresholds_num[0,:], 'g', label='$\Delta$=%.2f' %starts[0])
#plot(ais_starts/um, thresholds_pred[:,0], '--g', labe$l='approximation')
#plot(ais_starts/um, thresholds_pt_ais[:,0], '.g', label='equivalent point AIS in $\Delta$+L/2')
plot(ais_lengths/um, thresholds_num[int(n/2),:], 'b', label='$\Delta$=%.2f' %starts[int(n/2)])
#plot(ais_starts/um, thresholds_pred[:,int(n/2)], '--b')
#plot(ais_starts/um, thresholds_pt_ais[:,int(n/2)], '.b')
plot(ais_lengths/um, thresholds_num[n-1,:], 'r', label='$\Delta$=%.2f' %starts[n-1])
#plot(ais_starts/um, thresholds_pred[:,n-1], '--r')
#plot(ais_starts/um, V0_point_ais_prediction(ais_starts+20*um)/mV, '--k', label='prediction point AIS')
#plot(ais_starts/um, thresholds_pt_ais[:,n-1], '.r')
ylabel('Threshold (mV)')
xlabel('L (um)')
ylim(-70,-45)
legend(frameon=False)

subplot(338)
plot(ais_lengths/um, thresholds_num_ais[0,:], 'g', label='$\Delta$=%.2f' %starts[0])
#plot(ais_starts/um, thresholds_pred[:,0], '--g', labe$l='approximation')
#plot(ais_starts/um, thresholds_pt_ais[:,0], '.g', label='equivalent point AIS in $\Delta$+L/2')
plot(ais_lengths/um, thresholds_num_ais[int(n/2),:], 'b', label='$\Delta$=%.2f' %starts[int(n/2)])
#plot(ais_starts/um, thresholds_pred[:,int(n/2)], '--b')
#plot(ais_starts/um, thresholds_pt_ais[:,int(n/2)], '.b')
plot(ais_lengths/um, thresholds_num_ais[n-1,:], 'r', label='$\Delta$=%.2f' %starts[n-1])
#plot(ais_starts/um, thresholds_pred[:,n-1], '--r')
#plot(ais_starts/um, V0_point_ais_prediction(ais_starts+20*um)/mV, '--k', label='prediction point AIS')
#plot(ais_starts/um, thresholds_pt_ais[:,n-1], '.r')
ylabel('Threshold AIS (mV)')
xlabel('L (um)')
ylim(-70,-45)
legend(frameon=False)

subplot(339)
plot(ais_lengths/um, thresholds_num_ais[0,:]-thresholds_num[0,:], 'g', label='$\Delta$=%.2f' %starts[0])
#plot(ais_starts/um, thresholds_pred[:,0], '--g', labe$l='approximation')
#plot(ais_starts/um, thresholds_pt_ais[:,0], '.g', label='equivalent point AIS in $\Delta$+L/2')
plot(ais_lengths/um, thresholds_num_ais[int(n/2),:]-thresholds_num[int(n/2),:], 'b', label='$\Delta$=%.2f' %starts[int(n/2)])
#plot(ais_starts/um, thresholds_pred[:,int(n/2)], '--b')
#plot(ais_starts/um, thresholds_pt_ais[:,int(n/2)], '.b')
plot(ais_lengths/um, thresholds_num_ais[n-1,:]-thresholds_num[n-1,:], 'r', label='$\Delta$=%.2f' %starts[n-1])
#plot(ais_starts/um, thresholds_pred[:,n-1], '--r')
#plot(ais_starts/um, V0_point_ais_prediction(ais_starts+20*um)/mV, '--k', label='prediction point AIS')
#plot(ais_starts/um, thresholds_pt_ais[:,n-1], '.r')
ylabel('$\Delta$V (mV)')
xlabel('L (um)')
ylim(4.5,6)
legend(frameon=False)

tight_layout()

show()

