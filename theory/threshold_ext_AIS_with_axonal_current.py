"""
A script to solve the equation for the threshold in an extended AIS that starts at the soma,
with a current injected at the AIS end.
"""

from brian2 import *
from scipy.optimize import fsolve
from shared.models.params_all import *

ra = 4*Ri/(pi*axon_diam**2)

# Equation to be solved for y (attention: I/2y must be smaller than one)
func = lambda y, I : y*tanh(y) - (I/(2.*y)) * tanh(y) * ( 1./(1.-(I/(2.*y))**2) + 1.) - (1./((2.*y/I)**2 - 1)) - 1 + I/2

## The threshold is a function of y
nu0_func = lambda y, I: 2.*log(sqrt(2.)*y) - 2.*log(cosh( y + arctanh(I/(2.*y)) ) ) 
nuL_func = lambda y, I: 2.*log(sqrt(2.)*y) + log(1.- (I/(2.*y))**2 ) 
V0_func = lambda y, I, L: Va + Ka*nu0_func(y, I) - Ka*log(L*ra*Gna*(ENa-Va)/Ka)
VL_func = lambda y, I, L: Va + Ka*nuL_func(y, I) - Ka*log(L*ra*Gna*(ENa-Va)/Ka)
    
# We solve for y to calculate the threshold
y_initial_guess = y0

n = 101
currents = linspace(-0.5, 0., n)*nA
ais_lengths = linspace(1, 40., n)*um

y_solutions = zeros((n,n))
nu0s_num = zeros((n,n))
nuLs_num = zeros((n,n))
thresholds_num = zeros((n,n))
thresholds_num_ais = zeros((n,n))

for i in range(n):
    curr = currents[i]
    y_lim = y_initial_guess # when I=0, y=y0 for all L and when L=0, y is also close to y0 for all I
    for j in range(n):
        length = ais_lengths[j]
        
        true_curr =  (ra*length*curr)/Ka # in original units
        
        if abs(curr) < 2*y_lim*Ka/(ra*length): 
            y_solution = fsolve(func, y_initial_guess, args=true_curr)
            y_solutions[i,j] =  y_solution
            nu0s_num[i,j] = nu0_func(y_solution, true_curr)
            nuLs_num[i,j] = nuL_func(y_solution, true_curr)
                
            v0 = V0_func(y_solution, true_curr , length)
            vL = VL_func(y_solution, true_curr , length)
            thresholds_num[i,j] = v0/mV
            thresholds_num_ais[i,j] = vL/mV
        else:
            y_solutions[i,j] =  nan
            nu0s_num[i,j] = nan
            nuLs_num[i,j] = nan
            thresholds_num[i,j] = nan
            thresholds_num_ais[i,j] = nan
            
        y_lim = y_solution

# Plots

figure('BVP with current', figsize=(11,9))
lengths = ais_lengths/um
currs = currents/nA

subplot(331)
plot(currents/nA, y_solutions[:,n/4], '.', label='L=%.2f' %lengths[n/4])
plot(currents/nA, y_solutions[:, n/2], '.', label='L=%.2f' %lengths[n/2])
plot(currents/nA, y_solutions[:,n-1], '.', label='L=%.2f' %lengths[n-1])
ylabel('y')
xlabel('I')
legend(frameon=False)
#ylim(1,1.5)

subplot(332)
plot(currents/nA, nu0s_num[:,n/4], label='L=%.2f' %lengths[n/4])
plot(currents/nA, nu0s_num[:,n/2], label='L=%.2f' %lengths[n/2])
plot(currents/nA, nu0s_num[:,n-1], label='L=%.2f' %lengths[n-1])
ylabel('v0')
xlabel('I')
legend(frameon=False)
#ylim(-0.5,1.5)

subplot(333)
plot(currents/nA, nuLs_num[:,n/4], label='L=%.2f' %lengths[n/4])
plot(currents/nA, nuLs_num[:,n/2], label='L=%.2f' %lengths[n/2])
plot(currents/nA, nuLs_num[:,n-1], label='L=%.2f' %lengths[n-1])
ylabel('vL')
xlabel('I')
legend(frameon=False)
#ylim(-0.5,1.5)

subplot(334)
plot(ais_lengths/um, thresholds_num[0,:], label='I=%.2f' %currs[0])
plot(ais_lengths/um, thresholds_num[n/2,:], label='I=%.2f' %currs[n/2])
plot(ais_lengths/um, thresholds_num[n-1,:], label='I=%.2f' %currs[n-1])
ylabel('Threshold soma (mV)')
xlabel('AIS length (um)')
legend(frameon=False)
ylim(-70,-40)

subplot(335)
plot(ais_lengths/um, thresholds_num_ais[0,:], label='I=%.2f' %currs[0])
plot(ais_lengths/um, thresholds_num_ais[n/2,:], label='I=%.2f' %currs[n/2])
plot(ais_lengths/um, thresholds_num_ais[n-1,:], label='I=%.2f' %currs[n-1])
ylabel('Threshold AIS end (mV)')
xlabel('AIS length (um)')
legend(frameon=False)
ylim(-70,-40)

subplot(336)
plot(ais_lengths/um, thresholds_num_ais[0,:]-thresholds_num[0,:], label='I=%.2f' %currs[0])
plot(ais_lengths/um, thresholds_num_ais[n/2,:]-thresholds_num[n/2,:], label='I=%.2f' %currs[n/2])
plot(ais_lengths/um, thresholds_num_ais[n-1,:]-thresholds_num[n-1,:], label='I=%.2f' %currs[n-1])
ylabel('$\Delta$V (mV)')
xlabel('AIS length (um)')
legend(frameon=False)

subplot(337)
plot(currents/nA, thresholds_num[:,n/4], label='L=%.2f' %lengths[n/4])
plot(currents/nA, thresholds_num[:,n/2], label='L=%.2f' %lengths[n/2])
plot(currents/nA, thresholds_num[:,n-1], label='L=%.2f' %lengths[n-1])
ylabel('Threshold soma (mV)')
xlabel('I')
legend(frameon=False)
ylim(-70,-40)

subplot(338)
plot(currents/nA, thresholds_num_ais[:,n/4], label='L=%.2f' %lengths[n/4])
plot(currents/nA, thresholds_num_ais[:,n/2], label='L=%.2f' %lengths[n/2])
plot(currents/nA, thresholds_num_ais[:,n-1], label='L=%.2f' %lengths[n-1])
ylabel('Threshold AIS end (mV)')
xlabel('I')
legend(frameon=False)
ylim(-70,-40)

subplot(339)
plot(currents/nA, thresholds_num_ais[:,n/4]-thresholds_num[:,n/4], label='L=%.2f' %lengths[n/4])
plot(currents/nA, thresholds_num_ais[:,n/2]-thresholds_num[:,n/2], label='L=%.2f' %lengths[n/2])
plot(currents/nA, thresholds_num_ais[:,n-1]-thresholds_num[:,n-1], label='L=%.2f' %lengths[n-1])
ylabel('$\Delta$V (mV)')
xlabel('I')
legend(frameon=False)

tight_layout()

show()





























