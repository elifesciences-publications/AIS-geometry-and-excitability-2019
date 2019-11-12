

"""
AIS geometry and excitability, "Methods" (figure 13).

Exact solution of the bifurcation problem. 
"""

from brian2 import *
from scipy.optimize import fsolve
from shared import *
from mpl_toolkits.axes_grid1 import SubplotDivider, LocatableAxes, Size

ra = 4*Ri/(pi*axon_diam**2)

# Equation to be solved for y 
func = lambda y, delta : (1+delta)*y*tanh(y) + delta*y**2 * (1-tanh(y)**2) - 1

# The threshold is a function of y
nu0_func = lambda y, delta: 2*log(sqrt(2)*y) - 2*log(cosh(y)) - 2*delta*y*tanh(y) 
V0_func = lambda y, delta, L: Va + Ka*nu0_func(y, delta) - Ka*log(L*ra*Gna*(ENa-Va)/Ka)
corrective_term = lambda y, delta: nu0_func(y, delta) + log(delta+0.5) + 1

# We solve for y to calculate the threshold
y_initial_guess = 1.2 # corresponding to an extended AIS starting at the soma

n = 101
deltas = linspace(0., 5., n) # AIS start position divided by AIS length

y_solutions = zeros(n) # solutions of func
nu0s_num = zeros(n) # the nu_O term in the trheshold formula
corrections_num = zeros(n) # the corrective term F(delta)

for i in range(n):
    delta = deltas[i]    
    y_solution = fsolve(func, y_initial_guess, args=delta)    
    y_solutions[i] =  y_solution
    nu0s_num[i] = nu0_func(y_solution, delta)
    corr_num = corrective_term(y_solution, delta)
    corrections_num[i] = corr_num

# Plots

fig_F = figure(1, figsize=(4,3))
plot(deltas, corrections_num, 'k')
ylabel('$F(\Delta/L)$', fontsize=14)
xlabel('$\Delta/L$', fontsize=14)
ylim(0,0.2)
legend(frameon=False)

tight_layout()

show()






