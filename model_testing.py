


'''
Just a demo of how to use the model
'''

from __future__ import print_function
from brian2 import *
from shared import params_model_description, model_Na_Kv1, model_Na_Kv1_myelin, model_Na_Kv1_axodendritic, measure_current_threshold, measure_voltage_threshold
from joblib import Parallel, delayed
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.ticker as ticker
from scipy import stats


defaultclock.dt = 0.005*ms

params = params_model_description
start =5.*um
end = 35.*um

# neuron = model_Na_Kv1_myelin(params, -75.*mV, start, end, density=False)
neuron = model_Na_Kv1_axodendritic(params, -75.*mV, start, end, density=False)


M = StateMonitor(neuron, ('v','I_VC', ), record = 0)
M_AIS = StateMonitor(neuron, ('v'), record = neuron.morphology.stem.axon[end])


# Protocol to elicit an AP with fixed initial voltage
neuron.V_VC[0] = -75.*mV
neuron.VC_on[0] = 0
run(20*ms)
neuron.VC_on[0] = 0
neuron.I_CC[0] = 2.*nA
run(2*ms)
neuron.VC_on[0] = 0
neuron.I_CC[0] = 0*amp 
run(20*ms)


# Plotting
figure('AP', figsize=(10,4))

subplot(111)
plot(M.t/ms, M.v[0]/mV, 'k', label='soma') 
plot(M_AIS.t/ms, M_AIS.v[0]/mV, 'r', label='AIS')
ylabel('V (mV)', fontsize=14)
xlabel('Time (ms)', fontsize=14)
xlim(19, 25)
legend(frameon=False)

# subplot(122)
# plot(neuron.axon.distance[int(end/um):]/um, M_AIS.v[:,4200]/mV, 'r', label='AIS')
# ylabel('V (mV)', fontsize=14)
# xlabel('Time (ms)', fontsize=14)
# legend(frameon=False)

tight_layout()

show()