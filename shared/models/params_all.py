
"""
Parameters for the biophysical model of an action potential.
"""

from brian2.units import *

### Morpho
axon_diam = 1. * um # axon diameter

### Passive parameters
EL = -75. * mV # resting membrane potential
Cm = 0.9* uF / cm ** 2 # capacitance per unit area of membrane
gL = 1./(15000. * ohm * cm**2) # leak conductance density
Ri = 100. * ohm * cm # intracellular resistivity

### Channels 
T = 33. # simulation temperature
factor = (1. / 2.8) ** ((T - 23.) / 10.) # correction for the temperature: we multiply tau by Q10^-(T-23)/10 because it is equivalent to multpily the rates by  Q10^(T-23)/10

## AIS
# Na channels 
ENa = 70. * mV # reversal potantial of Na conductance
gna_dens = 3500. * (siemens / meter ** 2) # Na conductance surfacic density in the AIS
Gna = 350.* nS # total Na conductance in the AIS

Va = -35. * mV  # half-activation voltage
Ka = 5. * mV  # slope of the activation curve
Taum_max = factor * 0.15 * ms  # maximal time constant of activation
Vh = -65.*mV # half-inactivation voltage
Kh = 5. * mV  # slope of the inactivation curve
Tauh_max =  factor * 5. * ms # maximal time constant of inactivation

# K1 channels 
EK = -90. * mV # reversal potantial of K conductance
gk_dens = 1500. * (siemens / meter ** 2) # K conductance surfacic density in the AIS
Gk = 150.*nS # total K conductance in the AIS

Vn = -70.00001 * mV  # # half-activation voltage (v=nan when V0 = -70mV = Vn)
Kn = 20. * mV  # slope of the activation curve
Taun_max = 1. * ms # maximal time constant of activation

# K7 channels 
EM = -90. * mV # reversal potantial of K7 conductance
gm_dens = 150. * (siemens / meter ** 2) # K7 surfacic conductance density in soma

## Soma
# Na channels 
gna_soma = 500. * (siemens / meter ** 2) # Na surfacic conductance density in soma
gna_dend = 20. * (siemens / meter ** 2) # Na surfacic conductance density in dendrites and axon outside the AIS

Va_soma = -30.*mV # half-activation voltage
Vh_soma = -60.*mV  # half-inactivation voltage

# K channels 
gk_soma = 500. * (siemens / meter ** 2) # K surfacic conductance density in soma
gk_dend = 20. * (siemens / meter ** 2) # K surfacic conductance density in dendrites and axon outside the AIS

### Prediction
y0 = 1.19968


