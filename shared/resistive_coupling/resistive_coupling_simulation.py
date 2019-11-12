'''
Measurement of passive responses to current pulses in a ball-and-stick model.
Input resistances are calculated.

The neuron is a sphere + a long cylindrical axon
'''
from brian2 import *
from tempfile import mkdtemp
from joblib import Memory

cachedir = mkdtemp()
memory = Memory(cachedir=cachedir)

## Parameters

defaultclock.dt = 0.1*ms

EL = -75*mV # resting membrane potential
gL = 1./(15000 * ohm * cm**2) # leak conductance density
#print ("Time constant=", 0.9*uF/cm**2/gL)

axon_diameter = 1*um
Ri = 100.*ohm*cm # intracellular resistivity

def pulse_simulation(soma_diameter, distance, current = 1*nA):
    '''
    Simulates a ball-and-stick model with a 200 ms current pulse injected at the specified distance.

    Returns: t, x, V(soma,t), V(axon,t), V(x,end of pulse)
    '''
    morpho = Soma(soma_diameter)
    morpho.axon = Cylinder(diameter=axon_diameter, length=2000*um, n=1000)

    eqs='''
    Im = gL*(EL-v) : amp/meter**2
    I : amp (point current)
    '''

    neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=0.9*uF/cm**2, Ri=Ri, method="exponential_euler")
    neuron.v = EL
    injection_site = morpho.axon[distance]

    Msoma = StateMonitor(neuron, 'v', record = 0)
    Maxon = StateMonitor(neuron, 'v', record = injection_site)

    run(20*ms)
    neuron.I[injection_site] = current
    run(200*ms)
    v = neuron.v[:].copy()
    neuron.I[injection_site] = 0*nA
    run(100*ms)

    x = neuron.distance

    return Msoma.t,x,Msoma.v[0],Maxon.v[0],v


#@memory.cache # doesn't work!
def input_resistance(soma_diameter, distance, current = 1*nA):
    '''
    Calculates the input resistance at 300 us and 200 ms for an injection at a given distance
    from the soma. The resistance is calculated both at axonal injection site and soma.
    '''
    t, _, vs, va, _ = pulse_simulation(soma_diameter, distance, current)

    Ra300us = (va[int((20*ms+1000*us)/defaultclock.dt)]-EL) / current
    Rs300us = (vs[int((20*ms+1000*us)/defaultclock.dt)]-EL) / current
    Ra200ms = (va[int((20*ms+200*ms)/defaultclock.dt)]-EL) / current
    Rs200ms = (vs[int((20*ms+200*ms)/defaultclock.dt)]-EL) / current

    return Ra300us,Rs300us,Ra200ms,Rs200ms

