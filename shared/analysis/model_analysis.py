
"""

Functions to measure electrical properties in neuron models.

"""
from brian2 import *

__all__ = ['calculate_resting_state', 'measure_input_resistance', 'measure_voltage_threshold',  'measure_current_threshold',  \
           'measure_voltage_threshold_no_VC',  'measure_current_threshold_no_VC','measure_threshold_in_vc']

defaultclock.dt = 0.005*ms

def calculate_resting_state(neuron, Na_start=5.*um, Na_end=30.*um, hyperpol_current = None):
    '''
    Measure the resting membrane potential of a neuron model, at the soma and at the AIS end.
    '''
    M = StateMonitor(neuron, ('v'), record = 0)
    M_AIS = StateMonitor(neuron, ('v'), record = neuron.morphology.axon[Na_end]) 
    
    if hyperpol_current is None :
        run(150.*ms)
    else: 
        ais_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int((Na_end)/um)
        neuron.I_CC[ais_idx] = hyperpol_current
        run(150*ms)
        
    return M.t/ms, M.v[0]/mV, M_AIS.t/ms, M_AIS.v[0]/mV

def measure_input_resistance(params, neuron, resting_vm=-75.*mV):
    '''
    Measure the input resistance at the soma, as the voltage gradient at the end of a current pulse divided by the amplitude of the current pulse.
    '''
    
    M = StateMonitor(neuron, ('v','I_VC'), record = 0)
    
    current = 0.1*nA
    
    neuron.V_VC[0] = resting_vm 
    neuron.VC_on[0] = 1
    run(20*ms)
    neuron.VC_on[0] = 0
    neuron.I_CC[0] = current
    run(200.*ms)
    neuron.I_CC[0] = 0*amp 
    run(20*ms)
    
    Rin_soma = (M.v[0][int((20*ms+200*ms)/defaultclock.dt)]-resting_vm) / current
    
    return Rin_soma, M.v[0]/mV, M.t/ms
    
def measure_current_threshold(params, neuron, resting_vm=-75.*mV, ais_start=5.*um, ais_end=30.*um, pulse_length = 2.*ms, i_inj = 0, si_site = False, latency = 20.*ms, plot_v = True, color='k'):
    '''
    Measure the current threshold as the minimal current that elicits a spike, for variable duration of the current pulse,
    with the bissection method. The soma is voltage-clamped until the injection of the current pulse.
    
    Returns the minimal current amplitude that elicits a spike at the AIS.    
    '''
    
    if ais_start > ais_end:
        print ('The AIS starts after it ends.')
        i_current = nan
    else:      
        i_max = 3.*nA
        i_min = 0.*nA
        i_current = 1.5*nA
        spike = False
        
        M = StateMonitor(neuron, ('v','I_VC'), record = 0)
        
        if si_site:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[si_site])
        else:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_end])
            
        store()
        
        # THe location of current injection
        inj_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int(ais_end/um) 
       
        while True:
            print (i_current)
            
            restore()
            
            neuron.V_VC[0] = resting_vm 
            neuron.VC_on[0] = 1
            neuron.I_CC[inj_idx] = i_inj
            run(latency)
            neuron.VC_on[0] = 0
            neuron.I_CC[0] = i_current
            neuron.I_CC[inj_idx] = i_inj
            run(pulse_length)
            neuron.I_CC[0] = 0*amp 
            neuron.I_CC[inj_idx] = i_inj
            run(20*ms)
            
            figure('measure current threshold', figsize=(8,7))
            plot(M.t/ms, M.v[0]/mV, color=color)
            plot(M_AIS.t/ms, M_AIS.v[0]/mV, '--', color=color)
            
            m_max_ais = max(M_AIS.m[0])

            if m_max_ais >= 0.5 and abs(i_current - i_min) <= 0.1*pA and spike == False :
                print ('stop')
                plot(M.t/ms, M.v[0]/mV, 'r', label='V soma at threshold')
                plot(M_AIS.t/ms, M_AIS.v[0]/mV, '--r', label='V AIS at threshold')
                legend()
                break
            if m_max_ais <= 0.5:
                i_min = i_current
                spike = False
            else: 
                i_max = i_current
                spike = True
            
            i_current = 0.5*i_max + 0.5*i_min
                
    return i_current

def measure_voltage_threshold(params, neuron, resting_vm=-75.*mV, ais_start=5.*um, ais_end=30.*um, i_rheo = 0.5*nA, pulse_length = 2.*ms, i_inj = 0, latency = 20.*ms, si_site = False, plot_v = True, color = 'k'):
    '''
    Measures the voltage threshold as the maximal membrane potential reached during the largest non-spiking situation.
    The soma is voltage-clamped at the resting membrane potential, 
    then a current pulse of an amplitude just below the current threshold is injected at the soma.
    
    Returns: voltage threshodl at the soma and at the AIS end.
    '''
    if ais_start > ais_end:
        print ('The AIS starts after it ends.')
        vth_soma = nan
        vth_ais = nan
    else:   
        M = StateMonitor(neuron, ('v','I_VC'), record = 0)
        
        if si_site:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[si_site])
        else:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_end])
        
        inj_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int(ais_end/um)
        
        neuron.V_VC[0] = resting_vm 
        neuron.VC_on[0] = 1
        neuron.I_CC[inj_idx] = i_inj
        run(latency)
        neuron.VC_on[0] = 0
        neuron.I_CC[inj_idx] = i_inj
        neuron.I_CC[0] = i_rheo * 0.999
        run(pulse_length)
        neuron.I_CC[inj_idx] = i_inj
        neuron.I_CC[0] = 0*amp 
        run(20*ms)
        
        vth_soma = max(M.v[0]) 
        vth_ais = max(M_AIS.v[0]) 
                     
    return vth_soma, vth_ais, M.v[0], M_AIS.v[0]

def measure_current_threshold_no_VC(params, neuron, resting_vm=-75.*mV, ais_start=5.*um, ais_end=30.*um, pulse_length = 2.*ms, i_inj = 0, si_site = False, plot_v = True, color='k'):
    '''
    Same as 'measure_current_threshold' without voltage-clamping the soma before the current pulse.
    '''
    
    if ais_start > ais_end:
        print ('break')
        i_current = nan
    else:      
        i_max = 4.*nA
        i_min = 0.*nA
        i_current = 2.*nA
        spike = False
        
        M = StateMonitor(neuron, ('v','I_VC'), record = 0)
        
        if si_site:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[si_site])
        else:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_end])
            
        store()
        
        inj_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int(ais_end/um) 
       
        while True:
            print (i_current)
            
            restore()
            
            neuron.VC_on[0] = 0
            neuron.I_CC[inj_idx] = i_inj
            run(80.*ms)
            neuron.VC_on[0] = 0
            neuron.I_CC[0] = i_current
            neuron.I_CC[inj_idx] = i_inj
            run(pulse_length)
            neuron.I_CC[0] = 0*amp 
            neuron.I_CC[inj_idx] = i_inj
            run(20*ms)
            
            figure('measure current threshold', figsize=(8,7))
            plot(M.t/ms, M.v[0]/mV, color=color)
            plot(M_AIS.t/ms, M_AIS.v[0]/mV, '--', color=color)
            
            m_max_ais = max(M_AIS.m[0])

            if m_max_ais >= 0.5 and abs(i_current - i_min) <= 0.1*pA and spike == False :
                print ('stop')
                plot(M.t/ms, M.v[0]/mV, 'r', label='V soma at threshold')
                plot(M_AIS.t/ms, M_AIS.v[0]/mV, '--r', label='V AIS at threshold')
                legend()
                break
            if m_max_ais <= 0.5:
                i_min = i_current
                spike = False
            else: 
                i_max = i_current
                spike = True
            
            i_current = 0.5*i_max + 0.5*i_min
                
    return i_current

def measure_voltage_threshold_no_VC(params, neuron, resting_vm=-75.*mV, ais_start=5.*um, ais_end=30.*um, i_rheo = 0.5*nA, pulse_length = 2.*ms, i_inj = 0, si_site = False, plot_v = True, color = 'k'):
    '''
    Same as 'measure_voltage_threshold' without voltage-clamping the soma before the current pulse.
    '''
    if ais_start > ais_end:
        print ('break')
        th_rheo = nan
    else:   
        M = StateMonitor(neuron, ('v','I_VC'), record = 0)
        
        if si_site:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[si_site])
        else:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_end])
        
        inj_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int(ais_end/um)
        
        neuron.VC_on[0] = 0
        neuron.I_CC[inj_idx] = i_inj
        run(80.*ms)
        neuron.VC_on[0] = 0
        neuron.I_CC[inj_idx] = i_inj
        neuron.I_CC[0] = i_rheo * 0.999
        run(pulse_length)
        neuron.I_CC[inj_idx] = i_inj
        neuron.I_CC[0] = 0*amp 
        run(20*ms)
        
        th_rheo_soma = max(M.v[0]) #int(22.*ms/defaultclock.dt)])
        th_rheo_ais = max(M_AIS.v[0]) 
        v0_soma = mean(M.v[0][int(60.*ms/defaultclock.dt):int(80.*ms/defaultclock.dt)])
        v0_ais = mean(M_AIS.v[0][int(60.*ms/defaultclock.dt):int(80.*ms/defaultclock.dt)])
                     
    return th_rheo_soma, th_rheo_ais, M.v[0], M_AIS.v[0], v0_soma, v0_ais

def measure_threshold_in_vc(params, neuron, ais_start=5.*um, ais_end=30.*um, hyperpol_current = False, plot_v = True, color='k'):
    
    """
    Measures the threshold in voltage-clamp as the lowest somatic membrane potential at the soma that triggers a spike at AIS end. 
    """

    if ais_start > ais_end:
        print ('break')
        v_current = nan
    else:
        
        peak = 0.*mV
        v_max = -25.*mV
        v_min = -75.*mV
        v_current = -50.*mV
        v_previous = v_current
        spike = False
        
        M = StateMonitor(neuron, ('v','I_VC'), record = 0)
        M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_end]) 
        
        store()
        
        while True:
            restore()
            neuron.V_VC[0] = params.EL 
            neuron.VC_on[0] = 1
            run(20*ms)
            neuron.V_VC[0] = v_current
            neuron.VC_on[0] = 1
            run(50.*ms)
            neuron.VC_on[0] = 0
            run(20*ms)

            v_max_ais = max(M_AIS.v[0][int(20.*ms/defaultclock.dt):int(70.*ms/defaultclock.dt)])
            
            if v_max_ais >= peak and abs(v_current - v_min) <= .1*mV and spike == False :
                print( 'spike:', v_current)
                break
            if v_max_ais < peak:
                v_min = v_current
                spike = False
            else: 
                v_max = v_current
                spike = True
            
            v_previous = v_current
            v_current = 0.5*v_max + 0.5*v_min
        
    return v_current, v_previous


    
    
    
    
    
    
    
    
    
    
    
    

