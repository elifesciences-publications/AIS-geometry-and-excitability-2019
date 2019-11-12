## Models


* `model_spike_initiation.py`: a simple model of spike initiation in a ball-and-stick neuron, that considers only Na channels activation. 
Na channels are present only in the AIS.

* `model_Na_Kv1.py`: a simple biophysical model of an action potential,
with Na channels activation and inactivation, and Kv1-like channels activation.
The AIS contains a high density of both channels and the rest of the neuron a much lower density of both channels.

* `model_Na_Kv1_Kv7.py`: a simple biophysical model of an action potential,
with Na channels activation and inactivation, Kv1-like channels activation and Kv7-like channels activation.
The AIS contains a high density of both channels and the rest of the neuron a much lower density of both channels.
Kv7 channels are present only at the distal half of the AIS. 

* `params_model_description`: the parameters used to present the biophysical model of the AP in the section "Measuring excitability"

* `params_all`: the parameters used for all the other figures. 
The only parameter values that change are the Na and Kv1 density at the soma and the dendrites. 
