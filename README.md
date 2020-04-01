# Theoretical relation between axon initial segment geometry and excitability
Code corresponding to the paper
"[Theoretical relation between axon initial segment geometry and excitability](https://elifesciences.org/articles/53432)"
by Sarah Goethals and Romain Brette.

### Requirements

In addition to standard scientific packages (scipy, matplotlib), this code requires:
* joblib
* [Brian 2](http://briansimulator.org)

For Figure 3, electrophysiological data from Hu and Bean (2018)
(DOI: 10.5281/zenodo.3539297) must be
[downloaded](https://zenodo.org/record/3539297#.Xc0WbjJKhBw) into the folder
`shared/resistive_coupling/hu-bean-2018`.

## Passive properties

* `figure_2.py`: We show in simulations the difference between the small soma situation (=sealed end) and
the large soma situation (=killed end). The large soma situation is the starting point of
resistive coupling theory.
* `figure_3.py`: We show from Hu & Bean's data (2018) that cortical neurons follow resistive coupling theory.
This requires downloading electrophysiological data [link to be added].
* `figure_4.py`: We show that electrical properties are conserved with a power law scaling of axon diameter vs soma diameter. 
We show with data from various sources in the literature that this is roughly the
observed scaling. In other words, the geometry of neurons is consistent with resistive coupling theory.

## Measuring excitability

In this part, we show that excitability changes due to the AIS are reported in the somatic voltage threshold
(not the axonal threshold or the current threshold).

### A simple model of spike initiation

* `figure_5.py`: We present a simple biophysical model of spike initiation in the AIS, and illustrate how to use it.

### How to measure excitability?

* `figure_6.py`: We show that the measure of excitability that best captures the effect of the axon initial segment is the somatic voltage threshold, 
and not the current threshold or the voltage threshold at the AIS. It is possible to either run the simulations, 
or only plot the figure based on previously saved results (`figure_6.npz`).

## Relation between excitability and AIS geometry

Here we show how the somatic threshold depends on AIS position and length.
We also show that it depends on whether the AIS is hyperpolarized (eg by Kv channels).

### A point AIS

We show that either shifting a point AIS away from the soma or increasing the total Na conductance in the AIS decreases the voltage threshold.

* ` figure_7.py`: We compare the theoretical threshold with the threshold measured in a simple model of spike initiation and in a simple biophysical model of an AP.
It is possible to either run the simulations, or only plot the figure based on previously saved results (`figure_7.npz`).

### A spatially extended AIS

We show that threshold variations can be separated in three independent contributions: 
AIS length, AIS middle position and density of Nav channels in the AIS.

* `figure_8.py`: Relation between threshold and AIS start position and length in the biophysical model, for a fixed Na conductance density in the AIS. 
It is possible to either run the simulations, or only plot the figure based on previously saved results (`figure_8.npz`).

* `figure_9.py`: Effect of compressing the AIS around its middle position, for a fixed total Na conductance at the AIS. 
It is possible to either run the simulations, or only plot the figure based on previously saved results (`figure_9.npz`).

* `figure_10.py`: Dependence of the voltage threshold on AIS length, AIS middle position and density of Nav channels in the AIS.
It is possible to either run the simulations, or only plot the figure based on previously saved results (`figure_10.npz`).

### Non-sodium axonal currents
    
We show that if a strong hyperpolarizing current is present a the AIS, the threshold can increase when the AIS is shifted away from the soma.

* `figure_11.py`: We measure the voltage threshold in the biophysical model when a hyperpolarizing conductance (Kv7-like conductance) 
is added in the distal half of the AIS.

* `figure_12.py`: We compare the voltage threshold at the soma and at the AIS with the resting membrane potential at the soma and the AIS, 
for different point AIS positions and injected current in the AIS.

## Methods

* `figure_13.py`: This script computes the exact solution to the bifurcation problem to calculate the threshold. 
It illustrates how the threshold for an extended AIS differs from that of a point AIS located in the middle of it,
for the same number of total Na conductance.

In addition, two scripts illustrate in more details how the exact solution of the bifurcation problem depends on 
the AIS start and end position:

* `theory/threshold_ext_AIS.py` for an extended AIS starting away from the soma;
* `theory/threshold_ext_AIS_with_axonal_current.py` for an extended AIS starting at the soma, with a current injected at AIS end.

