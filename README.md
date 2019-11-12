# AIS-geometry-and-excitability-2019
Code corresponding to the paper "Theoretical relation between axon initial segment geometry and excitability" by Sarah Goethals and Romain Brette

## Geometry of the spike initiation system

First part of the paper:

Figure 2: we show in simulations the difference between the small soma situation (=sealed end) and
the large soma situation (=killed end). The large soma situation is the starting point of
resistive coupling theory.
Figure 3: we show from Hu & Bean's data that cortical neurons follow resistive coupling theory.
Figure 4: we show that electrical properties are conserved with a power law scaling
of axon diameter vs soma diameter. We show with data from various sources in the literature that this is roughly the
observed scaling. In other words, the geometry of neurons correspond to resistive coupling theory.

## Measuring excitability: the right threshold parameter is the somatic voltage threshold

In this part, we establish that what determines excitability for a somatic or dendritic input is
the somatic voltage threshold; as opposed to axonal threshold or somatic current threshold.

    A simple model of spike initiation:
Figure 5: we present a simple biophysical model of spike initiation in the AIS, and illustrate how to use the model.

    How to measure excitability?
Figure 6: we show that the measure of excitability that captures the best the effect of the axon initial segment is the somatic voltage threshold, 
and not the current threshold or the voltage threshold at the AIS. It is possbile to either run the simulations, 
either only plot the figure based on previously saved results (figure_6.npz).

## How the spike threshold changes with AIS plasticity

Here we show the somatic threshold depends on AIS position and length.
We also show that it depends on whether the AIS is hyperpolarized (eg by Kv).

    A point AIS

We show that either shifting a point AIS away from the soma or increasing the total Na conductance in the AIS decreases the voltage threshold.

Figure 7: we compare the theoretical threshold with the threshold measured in a simple model of spike initiation and in a simple biophysical model of an AP.
It is possbile to either run the simulations, either only plot the figure based on previously saved results (figure_7.npz).

    A spatially extended AIS

We show that the threshold variations can be separated in three independant contributions: 
the AIS length, the AIS middle position and the density of Na channels in the AIS.

Figure 8: relation between threshold and AIS start position and length in the biophysical model, for a fixed Na conductance density in the AIS. 
It is possbile to either run the simulations, either only plot the figure based on previously saved results (figure_8.npz).

Figure 9: effect of compressing the AIS around its middle position, for a fixed total Na conductance at the AIS. 
It is possbile to either run the simulations, either only plot the figure based on previously saved results (figure_9.npz).

Figure 10: dependence of the voltage threshold on the AIS length, the AIS middle position and the density of Na channels in the AIS.
It is possbile to either run the simulations, either only plot the figure based on previously saved results (figure_9.npz).

    Non-sodium axonal currents
    
We show that if strong hyperpolarizing current is present a the AIS, the trheshold can decrease when the AIS is shifted awway from the soma.

Figure 11: we measure the voltage threshold in the biophysical model when a hyeprpolarizing conductance (Kv7-like conductance) 
is added in the distal half of the AIS.

Figure 12: we compare the voltage threshold at the soma and at the AIS with the resting membrane potential at the soma and the AIS, 
for different point AIS positions and injected current in the AIS.

    Methods

Figure 13: computes the exact solution to the bifurcation problem to calculate the threshold. 
It illustrates how the threshold for an extended AIS differs from that of a point AIS located in the middle of it,
for the same number of total Na conductance.

In addition, two scripts illustrate in more details how the exact solution of the bifurcation problem depends on 
the AIS start and end position:

theory/threshold_ext_AIS: for an extended AIS starting away from the soma
theory/threshold_ext_AIS_with_axonal_current: for an extended AIS starting at the soma, with a current injected at AIS end.