'''

A script to plot the relationship between the AIS diameter and the soma diameter.
Measures of AIS diameter come from electron microscopy studies. AIS diameter is the 

'''

from brian2 import *
from scipy.stats import linregress
from numpy.random import choice

import matplotlib.ticker

### Importing digitized datasets
a_soma_sm, d_axon_sm = genfromtxt('data/figure_4/sasaki-murayama-1992-Fig2.csv', delimiter=',', skip_header=1).T
d_soma_sp, d_axon_sp = loadtxt('data/figure_4/slower-powell-1979-Fig71.txt', skiprows = 1).T
d_axon_cr, d_soma_cr = genfromtxt('data/figure_4/conradi-ronnevi-1977.csv', delimiter=',', skip_header=1).T

# Sazaki-Murayama                                 
d_soma_sm = sqrt(a_soma_sm/pi) # assuming a_soma_sm is indeed the area of the somatic membrane, as a sphere

### Data for population averages
# De Zeeuw et al, 1990 (cat olivary cells)
d_axon_dz = 1.1
# Ruigrok et al, 1990 (cat olivary cells)
d_soma_dz = 21.7

# Kosaka 1980 (rat CA3)
d_axon_k = 1.24
# Buckmaster 2012 (rat CA3) 
d_soma_k = 20.9

# Somogyu Hamori 1976 (rat Purkinje cell) 
d_axon_sh_rat = 0.73
# Takacs Hamori 1994 (rat Purkinje cell)
d_soma_sh_rat = 21.88

# Palay and Chan-Palay 2012 (Cerebellar granule cells)
d_axon_cgc = 0.2
# Delvendahl et al 2015 (Cerebellar granule cells)
d_soma_cgc = 5.9

### Power law fits
d_axon = array(list(d_axon_sm)+list(d_axon_sp)+list(d_axon_cr)+[d_axon_dz, d_axon_k, d_axon_sh_rat,d_axon_cgc])
d_soma = array(list(d_soma_sm)+list(d_soma_sp)+list(d_soma_cr)+[d_soma_dz, d_soma_k, d_soma_sh_rat,d_soma_cgc])
# Below: inter-set fits (just one point per data set)
#d_axon = array([mean(d_axon_sm),mean(d_axon_sp),mean(d_axon_cr),d_axon_dz, d_axon_k, d_axon_sh_rat,d_axon_cgc])
#d_soma = array([mean(d_soma_sm),mean(d_soma_sp),mean(d_soma_cr),d_soma_dz, d_soma_k, d_soma_sh_rat,d_soma_cgc])

d = linspace(5,50,100)
offset = mean(log(d_axon)-4./3*log(d_soma))
print ("Regression: d_AIS = ",exp(offset),"d_soma^(4/3)")

### Plot
fig_diam = figure(figsize=(5,4))
axes = subplot(111)
# full datasets
loglog(d_soma_sm,d_axon_sm,'o', color = 'darkorange',label='human motoneurons')
loglog(d_soma_sp,d_axon_sp,'o', color='forestgreen', label='primate pyramidal and stellate cells')
loglog(d_soma_cr,d_axon_cr,'o', color='darkblue', label='cat motoneuron')
# averages
loglog(d_soma_dz, d_axon_dz, 'X', color='deepskyblue', label='cat olivary cells', markersize=7)
loglog(d_soma_k, d_axon_k, 's',  color='deepskyblue',label='rat CA3 pyramidal cells', markersize=7)
loglog(d_soma_sh_rat, d_axon_sh_rat, '^',  color='deepskyblue',label='rat Purkinje cells', markersize=7)
loglog(d_soma_cgc, d_axon_cgc, 'D',  color='deepskyblue',label='mouse granule cells', markersize=7)
# electrical equivalence
loglog(d, exp(offset + 4./3*log(d)),'k--', label='electrical equivalence')
ylabel('AIS diameter ($\mu$m)')#, fontsize=20)
xlabel('Soma diameter ($\mu$m)')#, fontsize=20)
axes.xaxis.set_major_formatter(ScalarFormatter())
axes.yaxis.set_major_formatter(ScalarFormatter())
yticks([1, 2, 3, 4])#, fontsize = 16)
xticks([10,20,30,40,50])#, fontsize = 16)

# legends
lines = axes.get_lines()
legend1 = legend([lines[i] for i in [0,1,2]], ['human motoneurons', 'primate pyramidal and stellate cells', \
                 'cat motoneuron'], loc='upper left', frameon=False)
legend2 = legend([lines[i] for i in [3,4,5,6,7]], ['cat olivary cells', 'rat CA3 pyramidal cells', 'rat Purkinje cells'\
                 , 'mouse granule cells', 'electrical equivalence'], loc='lower right', frameon=False)
axes.add_artist(legend1)
axes.add_artist(legend2)

tight_layout()

show()
