'''
A script to read the data (AIS diameter and soma diameter) from Conradi and Ronnevi 1997.
'''

from brian2 import *
from scipy.stats import linregress

d_axon, d_soma = genfromtxt('conradi-ronnevi-1977.csv', delimiter=',', skip_header=1).T

slope, intercept, r_value, p_value, std_err = linregress(log(d_soma), log(d_axon))
print('Linear regression slope (log space): {}'.format(slope))

# fit with slope 4/3
offset2 = mean(log(d_axon)-4./3*log(d_soma))

figure()       
loglog(d_soma,d_axon,'ko')
           
d = linspace(20,60,100)
loglog(d, exp(intercept + slope*log(d)),'b', label='linear regression')
loglog(d, exp(offset2 + 4./3*log(d)),'r', label='fit with slope 4/3')

ylabel('AIS diameter (um)')
xlabel('Soma diameter (um)')
legend(frameon=False)
tight_layout()

show()

