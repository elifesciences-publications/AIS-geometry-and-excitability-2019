'''

A script to read the data (AIS diameter and soma diameter) from Sloper and Powell 1979. 

'''


from brian2 import *
from scipy.stats import linregress
from numpy.random import choice

d_soma, d_axon = loadtxt('slower-powell-1979-Fig71.txt', skiprows = 1).T

# Regression in log plot
slope, intercept, r_value, p_value, std_err = linregress(log(d_soma), log(d_axon))
print slope

# fit with slope 4/3
d = linspace(10,40,100)
offset = mean(log(d_axon)-4./3*log(d_soma))

figure()
loglog(d_soma, d_axon, 'ko')
loglog(d_soma, exp(intercept + slope*log(d_soma)),'r', label='linear regression in log')
loglog(d, exp(offset + 4./3*log(d)),'b', label='fit with slope 4/3')

ylabel('AIS diameter (um)')
xlabel('Soma diameter (um)')
legend(frameon=False)
tight_layout()

show()


## Bootstrap analysis
#slopes = zeros(1000)
#for i in range(1000):
#    bootstrap = choice(arange(len(d_soma)), len(d_soma))
#    slope, intercept, r_value, p_value, std_err = linregress(log(d_soma[bootstrap]), log(d_axon[bootstrap]))
#    slopes[i] = slope
#print mean(slopes),std(slopes)
#
#xlim(0,40)
#ylim(0.3,5)
#
#figure()
#plot(d_soma, d_axon, '.k')
#plot(d, exp(offset + 4./3*log(d)),'k')
#plot(d, exp(offset2 + 3./2*log(d)),'b')
## Best linear fit
##linear_slope = sum(d_axon*d_soma)/sum(d_soma*d_soma)
## Best linear fit, calculated in log space (ie for relative precision)
#plot(d, linear_slope*d,'k--')
#xlim(0,40)
#