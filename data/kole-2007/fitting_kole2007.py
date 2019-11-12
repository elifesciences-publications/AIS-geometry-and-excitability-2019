

'''
Fitting the activation curve of Kv1 channels from Kole et al, Neuron 2007. 
Data points where digitized from figure 2C. 
'''
from brian2 import *
from pylab import *
from scipy.optimize import curve_fit
import pandas as pd

data = pd.read_csv("kole2007fig2C.csv", decimal = ",") 

vm = data.x.values
ninf8 = data.Curve1.values

# Boltzmann fit
f = lambda x,vn,k : (1/(1+exp((vn-x)/k)))**8

subplot(111)
popt,_ = curve_fit(f, vm, ninf8, p0=[-70.,20.])
print popt
vn,k = popt

plot(vm,ninf8,'k')
xlim(-100, 50)
ylim(0,1)
plot(vm,f(vm,vn,k),'r')
xlabel('V (mV)')
ylabel('ninf8')

show()

