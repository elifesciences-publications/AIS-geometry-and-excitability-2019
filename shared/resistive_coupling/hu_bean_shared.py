'''
Shared analysis functions for Hu and Bean's data (2018)
'''
from __future__ import print_function
import os
from .abf import ABF
from numpy import zeros,array,median
import re

__all__ = ['load_all_files','load_file','path']

path = "/Users/sarah/Documents/Data/Bruce Bean/Fig 7 Data Soma-Axon Coupling/" # to be changed!

def load_file(filename):
    '''
    Returns distance, time, somatic voltage, axonal voltage, somatic current, axonal current
    Voltage in mV, current in nA
    '''
    # Distance
    m = re.search('\s(\d+\.?\d*)\s*um', filename)
    try:
        distance = float(m.group(1))
    except:
        print("problem with", filename)

    abf = ABF(filename)
    abf.setSweep(0,channel=0)
    Vs = abf.dataY
    abf.setSweep(0,channel=1)
    Is = abf.dataY * 0.020 # 1 = 20 pA
    abf.setSweep(0,channel=2)
    Va = abf.dataY
    abf.setSweep(0,channel=3)
    Ia = abf.dataY * 0.020
    t = abf.dataX
    if filename == '140107 1-1 48um.abf': # axon and soma were reversed
        Vs,Va = Va,Vs
        Is,Ia = Ia,Is

    # Find the current zero and fix it
    Is0 = median(Is)
    Ia0 = median(Ia)
    Is=Is-Is0
    Ia=Ia-Ia0

    return distance,t,Vs,Va,Is,Ia

def load_all_files():
    '''
    Returns the data of each file one by one,
    using generator syntax. To be used in a for loop.

    The output is a list of distance, t, Vs, Va, Is, Ia
    '''
    files = [f for f in os.listdir(path)]

    for i,filename in enumerate(files):
        if filename[-4:] == '.abf':
            #print("Loading",filename)
            try:
                yield load_file(path + '/' + filename)
            except:
                print("Could not load",filename)

if __name__ == '__main__':
    from pylab import *

    distance, t, Vs, Va, Is, Ia = load_file(path+'140107 1-1 48um.abf')
    print("dt=",t[1]-t[0])

    subplot(211)
    plot(t,Vs,'r')
    plot(t,Va,'b')
    subplot(212)
    plot(t,Is,'r')
    plot(t,Ia,'b')
    show()
