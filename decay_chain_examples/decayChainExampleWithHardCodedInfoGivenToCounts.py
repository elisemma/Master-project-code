import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import curie as ci

def plot_decay_chain(parent_isotope, units = 'd'):
    dc = ci.DecayChain(parent_isotope, A0=3E4, units=units)
    dc.counts = {'48V':[[5.0, 5.1, 6E5, 2E4], [6.0, 6.1, 7E5, 3E4]]} #This data is just some dummy data
    dc.fit_A0()
    print('A0: ',dc.A0)
    dc.plot()


plot_decay_chain('48V')