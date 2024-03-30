import curie as ci
import numpy as np
import matplotlib.pyplot as plt

# def basic_examples():

#     ### We will plot the same reaction from two different libraries
#     ### Passing f,ax to rx.plot allows multiple plots on the same figure

#     rx = ci.Reaction('90ZR(n,2n)89ZR', 'irdff')
#     f,ax = rx.plot(return_plot=True, label='library')
#     rx = ci.Reaction('90ZR(n,2n)89ZR', 'endf')
#     rx.plot(f=f,ax=ax, label='library')


#     ### Compare (n,2n) and (n,3n) for endf vs tendl
#     f, ax = None, None
#     for lb in ['endf','tendl']:
#         rx = ci.Reaction('226RA(n,2n)225RA', lb)
#         f, ax = rx.plot(f=f, ax=ax, return_plot=True, label='both', energy=np.arange(0,30,0.1))
#         rx = ci.Reaction('226RA(n,3n)224RA', lb)
#         f, ax = rx.plot(f=f, ax=ax, return_plot=True, label='both', energy=np.arange(0,40,0.1))

#     plt.show()

#     # ### Search the TENDL-2015 neutron library for reactions producing 225RA from 226RA
#     lb = ci.Library('tendl_n')
#     print(lb.search(target='226RA', product='225RAg'))



# basic_examples()


rx = ci.Reaction('*ZR(d,x)87Y')
# rx = ci.Reaction('*ZR(d,x)86Y')
# rx = ci.Reaction('*FE(d,x)*')
# rx = ci.Reaction('0LA(p,x)*')
# rx = ci.Reaction('0FE(d,x)*CO')



rx.search_exfor(plot_results=True,plot_tendl=True, show_legend=True, xlim=[None,None], verbose=True)
print(rx.plot_Dict)
# print(vars(rx))