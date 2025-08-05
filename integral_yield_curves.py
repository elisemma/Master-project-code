import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               # AutoMinorLocator)

import scipy.interpolate as inter
from scipy import interpolate

import curie as ci
from collections import OrderedDict
# params = {'text.usetex': True, 'mathtext.fontset': 'stix'}
# plt.rcParams["font.family"] = "Computer Modern"
# plt.rcParams.update(params)

import sys
# print(sys.version)


#  INPUT
########################################################################

el = ci.Element('Zr')     # Target material

particle='d'    # Incident beam

itp_list = ['86Y','87Y', '87Y', '88Y', '88Y', '90Ym']    # Product isotopes, in Curie notation
# i_list = ['86Y Cumulative','87Y Cumulative', '87Y Inependent', '88Y Cumulative', '88Y Independent',  '90mY Cumulative']    # Product isotopes, for the plot legend
A_list = [86, 87, 87, 88, 88, 90]
states = ['', '', '', '', '', 'm']
xs_types = ['Cumulative', 'Cumulative', 'Independent', 'Cumulative', 'Independent', 'Cumulative']

# Set up labels for plot legends
# prods = [r'$^{'+i[:-2]+r'}$'+i[-2:] for i in i_list]
# prods = [r'$^{\text{nat}}$Zr(d,x)$^{'+i[:-2]+r'}$'+i[-2:] for i in i_list]
# prods = [r'$^{\text{nat}}$Zr(d,x)$^{{A_list[i]}{states[i]}}$Y - {xs_types[i]}' for i in range(len(i_list))]
print(len(A_list))
print(len(xs_types))
print(len(states))


# prods = [f'$^{{\\text{{nat}}}}$Zr(d,x)$^{{{A_list[i]}{states[i]}}}$Y - {xs_types[i]}' for i in range(len(A_list))]
prods = [f'$^{{\\text{{nat}}}}$Zr(d,x)$^{{{A_list[i]}\\text{{{states[i]}}}}}$Y - {xs_types[i]}' for i in range(len(A_list))]



t_irr = 1   # irradiation length, in h  (only affects EOB thick target yield)

Z_beam = 1   # charge state of beam

show_all_plots = False  # Show all 4 yield plots, or only physical yield


xx = [np.array([6.346733794, 12.44547236, 19.24171826, 24.8961017, 26.39531264, 31.04155407,  36.351075, 41.11050957,  47.63516853]), 
      np.array([6.346733794, 12.44547236, 19.24171826, 24.8961017, 26.39531264, 31.04155407,  36.351075, 41.11050957,  47.63516853]),
      np.array([6.346733794, 12.44547236, 19.24171826, 24.8961017, 26.39531264, 31.04155407,  36.351075, 41.11050957,  47.63516853]),
	  np.array([6.346733794, 12.44547236, 19.24171826, 24.8961017, 26.39531264, 31.04155407,  36.351075, 41.11050957,  47.63516853]),
	  np.array([6.346733794, 12.44547236, 19.24171826, 24.8961017, 26.39531264, 31.04155407,  36.351075, 41.11050957,  47.63516853]),
	  np.array([6.346733794, 12.44547236, 19.24171826, 24.8961017, 26.39531264, 31.04155407,  36.351075, 41.11050957,  47.63516853])
	  ] # Energy [MeV]

yy = [np.array([0, 0, 0, 0, 0, 13.0295277, 33.1134087, 45.99303617,  39.07002542]),
	  np.array([0, 0, 10.00821349, 29.42720252, 28.86667139, 31.77747525, 24.35023946, 26.90301244, 39.3798355]),
	  np.array([0, 0, 3.1, 6.8, 7.0, 7.1, 6.0, 6.1, 7.3]),
	  np.array([0.429927518, 6.120008555, 9.135969156, 12.44828135, 15.59751062, 19.77940816, 80.31699393, 225.7576172, 353.4354775]),
	  np.array([0, 0, 0, 11.75, 0, 15.4, 21.4, 29.7, 45.4]),
	  np.array([0, 0.22851936, 0, 3.2853654, 3.44317019, 0, 5.271888315, 5.426059907, 6.20271661])
	  ] # [mb]
# xx = [np.array([47.63516853]), 
#       np.array([47.63516853]),
#       np.array([47.63516853]),
# 	  np.array([47.63516853]),
# 	  np.array([47.63516853]),
# 	  np.array([47.63516853])
# 	  ] # Energy [MeV]

# yy = [np.array([39.07002542]),
# 	  np.array([39.3798355]),
# 	  np.array([7.3]),
# 	  np.array([353.4354775]),
# 	  np.array([45.4]),
# 	  np.array([6.20271661])
# 	  ] # [mb]



# OPTIONAL
########################################################################
line_style = ['solid',(0,(5,6)), (0,(10,6)), (0,(5,4,1,4,1,4)), ':', 'solid']
# '-' or 'solid'	solid line
# '--' or 'dashed'	dashed line
# '-.' or 'dashdot'	dash-dotted line
# ':' or 'dotted'

line_weight = [2, 1, 1, 1.3, 2, 1] 

colors = ['hotpink', 'cornflowerblue', 'mediumaquamarine', 'mediumorchid', 'deepskyblue', 'grey']

# Above are plot options for the lines 



########################################################################
def read_exfor(fname):
	widths = np.array([0, 15, 28, 41, 54, 57, 82, 85, 95, 104])
	widths = np.diff(widths)
	names = ["x", "dx", "y", "dy", "#1",
             "year, auth", "#2", "id", "com"]

	dtype = ""
	for i, w in enumerate(widths):
		if i < 4:
			dtype += f"float,"
		else:
			dtype += f"S{w},"
	dtype = dtype[:-1]

	data = np.genfromtxt(fname, delimiter=widths,
                         skip_header=11, skip_footer=2,
                         dtype=dtype, names=names, deletechars="",
                         comments=None)

	a = [list(item) for item in data]
	data = pd.DataFrame(a, columns=names)
	data = data.drop(columns=["#1", "#2"])

	data["year, auth"] = data["year, auth"].str.decode('utf-8')
	data["com"] = data["com"].str.decode('utf-8')
	# data["id"] = data["id"].astype(int)
	return data

########################################################################
# data = read_exfor("natGe72Se_EXFOR.txt")     # Only if you wanted to overlay yield data from EXFOR, which is currently commented out


# show_all_plots = True

if show_all_plots:
	fig = plt.figure(figsize=(12, 9))

for n in np.arange(0,len(yy)):

	## If I was using my own experimental data instead of an EXFOR import I would have to just make sure the energy points were in increasing order
	# xx = np.concatenate((E_BNL,E_LANL))[::-1]
	# yy = np.concatenate((cs_BNL,cs_LANL))[::-1]
	# plt.plot(xx[n],yy[n])
	# plt.show()
	# print(i_list[n])


	smooth_xs = yy[n]*(1.0e-27)
	e_range = xx[n]

	MM = el.mass

	dc = ci.Isotope(itp_list[n]).decay_const('h',True)   #1/h
	# print('dc: ',dc[0])


	y = np.empty_like(e_range)
	S_pw= el.S(e_range,particle=particle, density=1E-3)*1E3   # MeV/(g/cm^2)
	# print(el.S([1e-3, 1e-2],particle=particle, density=1E-3)*1E3 )



	for i, item in enumerate(yy[n]): 
		# print(i)
		# print(e_range[0:i+1] )
		# print(np.concatenate(([], e_range[0:i+1])))

		y[i] = np.trapz( (np.concatenate(([0, 0], smooth_xs[0:i+1])) / (Z_beam*1.602e-19)) / (np.concatenate((el.S([1e-3, e_range[0]-0.2],particle=particle, density=1E-3)*1E3 , S_pw[0:i+1]))), x=np.concatenate(([0.1, e_range[0]-0.2], e_range[0:i+1]))   ) *(6.022e23/MM) / (1E6 * 1E6)
		#                     cm2             # particles/A      MeV/(g/cm^2)       MeV              # nuclei       MBq    uA

	# y[0]=0
	y = np.concatenate(([0], y))
	e_range = np.concatenate(([0], e_range))
	EOB_TTY = (1-np.exp(-dc[0]*t_irr))  * y
	saturation_TTY = y
	physical_TTY = dc[0] *y
	print(physical_TTY)
	
	
	# print(physical_TTY)


	if show_all_plots:
		plt.subplot(2, 2, 1)
		plt.plot(e_range,EOB_TTY,label=prods[n], linestyle=line_style[n], linewidth=line_weight[n])
		plt.ylabel('t$_{irr}$ ('+str(t_irr)+' hr) EOB Thick Target Yield (MBq/uA)')#,labelpad=7,fontsize=22)
		plt.xlabel('Beam Energy (MeV)')#,labelpad=7,fontsize=22)
		# plt.yscale('log')
		plt.xlim(left=0)
		plt.legend(loc='best')

		plt.subplot(2, 2, 2)
		plt.plot(e_range,saturation_TTY,label=prods[n], linestyle=line_style[n], linewidth=line_weight[n])
		plt.ylabel('Saturation Thick Target Yield (MBq/uA)')#,labelpad=7,fontsize=22)
		plt.xlabel('Beam Energy (MeV)')#,labelpad=7,fontsize=22)
		# plt.yscale('log')
		plt.xlim(left=0)
		plt.legend(loc='best')

		plt.subplot(2, 2, 3)

	plt.plot(e_range,physical_TTY/0.0036,label=prods[n], color=colors[n], linestyle=line_style[n], linewidth=line_weight[n])
	plt.ylabel('Physical Thick Target Yield (MBq/C)', fontsize=14)#,labelpad=7,fontsize=22)
	plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
	plt.yscale('log')
	plt.xlim(left=0)
	# plt.ylim(bottom=1e-3)
	# plt.ylim(top=1e7)
	plt.legend(loc='best', fontsize=11)
	plt.tick_params(axis='both', which='major', labelsize=12)
	# print(prods[n],physical_TTY/0.0036)

	
	# (Preferred units)

	if show_all_plots:
		plt.subplot(2, 2, 4)
		plt.plot(e_range,physical_TTY,label=prods[n], linestyle=line_style[n], linewidth=line_weight[n])
		plt.ylabel('Physical Thick Target Yield (MBq/uA h)')#,labelpad=7,fontsize=22)
		plt.xlabel('Beam Energy (MeV)')#,labelpad=7,fontsize=22)
		plt.yscale('log')
		plt.xlim(left=0)
		plt.legend(loc='best')
		# (Only for comparison against IAEA data)
		# To convert, 1 MBq/C = 0.0036 MBq/uAh


# print(physical_TTY[0])
# print(physical_TTY)
# print(physical_TTY[-1])


# plt.savefig('./Figures/PhD/Physical_Yields_Zr_log.pdf',dpi=600)
plt.show()



print('__________Qaims stuff:__________')
yield_86cum_from_50_40 = (22341-8846)*(1-np.exp(-np.log(2)/(14.74*3600)*3600))/(np.log(2)/(14.74*3600))*1e-6 #MBq
yield_87cum_from_50_40 = (4741-3095)*(1-np.exp(-np.log(2)/(79.89*3600)*3600))/(np.log(2)/(79.89*3600))*1e-6 #MBq
yield_88ind_from_50_40 = (111-45)*(1-np.exp(-np.log(2)/(106.626*24*3600)*3600))/(np.log(2)/(106.626*24*3600))*1e-6 #MBq
yield_90cum_from_50_40 = (15942-7704)*(1-np.exp(-np.log(2)/(64.046*3600)*3600))/(np.log(2)/(64.046*3600))*1e-6 #MBq
tot_yield_from_50_40 = yield_86cum_from_50_40+yield_87cum_from_50_40+yield_88ind_from_50_40+yield_90cum_from_50_40
print(f'86Y: {yield_86cum_from_50_40} MBq')
print(f'86Y: {yield_86cum_from_50_40/tot_yield_from_50_40*100} %')
print(f'87Y: {yield_87cum_from_50_40/tot_yield_from_50_40*100} %')
print(f'88Y: {yield_88ind_from_50_40/tot_yield_from_50_40*100} %')
print(f'90Y: {yield_90cum_from_50_40/tot_yield_from_50_40*100} %')

# __________Qaims stuff:__________
# 86Y: 47.45741460350038 MBq
# 86Y: 57.11440678979428 %
# 87Y: 7.10053808723608 %
# 88Y: 0.2859099145879196 %
# 90Y: 35.49914520838172 %
