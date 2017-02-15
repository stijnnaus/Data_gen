# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 14:23:22 2017

@author: Stijn

Computes a running average of the seasonal cycle in methane.
This should be proportionate to the amount of CH4 in the atmosphere
"""

import numpy as np
from numpy import array
import matplotlib.pylab as plt
import os
import time
import sys

cwd = os.getcwd()

date,ch4_s,ch4_s_e,ch4_ns,ch4_ns_e = [],[],[],[],[]
f = open('ch4_mm_gl.txt','r') # global methane 
for line in f.readlines():
    if line[0] == '#': continue
    lin         = line.split()
    datei        = float(lin[2])
    ch4_si       = float(lin[3]) # normal ch4
    ch4_si_e     = float(lin[4])
    ch4_nsi      = float(lin[5]) # deseasonalized ch4
    ch4_nsi_e    = float(lin[6])
    date.append(datei)
    ch4_s.append(ch4_si)  ; ch4_s_e.append(ch4_si_e)
    ch4_ns.append(ch4_nsi); ch4_ns_e.append(ch4_nsi_e)
date = array(date)
ch4_s , ch4_s_e  = array(ch4_s) , array(ch4_s_e)
ch4_ns, ch4_ns_e = array(ch4_ns), array(ch4_ns_e)

f2 = open('ch4_gr_gl.txt','r')
yr,growth,growth_e = [],[],[]
for line in f2.readlines():
    if line[0] == '#': continue
    lin       = line.split()
    yri       = float(lin[0])
    growthi   = float(lin[1])
    growthi_e = float(lin[2])
    yr.append(yri)
    growth.append(growthi); growth_e.append(growthi_e)    
growth, growth_e = array(growth), array(growth_e)

sty, edy = int(date[0]), int(date[-1])
ny = edy-sty+1
# removing the trend
seas = (ch4_s - ch4_ns).reshape(ny,12)
seas_e = np.sqrt(ch4_s_e**2 + ch4_ns_e**2).reshape(ny,12)

amps,amps_e = np.zeros(ny), np.zeros(ny)
for i in range(ny):
    ipeak , ilow = np.argmax(seas[i]), np.argmin(seas[i])
    peak  , low  = seas[i][ipeak]    , seas[i][ilow]
    peak_e, low_e= seas_e[i][ipeak]  , seas_e[i][ilow]
    dif          = peak-low
    dif_e        = np.sqrt(peak_e**2+low_e**2)
    amps[i]      = dif
    amps_e[i]    = dif_e

amps_3y,   amps_3y_e   = np.zeros(ny-2), np.zeros(ny-2)
growth_3y, growth_3y_e = np.zeros(ny-2), np.zeros(ny-2)
for i in range(ny-2):
    amps_3y[i]     = np.mean(amps[i:i+3])
    amps_3y_e[i]   = np.sqrt( np.sum( (amps_e[i:i+3]/3)**2 ) )
    growth_3y[i]   = np.mean(growth[i:i+3])
    growth_3y_e[i] = np.sqrt( np.sum( (growth_e[i:i+3]/3)**2 ) )
    
    
r_1y = np.corrcoef(amps,growth)[0,1]
r_3y = np.corrcoef(amps_3y,growth_3y)[0,1]

plt.figure()
plt.title('The seasonal amplitude of methane')
plt.ylabel('Amplitude (ppb)')
plt.plot(range(sty,edy+1),amps,'-',color='red',label='Yearly')
plt.fill_between(range(sty,edy+1), amps-amps_e, amps+amps_e,color='lightcoral')
plt.plot(range(sty+1,edy),amps_3y,'k-',linewidth=4.,label='3yr running average')
plt.fill_between(range(sty+1,edy), amps_3y-amps_3y_e, amps_3y+amps_3y_e,color='gray')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Seasonal amplitude CH4.png')

plt.figure()
plt.title('Growth rate of methane')
plt.ylabel('Growth rate (ppb)')
plt.plot(range(sty,edy+1),growth,'-',linewidth=4.,color='red',label='Yearly')
plt.fill_between(range(sty,edy+1), growth-growth_e, growth+growth_e,color='lightcoral')
plt.plot(range(sty+1,edy),growth_3y,'k-',linewidth=4.,label='3yr running average')
plt.fill_between(range(sty+1,edy), growth_3y-growth_3y_e, growth_3y+growth_3y_e,color='gray')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Growth rate CH4.png')


plt.figure()
plt.plot(growth,amps,'o',label='r='+str(round(r_1y,3)))
plt.plot(growth_3y,amps_3y,'o',label='r='+str(round(r_3y,3)))
plt.legend(loc='best')







