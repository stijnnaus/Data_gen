# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:38:20 2016

Reads all NOAA data: Both CH4 and MCF.

@author: naus010
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from math import *
import os

def calc_yav_event(data, stations):
    '''
    Compute yearly averages per station. 
    The data is eventwise, which is taken into account by weighing the average 
    of two subsequent measurements by the measurement interval.
    Returns the yearly averages per station.
    '''
    dts_all = []
    avs_all = []
    for i,stat in enumerate(stations):
        datai = np.array(data[i])
        sty,edy = int(datai[0][0]), int(datai[-1][0])
        print 'Averaging:',stat,'from',sty,'to',edy
        
        ts,mcf,mcfe = datai[:,0],datai[:,1],datai[:,2]
        av,avs = 0.,[]
        dts = []
        for t in range(len(ts)):
            ye_new, ye_prv = int( ts[t] ), int( ts[t-1] )
            if t == 0: # BC 1 (first year)
                dt = ts[0]- ye_new
                av += dt * mcf[0]
                dts.append(dt)
            elif ye_new != ye_prv: # Transition between years
                virt = np.interp( ye_new, [ts[t-1],ts[t]], [mcf[t-1],mcf[t]] )
                av += ( virt + mcf[t-1] ) / 2 * ( ye_new - ts[t-1] )
                avs.append( [ ye_prv, av ] )
                dts.append(ts[t-1]-ts[t])
                av = ( mcf[t] + virt ) / 2 * ( ts[t] - ye_new )
            else: # Normal
                dmcf = (mcf[t] + mcf[t-1]) / 2
                dt   = ts[t]  - ts[t-1]
                dts.append(dt)
                av += dmcf * dt
                
            if t+1 == len(ts): # BC 2 (last year)
                dt = ye_new - ts[t] + 1
                av += dt * mcf[t]
                dts.append(dt)
                avs.append( [ int(ts[t]), av ] )
        dts_all.append(np.array(dts))
        avs_all.append(np.array(avs))
    return np.array(avs_all),dts_all

def read_noaa_event(fil):
    '''
    Reads the event-wise NOAA flask data (MCF).
    Returns station names and data per station.
    '''
    # Reading the measurements
    f = open(fil, mode = 'r')
    stations = []
    data     = []
    i = 0 # measurement number
    for m,line in enumerate(f.readlines()):
        if line[0] == '#': continue
        lin = line.split()
        name = lin[0]
        timei = float(lin[1])
        if timei < styear or timei >= edyear: continue
        coni = float(lin[6]) # MCF concentration
        coni_e = float(lin[7]) # Error
        
        if i == 0:
            stations.append(name)
            datai = []
        elif name!=stations[-1]:
            stations.append(name)
            data.append( np.array(datai) )
            datai = []
            
        datai.append( [timei,coni,coni_e] )
        i+=1
    data.append( np.array(datai) )
    
    return stations, np.array(data)
    
def read_noaa_mm(fil_start, stations):    
    ''' 
    Returns monthly mean data per station.
    '''
    datas = []
    for ist,stat in enumerate(stations):
        data, datai = [], []
        f = open(fil_start+stat+'.txt', 'r')
        yrs = []
        for i,line in enumerate(f.readlines()):
            if line[0] == '#': continue
            lin = line.split()
            yr = int(lin[1])
            if len(yrs) == 0: yrs.append(yr)
            elif yr != yrs[-1]: yrs.append(yr); data.append(datai); datai = []
            mo = int(lin[2])-1
            con = float(lin[3])
            datai.append([yr,mo,con])
        datas.append( np.array(data) )
    return datas

    
cwd = os.getcwd()
# NOAA station characteristics
# Latitudinal distribution of the boxes:
boxes = [ [-90., -30.], [-30.,0.], [0.,30.], [30.,60.], [60.,90.] ]
nbox = len(boxes)
styear, edyear = 1992, 2016
mcf_yrs = range(styear, edyear)
ar1, ar2, ar3  = 0.5, 0.5*sqrt(3) - 0.5, 1 - 0.5*sqrt(3)
tar1,tar2,tar3 = ar1*1.5, ar2*1.25, ar3*1.
w_area = [ar3+ar2, ar1, ar1, ar2, ar3 ] # Weights per box based on area per box
w_trop = [tar2+tar2,tar1,tar1,tar2,tar3] # Weights per box based on area per box times tropopause height
stat_chars, box_dist = read_station_chars(cwd+'\\NOAA data\\NOAA_stations_ALL2.txt',boxes) # All possible stations

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                           Methyl Chloroform
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

stat_noaa, mcf_obs = read_noaa_event(cwd+'\\NOAA data\\MCF_NOAA_flask.txt')
yavs_mcf,dts_mcf = calc_yav_event( mcf_obs, stat_noaa )
glob_av_mcf1 = glob_av( stat_noaa, yavs_mcf, stat_chars, w_area, sty=styear, edy=edyear )
glob_av_mcf2 = glob_av( stat_noaa, yavs_mcf, stat_chars, w_trop, sty=styear, edy=edyear )

# Select the stations that Montzka (2011) used
stat_mon = np.array(['spo','cgo','smo','mlo','kum','nwr','lef','brw','alt'])
nmon = len(stat_mon)
stat_srt = []
yavs_mon  = []
for i,stat in enumerate(stat_noaa):
    if stat in stat_mon:
        stat_srt.append(stat)
        yavs_mon.append(yavs_mcf[i])
stat_mon = stat_srt
yavs_mon = np.array(yavs_mon)
glob_av_mon = glob_av( stat_mon, yavs_mon, stat_chars, w_area, sty=styear, edy=edyear )
mcf_noaa_yrs = glob_av_mon[:,0]+.5
mcf_noaa_glob = glob_av_mon[:,1]

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                           Methane
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

stym_mm, edym_mm = 1983, 2015 # Year in which the first measurements start
sty_mm, edy_mm = 1985, 2015 # Time interval of interest
nyrm_mm = edym_mm-stym_mm+1
nyr_mm = edy_mm-sty_mm+1
stats_mm = np.array(['spo', 'smo', 'cgo', 'alt', 'brw', 'kum', 'mlo', 'nwr', 'asc','mhd','thd'])
boxes_mm = [stat_chars[stat]['box'] for stat in stats_mm]
stats_mm = stats_mm[np.lexsort((stats_mm,boxes_mm))][::-1]
nst_mm = len(stats_mm)
fil_start = os.path.join(cwd,'NOAA data','monthly mean','ch4_mm_')
ch4_mm_raw = read_noaa_mm(fil_start, stats_mm)
ch4_mm_ipl = interp_mm_all(ch4_mm_raw, stations=stats_mm)
ch4_mm_yav = yav_mm_all(ch4_mm_ipl, stations=stats_mm)
ch4_mm_global = glob_av(stats_mm, ch4_mm_yav, stat_chars, w_area, sty=sty_mm, edy=edy_mm)

clrs = ['teal','gray','green','red','maroon','blue',\
        'steelblue','cyan','orange','lime','magenta']
box_markers = ['p', 'v', 's', '*', '^', 'd']
plt.figure(figsize=(15,15))
plt.title('The methane concentration, as derived from monthly mean NOAA data.\n\
    The monthly mean per station, the yearly mean per station, and the global yearly average\n\
    are given by the solid lines, by the colored markers, and by the black markers respectively.\n\
    The marker shape of the yearly averages indicates the box in which they were placed.')
plt.xlabel('Year')
plt.ylabel('CH4 (ppb)')
for ist,data in enumerate(ch4_mm_ipl):
    # Monthly means
    mis = 0
    l1,l2,l3 = data.shape
    filt, mdates = [], []
    dataf = np.reshape(data,(l1*l2,l3))
    for i,val in enumerate(dataf[:,2]): # Don't plot missing measurements
        y, m = dataf[i][0], dataf[i][1]
        if val != -999: 
            mdates.append(y+m/12.)
            filt.append(val)
        else: mis+=1
    # Yearly means
    yrs = ch4_mm_yav[ist][:,0]
    yavs = ch4_mm_yav[ist][:,1]
    stat = stats_mm[ist]
    ibox = stat_chars[stat]['box']
    plt.plot(mdates, filt, '-',color=clrs[ist])
    plt.plot(yrs+.5, yavs, box_markers[ibox], color=clrs[ist], label='station '+str(stat)+' ('+str(mis)+' missing data points)')
plt.plot(ch4_mm_global[:,0]+.5, ch4_mm_global[:,1], 'ko-', label='Global mean')
plt.tight_layout()
plt.legend(loc='best')
plt.savefig('CH4_from_NOAA_mm_full')

ch4_noaa_glob = ch4_mm_global[:,1]
ch4_noaa_yrs = ch4_mm_global[:,0]+.5




# Computing an approximate measurement frequency for MCF
nst = len(stat_noaa)
lens = array([len(mcf_obs[i]) for i in range(nst)])
yrs = array([24.,23.,24.,20.,20.,19.,18.,24.,24.,17.,24.,23.,12.,20.]) # number of years that measurements are taken at each station
measf = np.transpose((stat_noaa,lens/yrs)) # measurement frequency per station


















    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    