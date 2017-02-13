# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 13:18:38 2016

@author: naus010
"""
import os
import numpy as np
from numpy import array
import matplotlib.pyplot as plt

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                             AGAGE
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def glob_av_agage(stations, yearly_avs, box_dist, box_weights, sty=1995, edy=2005):
    '''
    Computes the global yearly average from yearly averages per station.
    Takes into account box weights
    Returns global yearly averages
    '''
    nbox = len(box_weights)
    glob_av_y = [] # Global means
    for i,year in enumerate(range(sty,edy)):
        mean = True
        # Put concentration at each station in the correct box
        box_vals = [[] for ii in range(nbox)]
        for j,avi in enumerate(yearly_avs):
            val = None
            for ele in avi:
                if ele[0] == year: # Select year
                    val = ele[1]
            if val == None: continue
            
            stat = stations[j]
            box  = stat_chars[stat]['box'] # Number of the box
            box_vals[ box ].append(val)
        # Compute the mean in each box
        box_means = np.zeros(nbox)
        for j,vals in enumerate(box_vals):
            if vals == []: # Corrects for empty box
                #print 'In year',year,'box',j,'is empty'
                mean = False
            else:
                box_means[j] = sum(vals)/len(vals)
        # Compute the global average from box means
        if mean:
            glob_av_i = np.average(box_means, weights = box_weights)
            glob_av_y.append(np.array([year, glob_av_i]))
    return np.array(glob_av_y)

def read_mm_gage(fil,stations):
    ''' Reads the monthly mean GAGE data product per station'''
    mcfs_all, ch4s_all = [], []
    for i,stat in enumerate(stations): 
        print 'Reading GAGE data from',stat,'...'
        try:
            f = open(fil%(stat.upper()))
        except IOError:
            print 'Could not access GAGE data from',stat
            continue
        ch4s,mcfs = [],[]
        for line in f.readlines():
            lin = line.split()
            try:
                yy_d=float(lin[0])
            except ValueError:
                continue
            except IndexError:
                continue
            mm = int(lin[1])-1; yy = int(lin[2])
            mcf = float(lin[13]); mcf_e = float(lin[14])
            ch4 = float(lin[25]); ch4_e = float(lin[26])
            if mcf == 0.0: # Check for missing MCF data point
                mcfs.append([yy,mm,-999,-999])
            else:
                mcfs.append([yy,mm,mcf,mcf_e])
            if ch4 == 0.0: # Check for missing CH4 data point
                ch4s.append([yy,mm,-999,-999])
            else:
                ch4s.append([yy,mm,ch4,ch4_e])
        mcfs_all.append(mcfs)
        ch4s_all.append(ch4s)
    return mcfs_all, ch4s_all

def read_mm_agage(fil, stations):
    '''Reads the monthly mean AGAGE data product per station'''
    mcfs_all, ch4s_all = [], []
    for i,stat in enumerate(stations): 
        print 'Reading AGAGE data from',stat,'...'
        try:
            f = open(fil%(stat.upper()))
        except IOError:
            print 'Could not access AGAGE data from',stat
            continue
        ch4s,mcfs = [],[]
        for line in f.readlines():
            lin = line.split()
            try:
                yy_d=float(lin[0])
            except ValueError:
                continue
            except IndexError:
                continue
            mm = int(lin[1])-1; yy = int(lin[2])
            mcf = float(lin[10]); mcf_e = float(lin[11])
            ch4 = float(lin[22]); ch4_e = float(lin[23])
            if mcf == 0.0: # Check for missing MCF data point
                mcfs.append([yy,mm,-999,-999])
            else:
                mcfs.append([yy,mm,mcf,mcf_e])
            if ch4 == 0.0: # Check for missing CH4 data point
                ch4s.append([yy,mm,-999,-999])
            else:
                ch4s.append([yy,mm,ch4,ch4_e])
        mcfs_all.append(mcfs)
        ch4s_all.append(ch4s)
    return mcfs_all, ch4s_all

def read_glob_agage(fil):
    '''Reads the global monthly mean data product of AGAGE'''
    mcf,mcf_e,mcf_dates = [],[],[]
    ch4,ch4_e,ch4_dates = [],[],[]
    f = open(fil)
    for n,line in enumerate(f.readlines()[6:]):
        lin = line.split()
        yr_dec = float(lin[0])
        mcfi,mcfi_e = float(lin[9]),float(lin[10])
        ch4i,ch4i_e = float(lin[15]),float(lin[16])
        mcf.append(mcfi)
        mcf_e.append(mcfi_e)
        mcf_dates.append(yr_dec)
        if ch4i!=0.:
            ch4_dates.append(yr_dec)
            ch4.append(ch4i)
            ch4_e.append(ch4i_e)
    return np.array(mcf_dates),np.array(mcf),np.array(mcf_e),\
            np.array(ch4_dates), np.array(ch4), np.array(ch4_e)
    
def yav_ag(yrs_d, data):
    ''' 
    Computes yearly global means from global monthly means and decimal years 
    (only if 12 months of data are available for a given year)
    '''
    pyr = int(yrs_d[0])
    yrs = [] # Years with means available
    valsy = [] # Means of all years
    vals = [] # Values of one year
    for i,yr in enumerate(yrs_d):
        nyr = int(yr)
        val = data[i]
        if nyr == pyr:
            vals.append(val)
        elif len(vals)==12: # Only include if data is available for all months
            yrs.append(pyr)
            valsy.append(np.mean(vals))
            vals = [val]
        else:
            print 'In year',nyr,'only',len(vals),'months of data are available.'
            vals = [val]
        if i == len(yrs_d)-1: # BC
            if len(vals) == 12:
                yrs.append(pyr)
                valsy.append(np.mean(vals))
            else:
                print 'In year',yr,'only',len(vals),'months of data are available.'
        pyr = int(yr)
    return np.array(yrs)+.5, np.array(valsy)

def partition(data):
    '''
    Divide a dataset of format [[yr,mo,d1,d2..],[yr,mo+1,d1,d2..]] in years
    '''
    pyr = data[0][0]
    datap = [] # partitioned data
    datay = []
    for i,datam in enumerate(data):
        nyr = data[i][0]
        if nyr == pyr:
            datay.append(datam)
        else:
            datap.append(datay)
            datay = [datam]
        pyr = nyr
    datap.append(datay) # last year
    return datap

def fuse_ga_ag(data_ga, data_ag, stat_ga, stat_ag):
    '''
    Takes a gage and agage dataset with specified stations, and sticks the
    two together. Wherever they overlap, agage data is preferred.
    '''
    stat_all = np.unique(stat_ga+stat_ag)
    data_comb = []
    for stat in stat_all:
        if stat in stat_ga and stat in stat_ag: 
            iga = stat_ga.index(stat)
            iag = stat_ag.index(stat)
            dati_ga = data_ga[iga]
            dati_ag = data_ag[iag]
            dati_comb = []
            yr_st,mo_st = dati_ag[0][0],dati_ag[0][1] # starting point agage data
            for ii,dat in enumerate(dati_ga): # start with gage data
                yr,mo = dat[0],dat[1]
                if yr >= yr_st and mo >= mo_st: # stop if agage data is available
                    break
                dati_comb.append(dat)
                if ii<10:
                    print dati_comb
            dati_comb += dati_ag
        elif stat in stat_ga:
            iga = stat_ga.index(stat)
            dati_comb = data_ga[iga]
        elif stat in stat_ag:
            iag = stat_ag.index(stat)
            dati_comb = data_ag[iag]
        data_comb.append(partition(dati_comb))
    return data_comb, stat_all

# (A)GAGE station characteristics
cwd = os.getcwd()
dirc = cwd+'\\AGAGE data\\'
stat_ga = ['cgo', 'mhd', 'smo', 'rpb', 'org'] # GAGE stations
stat_ag = ['cgo', 'mhd', 'smo', 'rpb', 'thd'] # AGAGE stations
boxes = [[-90,-30],[-30,0],[0,30],[30,90]]
box_weights = [.25,.25,.25,.25]
stat_chars_ag, box_dist = read_station_chars(dirc+'AGAGE stations.txt', boxes)

# From global means
mcf_yr, mcf, mcf_e, ch4_yr, ch4, ch4_e = read_glob_agage(dirc+'global_mean_agage.txt')
mcf_raw_agage_yrs, mcf_raw_agage, mcf_raw_agage_e = mcf_yr, mcf, mcf_e
ch4_raw_agage_yrs, ch4_raw_agage, ch4_raw_agage_e = ch4_yr, ch4, ch4_e
ch4_ag_yrs, ch4_ag_glob = yav_ag(ch4_yr,ch4)
mcf_ag_yrs, mcf_ag_glob = yav_ag(mcf_yr,mcf)

# From station means
dirc_mm = os.path.join(dirc, 'monthly mean')
mcf_mm_ag, ch4_mm_ag = read_mm_agage(dirc_mm+'\\%s-gcmd.mon',stat_ag)
mcf_mm_ga, ch4_mm_ga = read_mm_gage(dirc_mm+'\\%s-gage.mon',stat_ga)
ch4_mm, stat_gaag = fuse_ga_ag(ch4_mm_ga, ch4_mm_ag, stat_ga, stat_ag)
mcf_mm, stat_gaag = fuse_ga_ag(mcf_mm_ga, mcf_mm_ag, stat_ga, stat_ag)
mcf_mm_ipl = interp_mm_all(mcf_mm,stations=stat_gaag)
ch4_mm_ipl = interp_mm_all(ch4_mm,stations=stat_gaag)
mcf_mm_yav = yav_mm_all(mcf_mm_ipl, stations=stat_gaag)
ch4_mm_yav = yav_mm_all(ch4_mm_ipl, stations=stat_gaag)

mcf_mm_glob = glob_av(stat_gaag, mcf_mm_yav, stat_chars_ag, box_weights, sty=1990, edy=2016)
ch4_mm_glob = glob_av(stat_gaag, ch4_mm_yav, stat_chars_ag, box_weights, sty=1990, edy=2016)

mcf_ag_yrs_mm = mcf_mm_glob[:,0]+.5
mcf_ag_glob_mm = mcf_mm_glob[:,1]
ch4_ag_yrs_mm = ch4_mm_glob[:,0]+.5
ch4_ag_glob_mm = ch4_mm_glob[:,1]


clrs = ['teal','gray','green','red','maroon','blue',\
        'steelblue','cyan','orange','lime','magenta']
box_markers = ['p', 'v', 's', '*', '^', 'd']
plt.figure(figsize=(15,15))
plt.title('The methane concentration, as derived from monthly mean AGAGE data.\n\
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
    stat = stat_gaag[ist]
    ibox = stat_chars_ag[stat]['box']
    plt.plot(mdates, filt, '-',color=clrs[ist])
    plt.plot(yrs+.5, yavs, box_markers[ibox], color=clrs[ist], label='station '+str(stat)+' ('+str(mis)+' missing data points)')
plt.plot(ch4_ag_yrs_mm+.5, ch4_ag_glob_mm, 'ko-', label='Global mean')
plt.tight_layout()
plt.legend(loc='best')
plt.savefig('CH4_from_AGAGE_mm_full')



















