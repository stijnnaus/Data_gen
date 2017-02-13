# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 11:39:57 2016

@author: Stijn
"""
import numpy as np
from numpy import array

def read_station_chars(fil,boxes):
    '''
    Reads in the station characteristics and based on their latitude, it 
    distributes them between the specified latitudinal box distribution.
    Returns box distribution and a dict with characteristics of each station.
    '''
    f = open(fil, mode = 'r')
    box_dist = [[] for n in range(len(boxes))]
    stat_chars = {}
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for line in f.readlines():
        if line[0] not in alphabet: continue
        lin = line.split()
        code = lin[0].lower()
        lat  = float(lin[1])
        lon  = float(lin[2])
        alt  = float(lin[3])
        for j,box in enumerate(boxes):
            if lat >= box[0] and lat < box[1]:
                box_dist[j].append(code)
                box_no = j
        dic = {}
        dic['lat'] = lat; dic['lon'] = lon; dic['ele'] = alt; dic['box'] = box_no
        stat_chars[code] = dic
    return stat_chars, box_dist

def glob_av(stations, yearly_avs, stat_chars, box_weights, sty=1995, edy=2005, bm=False):
    '''
    Computes the global yearly average from yearly averages per station.
    Takes into account box weights
    Returns global yearly averages
    '''
    nbox = len(box_weights)
    glob_av_y = [] # Global means
    box_meanss = []
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
        box_meanss.append(box_means)
        # Compute the global average from box means
        if mean:
            glob_av_i = np.average(box_means, weights = box_weights)
            glob_av_y.append(array([year, glob_av_i]))
    if bm:
        return array(box_meanss), array(glob_av_y)
    return array(glob_av_y)

def yav_mm_all(data, stations=[]):    
    ''' 
    Computes the yearly average from monthly mean data for all stations.
    Only computes the yearly average if data for all months are present.
    '''
    yavs_all = []
    for i,datast in enumerate(data): # Iteratively average for all stations
        if len(stations) == len(data):
            yavs_st = yav_mm_stat(datast, stations[i])
        else: 
            yavs_st = yav_mm_stat(datast)
        yavs_all.append(yavs_st)
    return np.array(yavs_all)
        

def yav_mm_stat(data, station=None):
    ''' 
    Computes the yearly average from monthly mean data for 1 station.
    Only computes the yearly average if data for all months are present.
    '''
    means = []
    for iy, datay in enumerate(data):
        vals = datay[:,2]
        y = datay[0][0]
        if -999 not in vals: # Only if no month is missing
            means.append(np.array([y, np.mean(vals)]))
        elif station!=None:
            print 'No yearly mean for station',station,'in year',y
    return np.array(means)

def interp_mm_all(data, tech='yearly', tol=1, stations = []):
    ''' 
    Interpolates the data for all stations
    tol: The number of years the function will look for values
    tech: The interpolation technique. This is either interpolation between
    the next and prev year (yearly, default) or between the next and prev month
    (monthly)->NOT IMPLEMENTED YET!.
    Output: A round array, of shape nyearx12. Gaps that could not be filled are 
    indicated by -999.
    '''
    data_ipl = []
    for i,datai in enumerate(data):
        if len(stations) == len(data):
            print 'Interpolating station',stations[i]
        datai_ipl = interp_mm_stat(datai, tech='yearly', tol='1')
        data_ipl.append(datai_ipl)
    return data_ipl
        

def interp_mm_stat(data, tech='yearly', tol=1):
    '''
    Fills in missing data through linear interpolation for one station.
    data: A list of [ [year,month,value], ...], [year+1,month,value], ...] ] 
    with some months/years missing.
    '''
    try:
        assert tech == 'yearly'
    except AssertionError:
        print 'As of yet the only tech implemented is yearly. So use yearly.'
    ny = len(data)
    data_filled = fill_gaps(data)
    data_interp = np.copy(data_filled)
    for iy,datay in enumerate(data_filled):
        for im,datam in enumerate(datay):
            val = datam[2]
            if val != -999: continue
            if iy == 0 or iy+1 == ny: continue # Boundary conditions
            # The interpolation
            ipr,prev = 0,-999; ine,nex = 0,-999
            while prev == -999 and (iy-ipr-1) >= 0 and ipr <= tol:
                ipr+=1; prev = data_filled[iy-ipr][im][2]
            while nex  == -999 and (iy+ine+1) < ny and ine <= tol:
                ine+=1; nex  = data_filled[iy+ine][im][2]
            if prev == -999 or nex == -999: continue # No interpolation possible
            #print 'Interpolated at year', datam[0], 'and month', datam[1]
            data_interp[iy][im][2] = np.interp(iy, [iy-ipr, iy+ine], [prev, nex])
    return data_interp

                        
def fill_gaps(data):
    ''' Makes an array round by filling gaps with -999 values '''
    # Make the array round by filling in the missing data with -999
    ny = len(data)
    ndata = len(data[0][0]) # Number of data elements (excepting yr,mo)
    filled = np.zeros((ny, 12, ndata))
    for iy,datay in enumerate(data):
        y = datay[0][0]
        for im in range(12):
            filled[iy][im][0] = y
            filled[iy][im][1] = im
            val = -999 # True if y,m is present in the data
            for i in range(len(datay)):
                mi = datay[i][1]
                if mi == im: # Check if month is found in data
                    val = datay[i][2:ndata]
            filled[iy][im][2:ndata] = val
    return filled