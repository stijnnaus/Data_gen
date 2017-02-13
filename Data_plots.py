# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:37:39 2016

@author: naus010
"""

'''
First run Read_NOAA_Data and Read_CH4_Data.
Then run this file to make the corresponding plots.
'''
import os
from matplotlib.ticker import MultipleLocator as ML

def write_data(years, data, filename, header=None):
    '''
    Write yearly averages to a text file, with one column the years, the other
    the concentrations and an optional header.
    '''
    fileloc = os.path.join(os.getcwd(), 'Data output',filename)
    f = open(fileloc, 'w')
    if header!=None:
        f.write(header+'\n')
    for i,yr in enumerate(years):
        val = data[i]
        f.write('%i'%(int(yr)))
        f.write('\t')
        f.write('%.3f'%(val))
        f.write('\n')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                              AGAGE PLOTS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fig = plt.figure(figsize = (10,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
ax1.set_title('Global mean monthly AGAGE data')
ax1.set_ylabel(r'CH$_4$ (ppb)')
ax2.set_ylabel('MCF (ppt)')
ax1.errorbar(ch4_raw_agage_yrs, ch4_raw_agage, yerr = ch4_raw_agage_e, color = 'red', fmt='o', label = r'CH$_4$' )
ax2.errorbar(mcf_raw_agage_yrs, mcf_raw_agage, yerr = mcf_raw_agage_e, color = 'blue', fmt='o', label = 'MCF' )
ax1.legend(loc = 'lower left')
ax2.legend(loc = 'upper right')
plt.savefig('AGAGE_glob_mean_MCF_CH4')

fig3 = plt.figure(figsize = (10,20))
ax1 = fig3.add_subplot(211)
ax2 = fig3.add_subplot(212)
ax1.set_title('Global yearly averages')
ax2.set_xlabel('Year')
ax1.set_ylabel(r'CH$_4$ (ppb)')
ax2.set_ylabel('MCF (ppt)')
ax1.plot(ch4_ag_yrs, ch4_ag_glob, 'ro', label = r'AGAGE (from global monthly means)' )
ax1.plot(ch4_ag_yrs_mm, ch4_ag_glob_mm, 'bo', label = r'AGAGE (from station monthly means)') # chhh4 is from OH_GTB
ax1.plot(ch4_noaa_yrs, ch4_noaa_glob, 'go', label='NOAA (from station monthly means)')
ax2.plot(mcf_ag_yrs, mcf_ag_glob, 'ro', label = 'AGAGE (from global monthly means)' )
ax2.plot(mcf_ag_yrs_mm, mcf_ag_glob_mm, 'bo', label = 'AGAGE (from station monthly means)' ) # mcfff is from OH_GTB
ax2.plot(mcf_noaa_yrs, mcf_noaa_glob, 'go', label = 'NOAA (from event-wise data)') # glob_av_mon is from Read_NOAA_Data
ax1.legend(loc = 'best')
ax2.legend(loc = 'best')
plt.savefig('global_mean_comparison')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                              MCF PLOTS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#x = np.linspace(1992,2016)
#y = np.linspace(-90,90)
#colors = ['blue','red','lightgray','maroon','darkgreen','lime','dimgray','steelblue', \
#          'navy','brown', 'cyan','magenta','indigo','orange']
#nx = len(x)
#fig = plt.figure(figsize = (13,13))
#ax1 = fig.add_subplot(111)
#ax2 = ax1.twiny()
#ax1.set_title('A plot of the availability of NOAA data for MCF. \n \
#The position of the station name points out the location of the station (top x-axis). \n \
#The line depicts the time period that the data from that station covers (bottom x-axis)', y=1.05)
#for i,stat in enumerate(stations):
#    lati,loni,styi,edyi,boxi = stat_noaa[stat]
#    peri = [styi, edyi]
#    ax1.plot( peri, [lati]*2, 'o-', color = colors[i%2] )
#    
#    if stat == 'kum': # because it is so close to smo
#        ax2.text( loni, lati-3., stat, fontweight='bold', color = colors[i%2] )
#    else:
#        ax2.text( loni, lati   , stat, fontweight='bold', color = colors[i%2] )
#ax1.plot(x, nx*[90], 'k--')
#ax1.plot(x, nx*[60], 'k--')
#ax1.plot(x, nx*[30], 'k--')
#ax1.plot(x, nx*[0], 'k--')
#ax1.plot(x, nx*[-30], 'k--')
#ax1.plot(x, nx*[-60], 'k--')
#ax1.plot(x, nx*[-90], 'k--')
#ax1.plot([styear]*nx, y, 'k--')
#ax1.plot([edyear]*nx, y, 'k--')
#ax1.axis([1991,2017.5,-95,95])
#ax2.axis([0.,360.,-95,95])
#ax1.grid(True)
#ax1.set_xlabel('Data period covered (yr)')
#ax2.set_xlabel('Longitude measurement station (deg east)')
#ax1.set_ylabel('Latitude measurement station')
#plt.tight_layout()
#plt.savefig('NOAA_MCF_data_availability')
#
#fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True, figsize = (10,25))
#ax1.set_title('MCF concentrations per NOAA station')
#ax2.set_title('Yearly averages per NOAA station and global AGAGE average')
#ax1.set_ylabel('MCF (ppb)')
#ax2.set_ylabel('MCF (ppb)')
#ax3.set_ylabel('MCF (ppb)')
#ax3.set_xlabel('Year')
#ax1.axis([1986.,2017.,0.,180.])
##ax2.errorbar(np.array(range(1988,1988+len(mcfff)))+.5,mcfff, yerr = mcf_obs_e,fmt='o',color='black',label = '(A)GAGE')
#for i,stat in enumerate(stat_gage): # GAGE/AGAGE stations
#    for j in range(2): # GAGE and AGAGE
#        if j == 0:
#            ax2.plot(ts_agage_dec[i][j], mcfs_agage[i][j], '--', color = colors[i],label = stat+' from AGAGE')
#        else:
#            ax2.plot(ts_agage_dec[i][j], mcfs_agage[i][j], '--', color = colors[i])
#for i,yavi in enumerate(yavs_mcf):
#    # Data
#    ax1.plot(mcf_obs[i][:,0], mcf_obs[i][:,1], '.', color = colors[i], label = stations[i])
#    # Yearly averages per station
#    ax2.plot(np.array(yavi[:,0])+0.5, yavi[:,1], 'o', color = colors[i], label = stations[i])
#ax2.legend( loc='best', ncol = 2, numpoints = 1 )
## Global ,yearly averages
#ax3.set_title('Global, yearly averages for the NOAA and the AGAGE data. \n \
#For NOAA results of two different box weighting procedures are shown.')
#ax3.plot( mcf_noaa_yrs, mcf_noaa_glob, 'o', label = 'NOAA' )
##ax3.errorbar(np.array(range(1988,1988+len(mcfff)))+.5,mcfff, yerr = mcf_obs_e,fmt='o',color='gray',label = '(A)GAGE my average')
#ax3.plot( mcf_ag_yrs, mcf_ag_glob, 'o', label = '(A)GAGE their average' )
#ax3.legend(loc='best', numpoints = 1)
#plt.savefig('Yearly_glob_avs_AGAGE-NOAA')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                       WRITING DATA
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

write_data(ch4_noaa_yrs, ch4_noaa_glob, 'ch4_noaa_glob.txt', \
    header='# Global yearly mean methane concentrations from NOAA data\n# Year\tCH4(ppb)')
write_data(ch4_ag_yrs, ch4_ag_glob, 'ch4_agage_glob.txt', \
    header='# Global yearly mean methane concentrations from AGAGE data\n# Year\tCH4(ppb)')
write_data(mcf_noaa_yrs, mcf_noaa_glob, 'mcf_noaa_glob.txt', \
    header='# Global yearly mean MCF concentrations from NOAA data\n# Year\tMCF(ppt)')
write_data(mcf_ag_yrs, mcf_ag_glob, 'mcf_agage_glob.txt', \
    header='# Global yearly mean MCF concentrations from AGAGE data\n# Year\tMCF(ppt)')









plt.plot(mcf_noaa_yrs, mcf_noaa_glob,'o')




ch4_e=[3.]*len(ch4_ag_yrs)
fig = plt.figure(figsize = (10,10))
ax1 = fig.add_subplot(111)
ax1.set_title('Global mean methane from AGAGE data')
ax1.set_ylabel(r'CH$_4$ (ppb)')
ax1.plot(ch4_ag_yrs, ch4_ag_glob,'g-', label = r'CH$_4$',linewidth=3.0)
plt.savefig('AGAGE_CH4.png')

plt.rcParams.update({'font.size': 18})
major_ticks = ML(5)
minor_ticks = ML(1)
fig = plt.figure(figsize=(12,7))
ax1 = fig.add_subplot(111)
ax1.set_ylabel('MCF (ppt)')
ax1.grid(True)
ax1.get_yaxis().set_tick_params(which='both', direction='out')
ax1.get_xaxis().set_tick_params(which='both', direction='out')
ax1.xaxis.set_major_locator(major_ticks)
ax1.xaxis.set_minor_locator(minor_ticks)
ax1.tick_params(axis='x',which='minor')
ax1.tick_params(axis='y',which='minor')
ax1.plot(mcf_ag_yrs, mcf_ag_glob,'-', color='gray', linewidth=2.0)
ax1.plot(mcf_ag_yrs, mcf_ag_glob,'o',color='steelblue',markersize=11)
plt.savefig('AGAGE_MCF.png')









