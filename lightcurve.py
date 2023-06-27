#Antoine
#take the photometry data and plot it

import numpy as np
import matplotlib.pyplot as plt

def filter99(data, date=np.array([]), mask=-99.999):
    dataxx = np.where( data == mask, np.nan, data )
    datayy = np.where( dataxx == 15.440, np.nan, data)
    if len(date)>0:
        dateyy = np.where( datayy, np.nan, date )
        #dateyy = np.where( dat == 15.440, np.nan, date)
        return (datayy, dateyy)
    else:
        return datayy
def plot_lc(xdat,ydat,yerr,xlims='',ylims='',save='',plot_delta=''):
    global delx, delt,autokm
    if plot_delta:  fig,(az,ax) = plt.subplots(2,1,sharex=True,height_ratios=[2,6])
    else:           fig,ax = plt.subplots()
    fig.tight_layout()
    fig.set_size_inches(42,6.5)
    fig.dpi=100
    #
    xdat = xdat - 2459936.89779
    delx = delx - 2459936.89779
    delt = delt / autokm
    #
    ax.scatter(xdat,ydat,s=4.)
    #ax.scatter(xdat,idata25,s=4.,color='r' )
    ax.errorbar( xdat, ydat, yerr=yerr, ls='none', lw=1. )
    if plot_delta: az.plot(delx, delt, lw=1.,color='k')
    #
    ax.tick_params(axis='both',length=15.,width=5.,labelsize='larger')
    if ylims:   ax.set_ylim(ylims[0],ylims[1])
    if xlims:
        ax.set_xlim(xlims[0]-2459936.89779,xlims[1]-2459936.89779)
        ax.set_xticks(np.arange(4,25,1))
    ax.set_xlabel('Days from observation start')
    ax.set_ylabel('TESS Magnitude')
    ax.grid(axis='both')
    if plot_delta:
        az.set_ylabel('TESS-Didymos distance [AU]')
        #az.grid(axis='y')
        az.grid(axis='both')
        az.tick_params(axis='y',length=15.,width=5.,labelsize='larger')
    #
    ax.invert_yaxis()
    #
    if save: plt.savefig(save,dpi=300,bbox_inches='tight')
    plt.show()
    return
def time_binning(data, time, period=3.):
    dperi = period / 24.
    pp = dperi/2.
    ans = []
    for i in range(len(time)):
        b1 = time[i] - pp
        b2 = time[i] + pp
        #print(time,b1, b2)
        bounded = np.argwhere(np.logical_and(time > b1, time < b2))
        #
        #print()
        #print(bounded)
        #if len(bounded)==0:
        #ans.append(np.nan)
        #else:
        if period==0.:
            ans.append(np.array([i]))
        else:
            ans.append(np.squeeze(bounded))
        #print(bounded)
        #print(bounded)
    #b1 = time - pp
    #b2 = time + pp
    #print(data.shape, time.shape, type(period), type(b1), type(b2))
    #bounded = np.argwhere( np.logical_and( time > b1, time < b2))
    #print(bounded.shape)
    return ans
def time_smooth(bins, data, time):
    #binned_data is a list of bins, indexes
    ans = []
    for i in range(len(bins)):
        binn = bins[i]
        #print(binn)
        #if type(binn) == float:
        #ans.append(np.nan)
        #else:
        print(type(binn))
        print(binn)
        print(binn.size)
        if binn.size == 1:
            ans.append(data[binn])
        else:
            #print(binn)
            if len(binn) != 0:
                dato = data[ binn[0]:binn[-1]+1 ]
                n_nan = np.count_nonzero( np.isnan( dato ))
                nan_correct  = np.size(dato) / (np.size(dato) - n_nan)
                #print(dato)
                binmean = np.nanmean( dato )
                binfix = binmean * nan_correct
                ans.append(binfix)
            else:
                ans.append(np.nan)
            pass
        pass
    return ans


autokm = 1.496e8 #km in an au
ipath = '/run/user/1000/gvfs/sftp:host=copernicus.astro.umd.edu,user=antoine/'
        #'Volumes/TESS_05/AW_Data/Sector_60/S60_cam1_ccd4/Didymos/Didymos_comphot/comphot_mags.dat'
path1='Volumes/TESS_05/AW_Data/Sector_60/S60_cam1_ccd4/Didymos/Didymos_comphot/comphot_mags.dat'
path2 = 'Users/antoine/tess_output/Data/Sector_60/S60_cam1_ccd4/Didymos_files/Didymos/Didymos_comphot/comphot_mags.dat'
path3 = '/home/antoine/Desktop/codespace_tess/comphot_mags.dat'
jpath = '/home/antoine/Desktop/codespace_tess/didy_delta_eph.txt'
#idata = np.loadtxt(ipath+path2,dtype=float,skiprows=2,usecols=(1,4,5,6,7,8,9,10,11,12,13,14,15,16,17))
idata = np.loadtxt(path3,dtype=float,skiprows=2,usecols=(1,4,5,6,7,8,9,10,11,12,13,14,15,16,17))
idata2 = np.loadtxt(jpath,dtype=float,usecols=(2,11))
#print(idata[:,5])
delt = idata2[:,1].copy()
delx= idata2[:,0].copy()
#plt.plot( idata2[:,0], idata2[:,1] )
#plt.show()

#idata25 = np.where( idata[:,5] == -99.999, np.nan, idata[:,5] )
#idate25 = np.where( idata[:,5] == -99.999, np.nan, idata[:,0] )
#idata25, idate25 = filter99(idata[:,5], date=idata[:,0])
mask1 = np.logical_or( idata[:,5] == -99.999, idata[:,5] == 15.440 )
mask2 = idata[:,6] == -9.999
idata25 = np.where( mask1, np.nan, idata[:,5]  )
idate25 = np.where( mask1, np.nan, idata[:,0]  )

ierr25 = np.where( mask2, np.nan, idata[:,6]  )

print(idata25)
#ierr25 = filter99(idata[:,6],mask=-9.999)
#ierr25 = np.where(idata[:,6] == -9.999, np.nan, idata[:,6])
#print(ierr25-idata[:,6])
print(np.count_nonzero(idata25==15.440))

eee = time_binning(idata25, idate25, period=4.)
#print(eee)
fff = time_smooth(eee, idata25, idate25)
#print(fff)
smooth = np.array(fff,dtype=float)
#print(eee[10])
#print(smooth)
print(len(smooth),len(idate25),len(idata25))
#idate25-=2459936.89779
#plot_lc(idate25,idata25,ierr25,xlims=(2459940.,2459950.2),ylims=(12.5,13.5),save='s1.png')
#plot_lc(idate25,idata25,ierr25,xlims=(2459954.8,2459961.5),ylims=(12.75,14.25),save='s2.png')
plot_lc(idate25,smooth,ierr25,xlims=(2459940.,2459961.5),ylims=(12.5,14.25),plot_delta=True,save='s10.png')
#plot_lc(idate25,idata25,ierr25,xlims=(2459940.,2459961.5),ylims=(12.5,14.25),plot_delta=True)
#plot_lc(idate25,idata25,ierr25,xlims=(2459940.,2459961.5),ylims=(12.5,14.25),plot_delta=True,save='s8.png')

#fig,ax = plt.subplots()
#fig.figsize = (8,6)
#fig.dpi=140
#fig.tight_layout()
#
#ax.scatter(idate25,idata25,s=4.)
#ax.errorbar( idate25, idata25, yerr=ierr25, ls='none' )
#
#plt.gca().invert_yaxis()
#plt.show()

