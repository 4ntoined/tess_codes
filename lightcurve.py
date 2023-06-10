#Antoine
#take the photometry data and plot it

import numpy as np
import matplotlib.pyplot as plt

def filter99(data, date=np.array([]), mask=-99.999):
    dataxx = np.where( data == mask, np.nan, data )
    if len(date)>0:
        datexx = np.where( data == mask, np.nan, date )
        return (dataxx, datexx)
    else:
        return dataxx
def plot_lc(xdat,ydat,yerr,xlims='',ylims='',save=''):
    fig,ax = plt.subplots()
    fig.tight_layout()
    fig.set_size_inches(42,6)
    fig.dpi=100
    #
    xdat = xdat - 2459936.89779
    #
    ax.scatter(xdat,ydat,s=4.)
    ax.errorbar( xdat, ydat, yerr=yerr, ls='none', lw=1. )
    #
    if ylims:   ax.set_ylim(ylims[0],ylims[1])
    if xlims:
        ax.set_xlim(xlims[0]-2459936.89779,xlims[1]-2459936.89779)
        ax.set_xticks(np.arange(4,25,1))
    ax.set_xlabel('Days from observation start')
    ax.set_ylabel('TESS Magnitude')
    #
    plt.gca().invert_yaxis()
    #
    if save: plt.savefig(save,dpi=200,bbox_inches='tight')
    plt.show()
    return

ipath = '/run/user/1000/gvfs/sftp:host=copernicus.astro.umd.edu,user=antoine/'
        #'Volumes/TESS_05/AW_Data/Sector_60/S60_cam1_ccd4/Didymos/Didymos_comphot/comphot_mags.dat'
path1='Volumes/TESS_05/AW_Data/Sector_60/S60_cam1_ccd4/Didymos/Didymos_comphot/comphot_mags.dat'
path2 = 'Users/antoine/tess_output/Data/Sector_60/S60_cam1_ccd4/Didymos_files/Didymos/Didymos_comphot/comphot_mags.dat'
idata = np.loadtxt(ipath+path2,dtype=float,skiprows=2,usecols=(1,4,5,6,7,8,9,10,11,12,13,14,15,16,17))
#print(idata[:,5])

#idata25 = np.where( idata[:,5] == -99.999, np.nan, idata[:,5] )
#idate25 = np.where( idata[:,5] == -99.999, np.nan, idata[:,0] )
idata25, idate25 = filter99(idata[:,5], date=idata[:,0])
#ierr25 = np.where( idata[:,6] == -9.999, np.nan, idata[:,6])
ierr25 = filter99(idata[:,6],mask=-9.999)
#idate25-=2459936.89779
#plot_lc(idate25,idata25,ierr25,xlims=(2459940.,2459950.2),ylims=(12.5,13.5),save='s1.png')
#plot_lc(idate25,idata25,ierr25,xlims=(2459954.8,2459961.5),ylims=(12.75,14.25),save='s2.png')
plot_lc(idate25,idata25,ierr25,xlims=(2459940.,2459961.5),ylims=(12.5,14.25),save='s4.png')

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

