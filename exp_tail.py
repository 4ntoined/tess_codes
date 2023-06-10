#going to measure the tail brightness
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits

def square_ap(center,radius=1):
    """
    center: 2-tuple of (x,y) on frame
    radius: number of pixels from center to include
    will return a pair of splice values that will target a 2*rad+1
        square box centered on center
    """
    x,y = center
    xpart = (x-radius,x+radius+1)
    ypart = (y-radius,y+radius+1)
    return (xpart, ypart)
def photo(frame,radius=1):
    """
    will take in [frame], measure the brightness in central 9 pixels
    and then track and measure the tail brightness
    frame : 2d array where comet is centered
    """
    ny, nx = frame.shape
    centerx, centery = (nx-1)//2, (ny-1)//2
    #center
    aper = square_ap( (centerx,centery), radius=radius )
    center9 = frame.copy()[ aper[1][0]:aper[1][1],aper[0][0]:aper[0][1] ]
    boy = photo2(center9)
    #nan_catch = np.count_nonzero( np.isnan( center9 ))
    #sum_center = np.nansum( center_9 )
    #center_phot = (sum_center, nan_catch)
    #for these 401x401 images, center is index 200, for radius 1, we can go as far
    #as center,index=399, or size-2, or size-(radius+1), for 9 and radius 2, middle is 4, end is 6!
    #all
    liss = []
    for i in range( centerx, nx-(radius+1) ):
        #thisone = ( centery, i  )
        aper = square_ap( (i, centery), radius=radius )
        data9 = frame.copy()[ aper[1][0]:aper[1][1],aper[0][0]:aper[0][1] ]
        photom = photo2(data9)
        liss.append(photom)
    return liss
def photo2(aperture):
    nsize = aperture.size
    summ = np.nansum( aperture )
    nan_catch = np.count_nonzero(np.isnan( aperture ))
    correction = nsize / ( nsize - nan_catch )
    corrected = summ * correction
    return (summ, nan_catch, corrected)
def get_photo(data):
    #where data is 1d array of tail photometry
    #xax = np.arange( 200,399, 1,dtype=int)
    ns = np.size(data) 
    yax = data
    xax = np.arange(ns,dtype=int)
    #plot_photo(xax,yax)
    return (xax,yax)
def plot_photo(datax, datay, title=''):
    """
    datax, datay: 1d arrays of floats to have plotted
    """
    fig, ax = plt.subplots()
    fig.set_size_inches(6,6)
    fig.dpi=140
    #
    ax.plot(datax, datay)
    ax.set_yscale('log')
    ax.set_xlim(0,100)
    #ax.set_ylim(bottom=1.)
    #
    if title: ax.set_title(title)
    ax.set_xlabel('pixels from center')
    ax.set_ylabel('profile')
    #plt.savefig(savehere,dpi=200,bbox_inches='tight')
    plt.show()
    return
def photo_all(files,radius=1):
    saved = []
    for i in range(len(files)):
        ifile = fits.open(top_data_directory+files[i])
        dati = ifile[0].data
        phot1 = photo( dati, radius=radius )
        #print(phot1)
        phot2 = [ i[2] for i in phot1 ]
        phot3 = np.array(phot2,dtype=float)
        xx,yy = get_photo(phot3)
        #plot_photo(xx,yy,title=files[i])
        saved.append( (xx,yy,files[i]) )
    return saved
def sorter(text):
    return
if __name__ == '__main__':
    root1 = '/run/user/1000/gvfs/sftp:host=copernicus.astro.umd.edu,user=antoine/'
        #'Volumes/TESS_05/AW_Data/Sector_60/S60_cam1_ccd4/Didymos/Didymos_coadd/fits/'
    path1 = ''
    path2 = 'Users/antoine/Desktop/stash/rotated/'
    #archived = '/run/user/1000/gvfs/sftp\:host\=copernicus.astro.umd.edu\,user\=antoine/'+\
    #    'Volumes/TESS_05/AW_Data/Sector_60/archive_S60_cam1_ccd4/Didymos/Didymos_coadd/fits/'
    top_data_directory = root1 + path2
    file2 = []
    for a,b,c in os.walk(top_data_directory):
        #print(c)
        file2 = c
    #trying with one
    files = sorted(file2)
    ones = files[:17]
    five = files[17:22]
    thir = files[22:]

    #fill =  fits.open(top_data_directory + 'Didymos_30day_01.fits')
    #dat1 = fill[0].data
    #goo = photo( dat1 )
    #print(goo[0])
    #print(goo[1])
    #print(len(goo[1]))
    #pho1 = goo[1]
    #pho2 = [ i[0] for i in pho1]
    #pho3 = np.array(pho2,dtype=float)
    #xax = np.arange( 200,399, 1,dtype=int)
    #yax = pho3
    #plot_photo(xax,yax)
    #print(yax)
    #do all of them actually
    #for i in range(len(files)):
    #    ifile = fits.open(top_data_directory+files[i])
    #    dati = ifile[0].data
    #    phot1 = photo( dati, radius=2 )
    #    print(phot1)
    #    phot2 = [ i[2] for i in phot1 ]
    #    phot3 = np.array(phot2,dtype=float)
    #    xx,yy = get_photo(phot3)
    #    plot_photo(xx,yy,title=files[i])
    she = photo_all(five,radius=5)
    shd = photo_all(thir,radius=5)
    #print( she[-1][1] )
    for i in range(len(she)): plot_photo( she[i][0],she[i][1],title=five[i] )
    plot_photo( she[0][0],she[0][1],title=thir[0] )
    #plot_photo( she[-3][0],she[-3][1],title=five[-3] )
    #plot_photo( she[7][0],she[7][1],title=files[7] )
    pass
else:
    pass

