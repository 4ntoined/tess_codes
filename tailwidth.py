#going to measure the tail brightness
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelmax
from astropy.io import fits

def gauss(x, sigma, mean):
    ans =  np.exp( -(x-mean)**2/(2.*sigma**2.)) / (sigma * np.sqrt(2*np.pi))
    return ans
def liner(x, m, b):
    ans = m*x+b
    return ans
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
def aperture_scan(frame, scan_location, aperture_size=1):
    ##1d data
    return
def wider(frame,radius=0):
    """
    will take in [frame], measure the brightness in central 9 pixels
    and then track and measure the tail brightness
    frame : 2d array where comet is centered
    """
    window_1 = 20           #radius of intial aperture
    continuum_long = 11     #number of pixels per flank to use for cont'm
    window_2 = 15           #central range to search for peaks
    sho = ( window_1*2+1 - window_2 ) // 2
    ny, nx = frame.shape
    centerx, centery = (nx-1)//2, (ny-1)//2
    #center
    #aper = square_ap( (centerx,centery), radius=radius )
    #center9 = frame.copy()[ aper[1][0]:aper[1][1],aper[0][0]:aper[0][1] ]
    #boy = photo2(center9)
    #nan_catch = np.count_nonzero( np.isnan( center9 ))
    #sum_center = np.nansum( center_9 )
    #center_phot = (sum_center, nan_catch)
    #for these 401x401 images, center is index 200, for radius 1, we can go as far
    #as center,index=399, or size-2, or size-(radius+1), for 9 and radius 2, middle is 4, end is 6!
    #all
    liss = []
    peaks = []
    widths = []
    bounder = (centerx, nx) #xaxis
    for i in range( bounder[0], bounder[1] ):
        #print(i)
        #thisone = ( centery, i  )
        bounds = (centery-window_1,centery+window_1+1) #yaxis
        span = bounds[1]-bounds[0]
        if radius==0:
            vertline = frame.copy()[bounds[0]:bounds[1],i]
        else:
        #### turn frame into vertline using aperture
            vertt =  [ np.nanmean( frame.copy()[ii-radius:ii+radius+1,i-radius:i+radius+2] ) for ii in range(bounds[0]+radius,bounds[1]-radius) ]
            vertline= np.array(vertt,dtype=float)
        vertline_i = np.arange((1-span)//2, vertline.size-span//2, 1)
        #if i == 208:
        #plt.plot(vertline_i,vertline)
        #plt.show()
        contx =  np.concatenate( [vertline_i[:continuum_long], vertline_i[-continuum_long:]] )
        contin = np.concatenate( [vertline[:continuum_long], vertline[-continuum_long:]] )
        try:
            #fitting = curve_fit(gauss, vertline_i, vertline, p0=(1.,vertline.size//2))
            #find continuum
            fit_par, fit_err = curve_fit(liner, contx, contin, p0=(0.,0.))
            confit = liner(vertline_i, fit_par[0], fit_par[1])
            #remove continuum
            subtracted = vertline - confit
            #detect location of max
            #print(subtracted)
            peakarea = subtracted[sho:-sho].copy()
            peaki = np.arange(peakarea.size,dtype=int)
            #
            #print(peakarea)
            extras = argrelmax(peakarea, order=2)
            #print(extras)
            exi = extras[0]  #i's ref peakarea
            exv = [peakarea[i] for i in exi] #i's reference exi
            exv = np.array(exv,dtype=float)
            #print(exv)
            maxpeak = np.squeeze( np.argwhere( exv == np.max(exv) ))
            peaks.append( exv[maxpeak] )
            #print(exv)
            #maxpeak = np.max( exv )
            #print('maxes at:', extras[0]-7)
            #print('max max loc:',exi[maxpeak]-7,'| maxmax value:', exv[ maxpeak ] )
            #take half the max, plot hori line there, see where peak intersects
            halfmax = exv[maxpeak] / 2.
            #find lesser intersect
            verge = True
            check = exi[maxpeak]
            while 1:
                #see if lightcurve dips below halfmax
                if peakarea[check] < halfmax:
                    #verge=True
                    break
                else:
                    check-=1
                if check <= -1:
                    verge=False
                    break
            lowend = check
            #find greater intersect
            check = exi[maxpeak]
            while 1: #check < peaki[-1]+1:
                if peakarea[check] < halfmax:
                    break
                else:
                    check+=1
                if check >= peaki[-1]+1:
                    verge=False
                    break
            highend = check
            if verge==False:
                #looking for intersects failed!
                widths.append(highend-lowend)
                liss.append((exv[maxpeak],highend-lowend))
            else:
                #print(lowend-7,highend-7)
                ### lowend 2 ###
                y1, y2 = peakarea[lowend:lowend+2]
                x1, x2 = lowend, lowend+1
                x0_low = (halfmax-y1) * ( (x2-x1) / (y2-y1) ) + x1
                #print(y1, y2, x1, x2, x0_low)
                ### highend 2 ###
                y1, y2 = peakarea[highend-1:highend+1]
                x1, x2 = highend-1, highend
                x0_high = (halfmax-y1) * ( (x2-x1) / (y2-y1) ) + x1
                #print(y1, y2, x1, x2, x0_high)
                wid = x0_high-x0_low
                #print(wid)
                widths.append(x0_high-x0_low)
                liss.append((exv[maxpeak],wid))
                pass
        except ValueError:
            #print(f"Fitting failed: column {i}")
            #opt = (-99,-99)
            #cov = (-9.9,-9.9)
            peaks.append(-99)
            widths.append(4)
            pass
        else:
            #plt.plot(vertline_i,vertline,color='b')
            #plt.plot(vertline_i,confit,color='r')
            #plt.plot(vertline_i,subtracted,color='g')
            #plt.plot(contx, contin, color='k')
            #plt.show()
            pass
            #opt, cov = fitting
        pass
        #widths.append( (opt,cov) )
        #aper = square_ap( (i, centery), radius=radius )
        #data9 = frame.copy()[ aper[1][0]:aper[1][1],aper[0][0]:aper[0][1] ]
        #photom = photo2(data9)
        #liss.append(photom)
    return (peaks, widths)
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
def plot_photo(datax, datay, title='', savehere=''):
    """
    datax, datay: 1d arrays of floats to have plotted
    """
    fig, ax = plt.subplots()
    #fig.figsize=
    fig.dpi=140
    #
    ax.plot(datax, datay)
    #ax.set_yscale('log')
    #ax.set_xlim(0,50)
    #
    if title: ax.set_title(title)
    ax.set_xlabel('pixels from center')
    ax.set_ylabel('width of asteroid and tail')
    if savehere:
        plt.savefig(savehere+title+'.png')
    plt.show()
    return
def photo_all(files):
    saved = []
    c1=[]
    c2=[]
    cc=[]
    for i in range(len(files)):
        ifile = fits.open(top_data_directory+files[i])
        dati = ifile[0].data
        p1, p2 = wider( dati, radius=3 )
        #print(phot1)
        #photos = np.array(phot1, dtype=float)
        peaks = np.array(p1,dtype=float)
        bads = [ type(i) == int for i in p2 ]
        wides = np.where(bads,np.nan,p2)
        c1.append(peaks)
        c2.append(wides)
        cc.append((peaks,wides))
        #plt.plot( np.log(peaks) )
        #plt.yscale='log'
        plt.show()
        #print(phot1)
        #phot2 = [ i[0][0] for i in phot1 ]
        #phot_err = [ np.sqrt(np.diag(i[1]))[0] for i in phot1]
        #phot3 = np.array(phot2,dtype=float)
        #xx,yy = get_photo(phot3)
        #plot_photo(xx,yy,title=files[i])
        #saved.append( (xx,yy,files[i],phot_err) )
        pass
    #pack = (c1,c2)
    return cc
if __name__ == '__main__':
    solong = '/run/user/1000/gvfs/sftp:host=copernicus.astro.umd.edu,user=antoine/'+\
        'Volumes/TESS_05/AW_Data/Sector_60/S60_cam1_ccd4/Didymos/Didymos_coadd/fits/'
    alc = '/run/user/5317/gvfs/sftp:host=copernicus.astro.umd.edu/'+\
        'Volumes/TESS_05/AW_Data/Sector_60/S60_cam1_ccd4/Didymos/Didymos_coadd/fits/' 
    savinghere = '/run/user/1000/gvfs/sftp:host=copernicus.astro.umd.edu,user=antoine/'+\
        'Volumes/TESS_05/AW_Data/didy_widths/'
    #archived = '/run/user/1000/gvfs/sftp\:host\=copernicus.astro.umd.edu\,user\=antoine/'+\
    #    'Volumes/TESS_05/AW_Data/Sector_60/archive_S60_cam1_ccd4/Didymos/Didymos_coadd/fits/'
    top_data_directory = solong
    #top_data_directory = alc
    files = []
    for a,b,c in os.walk(top_data_directory):
        #print(c)
        files = c
    #trying with one
    print(files)
    ones = files[:17]
    five = files[17:22]
    thir = files[22:]
    #print(ones,'\n',five,'\n',thir)
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
    she = photo_all(five)
    #print(she[0])
    #print(she[1])
    fig,ax=plt.subplots()
    ax.plot( she[1][1] )
    ax.set_xlim((-1,50))
    plt.show()
    #she0 = she[0]
    #print( she[-1][1] )
    #plot_photo( she0[0],she0[1],title=files[0] )
    #plot_photo( she[-3][0],she[-3][1],title=files[-3] )
    #for i in range(len(she)):
    #    plot_photo( she[i][0],she[i][1],title=files[i] )
    #plot_photo( she[7][0],she[7][1],title=files[7] )
    #xx = np.arange(-10,10)
    #yy = gauss(xx,1.2,0.5)
    #plt.plot(xx,yy)
    #plt.show()
    #we = photo()

    pass
else:
    pass

