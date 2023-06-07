#going to read in two fits files, rotate one of them, divide to compare or something. pain
import numpy as np
from astropy.io import fits
from scipy import ndimage
from scipy import signal
#from PIL import Image
#import matplotlib
import matplotlib.pyplot as plt

def readFits(filename,ext=0):
    inn = fits.open(filename)
    hdr = inn[ext].header
    dat = inn[ext].data
    return (dat,hdr)
#load in data
dir_new = '/Volumes/TESS_05/AW_DATA/Sector_60/S60_cam1_ccd4/Didymos/Didymos_coadd/fits/'
dir_old = '/Volumes/TESS_05/AW_DATA/Sector_60/archive_S60_cam1_ccd4/Didymos/Didymos_coadd/fits/'
file1 = 'Didymos_30day_01.fits'
file2 = ''
im_old, h_old = readFits(dir_old+file1)
im_new, h_new = readFits(dir_new+file1)

#rotating older one
im_rot = ndimage.rotate(im_old,189,reshape=True)

#comparing
#diff = im_rot - im_new
diff= 1

#convolvo
kern = np.ones((3,3),dtype=float) / 9.
arr2 = ndimage.convolve( im_new, kern, mode='nearest' )
print(arr2)

#plotting
fig,(ax1,ax2,ax3) = plt.subplots(1,3)

ax1.imshow(im_rot)
#plt.show(block=False)
ax2.imshow(im_new)
#ax3.imshow(diff)
plt.show(block=True)


