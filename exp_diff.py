#take image, rotate it, compare it to old version
import os
import numpy as np
from astropy.io import fits
from scipy import ndimage
import matplotlib.pyplot as plt

def usweall():
    return
def rotcomp(im1, im2, angle):
    #will rotate im1 according to angle, then subtract it from im2
    im_1 = ndimage.rotate(im1,angle,reshape=False)
    im_diff = im2 - im_1
    return im_diff
if __name__ == '__main__':
    solong = '/run/user/1000/gvfs/sftp:host=copernicus.astro.umd.edu,user=antoine/'+\
        'Volumes/TESS_05/AW_Data/Sector_60/S60_cam1_ccd4/Didymos/Didymos_coadd/fits/'
    archive = '/run/user/1000/gvfs/sftp:host=copernicus.astro.umd.edu,user=antoine/'+\
        'Volumes/TESS_05/AW_Data/Sector_60/archive_S60_cam1_ccd4/Didymos/Didymos_coadd/fits/'
    saving_here = '/run/user/1000/gvfs/sftp:host=copernicus.astro.umd.edu,user=antoine/'+\
        'Volumes/TESS_05/AW_Data/didy_diffs/'
    for i,j,k in os.walk(solong): files = k
    ims1_filename = [ archive+i for i in files]
    ims2_filename = [ solong+i for i in files]
    #
    for i in range(len(ims2_filename)):
        im1f = fits.open( ims1_filename[i] )
        im2f = fits.open( ims2_filename[i] )
        im1  = im1f[0].data
        im2  = im2f[0].data
        differ = rotcomp(im1,im2,189)
        saving = fits.PrimaryHDU(differ)
        saving.writeto(saving_here+files[i])
    pass
else:
    pass
