#this will take my fits images, excise a certain center, save it as png
import os
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits



def readin(filename,xsize=121,ysize=61,savename=''):
    ifile = fits.open(filename)
    header=ifile[0].header
    img = ifile[0].data
    plt.imshow(img,origin='lower',vmin=0.,vmax=np.max(img)*0.05,norm='linear')
    plt.show()
    ##
    ny,nx = img.shape
    banky = (ny-ysize) // 2
    bankx = (nx-1) // 2 - 20
    ##
    center = img[banky:-banky,bankx:bankx+xsize]
    plt.imshow(center,origin='lower',vmin=0.,vmax=np.max(center)*0.05,norm='linear')
    plt.axis('off')
    if savename: plt.savefig(savename,dpi=100,bbox_inches='tight',pad_inches=0.0)
    plt.show()
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
    readin(top_data_directory+five[1],savename='s05_2.png')
    readin(top_data_directory+five[2],savename='s05_3.png')
    readin(top_data_directory+five[4],savename='s05_5.png')
    readin(top_data_directory+thir[0],savename='s05_1.png')
else:
    pass

