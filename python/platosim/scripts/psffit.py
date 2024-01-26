#!/usr/bin/env python
#
# Perform PSF fitting 
# This program was adapted from 'lightcurve' program
# 
# Warning : Input pre-processed imagettes should not be corrected for the sky-background
#
# R. Samadi,  14 March 2016
# $Id$

import os
import re
import sys
import math
import yaml
import h5py
import struct
import getopt

import numpy as np
from numpy.linalg import  LinAlgError
from astropy.io import fits

import pislib3 as pis
import spline2dbase
import psffitlib

import xml.dom.minidom as xdm
from pylab import *




arc2rad = math.pi/3600./180. # arc second to radian

# convert to radians
def deg2rad(d) :
    return d*math.pi/180

# convert to degrees
def rad2deg(r) :
    return r*180/math.pi
    
        

        
def readHdf5ImageFile(idir, prefix):

    """Function to read the HDF5 configuration file.
    """

    try:
        # Read only once the HDF5 file
        h5InputImagePathFilename = f"{idir}/{prefix}-pp.hdf5"
        h5InputImageFile = h5py.File(h5InputImagePathFilename, "r") 
        pis.trace("Read HDF5 file : " + h5InputImagePathFilename)

        # Fetch images
        h5InputArr = h5InputImageFile["Images"]

        # Fetch outlier flags of present
        for key in h5InputImageFile.keys():
            if(key == "Outliers"):
                pis.trace('outlier flags are present')
                OutliersArray = h5InputImageFile["Outliers"]
        if(OutliersArray is None):
            pis.trace('No information about outliers')
            OutliersArray = None
            
        return h5InputArr, OutliersArray
    
    except IOError as e:
            print(e)
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise


    

    
def centroid_offset(psf):

    """Function to find the PSF centroid offset.
 
    Read the point spread function and calculate its centroid, then return
    the offset of the centroid from the centre.
    """
    n,m = psf.shape
    x,y = np.meshgrid(np.arange(0,m)+0.5, np.arange(n)+0.5)

    sx = float(np.sum(psf*x))
    sy = float(np.sum(psf*y))
    sp = float(np.sum(psf))

    return sx/sp - m/2, sy/sp - n/2





def lightcurve(idir, odir, maxmag, Mpx, Nacc, PSFFile, PSFbsres,
               ron, prnuerr, poserr, sposerr, rflxerr, Delta_mag,
               subname, fix_cont_pos=False, fix_cont_amp=False, bgc=True,
               adapt_imagette_size=False, method=0):

    # Fetch paths
    idir = idir.rstrip('/')
    prefix = idir.rpartition('/')[2]

    pis.trace(f"Using prefix {prefix}")
    if(subname !=''):
        subname = '-'+subname
    
    # Read xml file for PSF filename
    dom = xdm.parse("%s/%s.xml" % (idir,prefix))
    
    if(ron is None):
        ron = float(pis.from_dom('ReadoutNoise',dom))
        
    PSFSubPixels = int(pis.from_dom('PSFSubPixels',dom))
    
    if PSFFile == '':
        if(pis.get_first_value('PSFbsres',dom) is not None):
            PSFbsres0 = int(pis.from_dom('PSFbsres',dom))
            if(PSFbsres0>0):
                PSFFile = "%s/%s_psfccd.vec" % (idir,prefix)
                PSFbsres = PSFbsres0
    if bgc:
        pis.trace('Assuming imagettes corrected for the background')
        t_exp    = float(pis.fromdom('IntegrationTime',dom))
        t_trans  = float(pis.fromdom('ChargeTransferTime',dom))
        bgclevel  = float(pis.fromdom('SkyBackground',dom)) * (t_exp+t_trans)
        pis.trace('Background level = %f [e-]' % bgclevel)
    else:
        pis.trace('Assuming imagettes non-corrected for the background')
        bgclevel  = 0.


    # 
    pis.trace(f'Readout noise level : {ron} [e-]')

    # Read Flat field
    path_flat = f'{idir}/{prefix}_ff.fits'
    hdulist = fits.open(path_flat)
    flat = hdulist[0].data
    hdulist.close()
    
    # Flatfield
    if prnuerr >= 0.: # Error on PRNU knowledge
        pis.trace(f'PRNU knowledge error : {prnuerr}')
        flat  +=  np.random.normal(size=flat.shape)*prnuerr/100.
    else: # PRNU not corrected
        flat = np.ones((flat.shape))

    # Read HDF5 file
    ImagesArray, OutliersArray = readHdf5ImageFile(idir, prefix)

    # Number of image
    ImagesKeys = list(ImagesArray.keys())
    N = len(ImagesKeys)
    if N == 0 :
        pis.trace('ERROR : No images found, from time file or hdf5')
        return
    pis.trace(f'each curve with {N:d} time-sample(s)')

    # Use the first fits file to find the size of the images
    img = ImagesArray[ImagesKeys[0]][()].astype(np.float64)
    nx, ny = img.shape[1], img.shape[0]
    pis.trace(f"images    : {nx:d} x {ny:d}")
    pis.trace(f"imagettes : {Mpx:d} x {Mpx:d}")

    # Read the true PSF used by the simulator
    if(PSFFile == ''):
        PSFFile = f'{idir}/{prefix}_psfccd.fits'
    if re.match('.*\.fits$', PSFFile):
        PSFbsres = 20
        # read psf
        pis.trace(f"loading {PSFFile}")
        # PSF CCD using fits format
        f = fits.open(PSFFile)
        psf = np.array(f[0].data,dtype=np.float)
        f.close()

        pis.trace(f'Sub-pixel resolution of the PSF : {PSFSubPixels:d}')
        psf = psf.view(pis.Psf)
        psf.init(subres=PSFSubPixels)
        psf.normalize()   
        pis.trace( psf.info())
        PSFSizex = psf.shape[1]/PSFSubPixels
        PSFSizey = psf.shape[0]/PSFSubPixels
        if( (PSFSizey > Mpx) or (PSFSizex > Mpx) ):
            dx = int( (PSFSizex-Mpx)/2*PSFSubPixels )
            dy = int( (PSFSizey-Mpx)/2*PSFSubPixels )
            psf =   np.ascontiguousarray(psf[dy:-dy,dx:-dx])
            psf = psf.view(pis.Psf)
            psf.init(subres=PSFSubPixels)
            psf.normalize()   
            pis.trace( psf.info())
            PSFSizex = Mpx
            PSFSizey = Mpx
        lx = PSFbsres*PSFSizex
        ly = PSFbsres*PSFSizey
              
        # b-spline representation of the PSF
        psfbs = spline2dbase.Pixel2Spline(psf, lx ,ly,ktype=psffitlib.ktype)

    elif ( re.match('.*\.vec',PSFFile) ):
        pis.trace("loading %s" % PSFFile)
        psfbs = pis.readbinvec(PSFFile)
        s = psfbs.shape[0]
        lx = int(math.sqrt(s))
        ly = lx
        psfbs = psfbs.reshape((ly,lx))
        PSFSizex = int(lx / PSFbsres)
        PSFSizey = int(ly / PSFbsres)
        # the b-splined decomposition is normalized
        img = spline2dbase.Spline2Imagette(psfbs,PSFbsres,PSFSizex,PSFSizey,ktype=psffitlib.ktype)
        bssum = img.sum()
        psfbs /= bssum
    else:
        print("Unknown PSF file type  %s" % PSFFile)
        raise
              
    pis.trace("PSF size X,Y: %f, %f" % (PSFSizex,PSFSizey))

    pis.trace("PSF b-spline resolution: %i" % PSFbsres)        

    # pixel representation of the PSF:
    psf = spline2dbase.Spline2Imagette(psfbs,PSFbsres,PSFSizex,PSFSizey,subres=PSFSubPixels,ktype=psffitlib.ktype)
    PSFcx,PSFcy = pis.barycenter(psf,subres=PSFSubPixels)
    pis.trace("PSF barycenter: %f , %f " % (PSFcx,PSFcy))
    psfsum = float(psf.sum())
    psf /= psfsum
    
    # get centroid offsets in units of subpixel,
    (dx,dy) = centroid_offset(psf)
    pis.trace("PSF centroid ({0:.5f},{1:.5f}) subpixels from PSF centre".format(dx,dy))
    # convert to units of pixel
    dx /= PSFSubPixels
    dy /= PSFSubPixels    # get centroid offsets in units of subpixel,

    PSF = (psfbs,PSFbsres,PSFcx,PSFcy)
    
    # read info file for target stars
    coordpath = "%s/%s_starcoord.dat" % (idir,prefix)
    starsinfo = np.loadtxt(coordpath)
    if(len(starsinfo.shape)==1):
        starsinfo = starsinfo.reshape((1,starsinfo.size))

    Nstars = starsinfo.shape[0]
    pis.trace('number of stars: %i' % Nstars)
    # Error on star centroids:
    if(poserr>0.):
        pis.trace('Random error on star centroid: %f '% (poserr))
        starsinfo[:,3:5] += np.random.normal(size=(Nstars,2))*poserr

    if(sposerr>0.):
        pis.trace('Systematic error on star centroid: %f '% (sposerr))
        starsinfo[:,3] += np.random.normal(size=(1))*sposerr
        starsinfo[:,4] += np.random.normal(size=(1))*sposerr
         
    # Error on star flux:
    if(rflxerr>0.):
        pis.trace('Error on relative star flux: %f '% (rflxerr))
        starsinfo[:,5] +=  -2.5*np.log10((1.+np.random.normal(size=Nstars)*rflxerr/100.))
    
    starsinfo[:,3] += dx
    starsinfo[:,4] += dy
    # id,RA,D,X,Y,M,Xccd,Yccd,Xfp,Yfp
    # process each stars
    for k in range(Nstars):
        if( (maxmag is None) or (starsinfo[k,5]<maxmag) ):
            ID = int(starsinfo[k,0])
            name = "%09i" % (ID)
            pis.trace("*********************************")
            pis.trace("Processing star ID %09i" % ID)
            mag = starsinfo[k,5]
            xc = starsinfo[k,3]
            yc = starsinfo[k,4] 
            # extracting  the information about the contaminants
            distance = np.sqrt( (starsinfo[:,3]-xc)**2 + (starsinfo[:,4]-yc)**2 )
            select = (np.array(starsinfo[:,0],dtype=int) != ID) & (starsinfo[:,5]< mag+Delta_mag) \
                & (distance<8.) 
            Contaminants = []
            Nc = select.sum()
            pis.trace("Number of contaminants %i" % Nc)
            for l in range(Nstars):
                IDc = int(starsinfo[l,0]) 
                if(select[l]):
                    Contaminants.append((IDc,starsinfo[l,3],starsinfo[l,4],starsinfo[l,5]))
            (i0,j0) = (int(round(z - float(Mpx)/2)) for z in (xc,yc))
            if( (i0<0) or (j0<0) or (i0+Mpx>nx) or (j0+Mpx>ny)):
                pis.trace('No enough room in the imagettes, we skip this star')
            else:
                
                Star = (name,ID,xc,yc,mag)
                StarWindow = (i0,j0,Mpx,Mpx)
                #if(True):
                try:
                    results,flagfree = psffitlib.ExtractStarPhotometry(
                        Star,StarWindow,Contaminants,ImagesArray,Nacc,PSF,
                        ron,flat,fix_cont_pos=fix_cont_pos,fix_cont_amp=fix_cont_amp,
                        bgclevel=bgclevel,adapt_imagette_size=adapt_imagette_size,method=method,OutliersArray=OutliersArray)
                    # write lightcurve in ascii dat file
                    lcpath = "%s/%s%s.dat" % (odir,name,subname)
                    header = "flx cx cy  bg flx_err cx_err cy_err bg_err chi2 iter lamb\n"
                    header += ("i0: %i\n" % (i0))
                    header += ("j0: %i\n" % (j0))
                    header += ("xc: %f\n" % (xc-i0))
                    header += ("yc: %f\n" % (yc-j0))
                    header += ("flag: %s\n" % (np.array2string(flagfree)))
                    np.savetxt(lcpath,results,header=header)
                #else:
                except:
                    pis.trace("ERROR: some errors occur, this target is ignored")
                    pis.trace( "%s" % sys.exc_info()[0])
            pis.trace("Done with star ID %09i" % ID)
            



def main() :

    def usage() :
        print("usage : psffit [options] <imagette dir>")
        print()
        print(" -h           : this help")
        print(" -m <mag>     : magnitude threshold")
        print(" -a <integer> : number of imagettes accumulated (default: 1)")
        print(" -n <pixels>  : number of pixels (square) in window")
        print(" -o <dir>     : output directory")
        print(" -f <path>    : use this PSF ")
        print(" -b <integer> : b-spline resolution of the working PSF (default: 20)")
        print(" -p <float>   : knowledge error of the prnu in %. unknown if negative value (default 0.)")
        print(" -r <float>   : readout noise level (in e-, by default taken for the xml file)")
        print(" -c <float>   : random error on the absolute star centroid (in pixel, default: 0)")
        print(" -s <float>   : systematic error on the absolute star centroid (in pixel, default: 0)")
        print(" -F <float>   : relative error (in %) on the star flux (default: 0)")
        print(" -D <float>   : Delta_magnitude (default: 3), contaminant stars fainter than P(target) + Delta_magnitude  are ignored. P(target) is the magnitude of the targets ")
        print(" -S <string>  : string added to the output file, '-S pfit' will name output file as '<starid>-pfit.dat' (default: '')")
        print(" -K <int>     : type of Knots (0-> 'simple', 1-> Dierckx's distribution) (default: 0)")
        print(" -M <int>     : minimization method : 0->Levenberg-Marquardt (default), 1-> BFGS")
        print(" --no-bg-corr : assumes that the background has not be corrected on the input images")
        print(" --adapt-imagette-size : adapt the imagette size to the star magnitude: 5x5 for 12<=mag<13, 4x4 for 13<=mag<15, and 3x3 for mag>=15")
        print(" --fix-relative-pos    : maintain fixed the relative positions of the contaminants (except for the first exposure")
        print(" --fix-intensity       : maintain fixed the intensity of each contaminant (except for the first exposure")
        print(" --seed <value>  : impose this seed value")

    try:
        opts,args = getopt.getopt(sys.argv[1:],"ha:m:n:o:B:b:f:p:r:c:F:s:D:S:K:M:",
                                  ["fix-cont-pos","fix-cont-amp","no-bg-corr",
                                   "adapt-imagette-size","fix-relative-pos","fix-intensity","seed="])

    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    maxmag  = None
    idir    = "."
    odir    = "."
    Mpx     = 6
    PSFbsres = 20
    PSFFile = ''
    prnuerr = 0.
    poserr = 0.
    sposerr = 0.
    rflxerr = 0.
    Delta_mag = 3.
    ron = None
    fix_cont_pos = False
    fix_cont_amp = False
    subname = ''
    nacc = 1 
    psffitlib.ktype = 0
    bgc = True # assumes the background is corrected 
    adapt_imagette_size = False
    method = 0
    seed = None
    for o, a in opts:
        if o == "-h" :
            usage()
            sys.exit()
        elif o == "-S" :
            subname = a
        elif o == "-K" :
            psffitlib.ktype = int(a)
        elif o == "-a" :
            nacc = int(a)
        elif o == "-m" :
            maxmag = float(a)
        elif o == "-M" :
            method = int(a)
        elif o == "-D" :
            Delta_mag = float(a)
        elif o == "-p" :
            prnuerr = float(a)
        elif o == "-c" :
            poserr = float(a)
        elif o == "-s":
            sposerr = float(a)
        elif o == "-F" :
            rflxerr = float(a)
        elif o == "-r" :
            ron = float(a)
        elif o == "-n" :
            Mpx = int(a)
        elif o == "-b" :
            PSFbsres = int(a)
        elif o == "--fix-cont-pos" :
            fix_cont_pos = True  
        elif o == "--fix-cont-amp" :
            fix_cont_amp = True  
        elif o == "-o" :
            odir = a
        elif o == "-f" :
            PSFFile = a
        elif o == "--no-bg-corr" :
            bgc = False
        elif o == "--adapt-imagette-size":
            adapt_imagette_size = True
        elif o == "--fix-relative-pos":
            fix_cont_pos = True
        elif o == "--fix-intensity":
            fix_cont_amp = True
        elif o == "--seed" :
            seed = int(a)
        else:
            print("unhandled option %s" % (o))
            sys.exit(1) 

    nargs = len(args)

    if nargs > 1 :
        print("too many arguments")
        usage()
        sys.exit()

    if nargs == 1 :
        idir = args[0]

    if nargs == 0 :
        print("imagette directory must be specified")
        usage()
        sys.exit()

    pis.trace("This is PSFfit")
    pis.trace("input directory : %s" % (idir))
    pis.trace("output directory : %s" % (odir))

    if maxmag :
        pis.trace("magnitude threshold {0:.3f}".format(maxmag))
    else :
        pis.trace("no magnitude threshold")

    np.random.seed(seed)

    lightcurve(idir,odir,maxmag,Mpx,nacc,PSFFile,PSFbsres,ron,prnuerr,poserr,sposerr,
               rflxerr,Delta_mag,subname,fix_cont_pos=fix_cont_pos,
               fix_cont_amp=fix_cont_amp,bgc=bgc,adapt_imagette_size=adapt_imagette_size,method=method)

    pis.trace("Done.")

if __name__ == "__main__":
    main()

