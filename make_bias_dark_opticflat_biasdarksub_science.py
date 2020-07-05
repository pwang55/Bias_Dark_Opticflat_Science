"""

Usage:

    In python script directory:
    $ python make_bias_dark_opticflat_biasdarksub_science.py path_to_file/filelist.txt type

    In target directory:
    $ python path_to_script/make_bias_dark_opticflat_biasdarksub_science.py filelist.txt type

    type: bias, dark, opticflat, biasdarksub, science, biasdarksub2science

This python script defines functions to make bias.fits, dark300.fits, dark600.fits, optic flats, 
biasdarksub_xxx.fits and science_xxx.fits. For unnecessary re-running this script, 
it will find and use existing bias/dark/optic/twilight flat... from environment variable $_90PRIME_BIAS_DARK_FLATS_DIR

Note: "science" here means corrected frame, not the final masked science frame for stacking


"""

import sys
import os
from numpy import *
from astropy.io import fits
from scipy.signal import savgol_filter
import re


# If input arguments arenn't exactly 2 of them (which makes len(sys.argv)=3), print docstring and end script, else take input arguments
if len(sys.argv) !=3:
    print(__doc__)
    exit()
else:
    inputfilename = sys.argv[1]     # filelist containing all the fits files to be combined, for example biaslist.txt
    reduction_type = sys.argv[2]    # bias, dark, opticflat, biasdarksub, science




if not((reduction_type =='bias') or (reduction_type == 'dark') or (reduction_type == 'opticflat') or (reduction_type == 'biasdarksub') or (reduction_type == 'science') or (reduction_type =='biasdarksub2science')):
    print('\n\tWrong Reduction Type. Choose from "bias", "dark", "opticflat", "biasdarksub", "science", "biasdarksub2science"\n')
    exit()



RN = 10.0


# Get absolute path of this python script
script_path = sys.path[0]
# Get directory path to bias_dark_flats
bias_dark_flats_dir = os.environ['_90PRIME_BIAS_DARK_FLATS_DIR']
# Get target directory path that end with /, if inputfilename is not full path but just biaslist.txt, target path set to empty
if len(inputfilename.split('/')) > 1:
    target_path = inputfilename[0:-len(inputfilename.split('/')[-1])]
elif len(inputfilename.split('/')) == 1:
    target_path = ""


# Reading input files into an empty list and put target_path in front of each of them
with open(inputfilename) as f:
    filelists = [ target_path+l.strip() for l in f ]



# Define function to find overscan region; these n number will be the exact pixel number (start 0) of where biassec start and finish, so if using these in argument, the end has to add 1
# 90prime data NAXIS1 = y, NAXIS2 = x, so in header [NAXIS1, NAXIS2], but when using loaded array, dat[naxis2, naxis1]
# default to use the end of naxis1 = y as overscan
def find_overscan(biassec, datasec):
    biasstr = re.split('[\[:,\]]', biassec)
    datastr = re.split('[\[:,\]]', datasec)
    nyb0 = int(biasstr[1])-1
    nyb1 = int(biasstr[2])-1
    nxb0 = int(biasstr[3])-1
    nxb1 = int(biasstr[4])-1
    nyd0 = int(datastr[1])-1
    nyd1 = int(datastr[2])-1
    nxd0 = int(datastr[3])-1
    nxd1 = int(datastr[4])-1
    return nxb0, nxb1, nyb0, nyb1, nxd0, nxd1, nyd0, nyd1




# Define function to make bias.fits
def make_bias():
    n = len(filelists)          # number of bias files
    with fits.open(filelists[0]) as f:
        namp = len(f)           # number of amplifiers+1 (including Primary HDU)
        ovsy = int(f[1].header['ovrscan1'])
        ovsx = int(f[1].header['ovrscan2'])
        ny = f[1].header['naxis1'] - ovsy
        nx = f[1].header['naxis2'] - ovsx

    dat = zeros((n, nx, ny), dtype='float32')
    hlist = []
    for j in range(1,namp):         # iterate each amplifier
        for i in range(n):          # iterate each bias files
            print("reading extension {} in {}".format(j, filelists[i]))
            with fits.open(filelists[i]) as b:
                nxb0, nxb1, nyb0, nyb1, nxd0, nxd1, nyd0, nyd1 = find_overscan(b[j].header['biassec'], b[j].header['datasec'])
                dat[i] = b[j].data[nxd0:nxd1+1, nyd0:nyd1+1]
                overscan0 = average(b[j].data[nxd0:nxd1+1, nyb0:nyb1+1], axis=1)
                overscan1 = savgol_filter(overscan0, 201, 2)
                dat[i] = dat[i] - overscan1[:,None]

        hduI = fits.ImageHDU()
        hduI.data = median(dat, axis=0)
        hduI.header = fits.getheader(filelists[0], j)
        hduI.header['BZERO'] = 0.
        hlist.append(hduI)

    hdu0 = fits.PrimaryHDU()
    hdu0.header = fits.getheader(filelists[0])
    hlist.insert(0, hdu0)
    hduA = fits.HDUList(hlist)
    hduA.writeto(target_path+'bias.fits')




# Define function to make dark300.fits and dark600.fits
def make_dark(exptime=300):
    darkfiles = [ files for files in filelists if (fits.getheader(files)['EXPTIME']==exptime) ]

    print('Making dark'+str(exptime)+'.fits with following files:\n')
    for fname in darkfiles:
        print(fname)
    print('\n')

    n = len(darkfiles)
    with fits.open(darkfiles[0]) as f:
        namp = len(f)           # number of amplifiers+1 (including Primary HDU)
        ovsy = int(f[1].header['ovrscan1'])
        ovsx = int(f[1].header['ovrscan2'])
        ny = f[1].header['naxis1'] - ovsy
        nx = f[1].header['naxis2'] - ovsx

    dat = zeros((n, nx, ny), dtype='float32')
    hlist = []
    for j in range(1,namp):     # iterate each amplifier
        for i in range(n):      # iterate each bias files
            print ("reading extension {} in {}".format(j, darkfiles[i]))
            with fits.open(darkfiles[i]) as d:
                nxb0, nxb1, nyb0, nyb1, nxd0, nxd1, nyd0, nyd1 = find_overscan(d[j].header['biassec'], d[j].header['datasec'])
                dat[i] = d[j].data[nxd0:nxd1+1, nyd0:nyd1+1]
                overscan0 = average(d[j].data[nxd0:nxd1+1, nyb0:nyb1+1], axis=1)
                overscan1 = savgol_filter(overscan0, 201, 2)
                dat[i] = dat[i] - overscan1[:,None] - fits.getdata(bias_dark_flats_dir+'/bias/bias.fits', j)  # each frame subtract bias.fits

        hduI = fits.ImageHDU()
        hduI.data = median(dat, axis=0)
        hduI.header = fits.getheader(darkfiles[0], j)
        hduI.header['BZERO'] = 0.
        hlist.append(hduI)

    hdu0 = fits.PrimaryHDU()
    hdu0.header = fits.getheader(darkfiles[0])
    hlist.insert(0, hdu0)
    hduA = fits.HDUList(hlist)
    hduA.writeto(target_path+'dark'+str(exptime)+'.fits')




# Define function to make optic flats
def make_flat(filter='ASU1'):
    flatfiles = [ files for files in filelists if (fits.getheader(files)['FILTER']==filter) ]
    print('Making optic_flat_'+filter.lower()+'.fits with following files:\n')
    for fname in flatfiles:
        print(fname)
    print('\n')

    n = len(flatfiles)
    with fits.open(flatfiles[0]) as f:
        namp = len(f)           # number of amplifiers+1 (including Primary HDU)
        ovsy = int(f[1].header['ovrscan1'])
        ovsx = int(f[1].header['ovrscan2'])
        ny = f[1].header['naxis1'] - ovsy
        nx = f[1].header['naxis2'] - ovsx

    dat = zeros((n, nx, ny), dtype='float32')
    hlist = []

    for j in range(1, namp):        # iterate each amplifier
        for i in range(n):          # iterate each bias files
            print ("reading extension {} in {}".format(j, flatfiles[i]))
            with fits.open(flatfiles[i]) as d:
                nxb0, nxb1, nyb0, nyb1, nxd0, nxd1, nyd0, nyd1 = find_overscan(d[j].header['biassec'], d[j].header['datasec'])
                im = d[j].data[nxd0:nxd1+1, nyd0:nyd1+1]
                overscan0 = average(d[j].data[nxd0:nxd1+1, nyb0:nyb1+1], axis=1)
                overscan1 = savgol_filter(overscan0, 201, 2)
                im = im - overscan1[:,None] - fits.getdata(bias_dark_flats_dir+'/bias/bias.fits', j)    # each frame subtract bias.fits, doesn't need to use dark since flats are short exposure
                dat[i] = im/percentile(im,99)

        hduI = fits.ImageHDU()
        datm = median(dat,axis=0)
        datm = datm/percentile(datm,99)
        hduI.data = datm
        hduI.header = fits.getheader(flatfiles[0], j)
        hduI.header['BZERO'] = 0.
        hlist.append(hduI)

    hdu0 = fits.PrimaryHDU()
    hdu0.header = fits.getheader(flatfiles[0])
    hlist.insert(0, hdu0)
    hduA = fits.HDUList(hlist)
    hduA.writeto(target_path+'optic_flat_'+filter.lower()+'.fits')





# Define function to make biasdarksub frame from raw science frame, only subtract bias and sub, file wil start with biasdarksub_
def make_biasdarksub():

    with fits.open(filelists[0]) as f:
        namp = len(f)                    # number of amplifiers+1 (including Primary HDU)
        ovsy = int(f[1].header['ovrscan1'])
        ovsx = int(f[1].header['ovrscan2'])
        ny = f[1].header['naxis1'] - ovsy
        nx = f[1].header['naxis2'] - ovsx

    n = len(filelists)

    for i in range(len(filelists)):
        print('Subtract bias and dark for '+filelists[i].split('/')[-1]+' '+str(i)+'/'+str(n))


        exptime = fits.getheader(filelists[i])['exptime']
        biasfile = bias_dark_flats_dir+'/bias/bias.fits'
        if exptime == 300.0:
                darkfile = bias_dark_flats_dir+'/dark/dark300.fits'
        elif exptime == 600.0:
                darkfile = bias_dark_flats_dir+'/dark/dark600.fits'
        else:
                print('no proper dark frame for this file!')
                exit()


        with fits.open(filelists[i]) as im0:
            hlisti = [] 
            for j in range(1, namp):
                nxb0, nxb1, nyb0, nyb1, nxd0, nxd1, nyd0, nyd1 = find_overscan(im0[j].header['biassec'], im0[j].header['datasec'])
                im = im0[j].data[nxd0:nxd1+1, nyd0:nyd1+1]

                overscan0 = average(im0[j].data[nxd0:nxd1+1, nyb0:nyb1+1], axis=1)
                overscan1 = savgol_filter(overscan0, 201, 2)

                b = fits.getdata(biasfile, j)
                d = fits.getdata(darkfile, j)

                dat = (im - overscan1[:,None]- b - d)

                gain = float(im0[0].header['GAIN'+str(j)])

                hduIi = fits.ImageHDU()
                hduIi.data = dat
                hduIi.header = im0[j].header
                hduIi.header['BZERO'] = 0.0
                hlisti.append(hduIi)

            hdu0i = fits.PrimaryHDU()
            hdu0i.header = im0[0].header
            hlisti.insert(0,hdu0i)

        hduAi = fits.HDUList(hlisti)
        hduAi.writeto(target_path+'biasdarksub_'+filelists[i].split('/')[-1])




# Define function to make complete science frame that start with corrected_
def make_science():

    with fits.open(filelists[0]) as f:
        namp = len(f)                    # number of amplifiers+1 (including Primary HDU)
        ovsy = int(f[1].header['ovrscan1'])
        ovsx = int(f[1].header['ovrscan2'])
        ny = f[1].header['naxis1'] - ovsy
        nx = f[1].header['naxis2'] - ovsx

    n = len(filelists)

    for i in range(len(filelists)):
        print('Making corrected frame for '+filelists[i].split('/')[-1]+' '+str(i)+'/'+str(n))

        filt = fits.getheader(filelists[i])['filter'].lower()
        exptime = fits.getheader(filelists[i])['exptime']
        month = fits.getheader(filelists[i])['date'].split('-')[1]
        if month == '02':
            month_dir = '/twilight_flats/Feb/'
            month_string = 'feb'
        elif month == '03':
            month_dir = '/twilight_flats/Mar/'
            month_string = 'mar'
        flatfile = bias_dark_flats_dir+month_dir+'flat_'+filt+'_'+month_string+'.fits'

        biasfile = bias_dark_flats_dir+'/bias/bias.fits'
        if exptime == 300.0:
                darkfile = bias_dark_flats_dir+'/dark/dark300.fits'
        elif exptime == 600.0:
                darkfile = bias_dark_flats_dir+'/dark/dark600.fits'
        else:
                print('no proper dark frame for this file!')
                exit()


        with fits.open(filelists[i]) as im0:
            hlisti = []
            for j in range(1, namp):
                nxb0, nxb1, nyb0, nyb1, nxd0, nxd1, nyd0, nyd1 = find_overscan(im0[j].header['biassec'], im0[j].header['datasec'])
                im = im0[j].data[nxd0:nxd1+1, nyd0:nyd1+1]

                overscan0 = average(im0[j].data[nxd0:nxd1+1, nyb0:nyb1+1], axis=1)
                overscan1 = savgol_filter(overscan0, 201, 2)

                b = fits.getdata(biasfile, j)
                d = fits.getdata(darkfile, j)
                f = fits.getdata(flatfile, j)

                dat = (im - overscan1[:,None]- b - d)/f

                gain = float(im0[0].header['GAIN'+str(j)])

                hduIi = fits.ImageHDU()
                hduIi.data = dat
                hduIi.header = im0[j].header
                hduIi.header['BZERO'] = 0.0
                hlisti.append(hduIi)

            hdu0i = fits.PrimaryHDU()
            hdu0i.header = im0[0].header
            hlisti.insert(0,hdu0i)

        hduAi = fits.HDUList(hlisti)
        hduAi.writeto(target_path+'corrected_'+filelists[i].split('/')[-1])



# Define function that takes biasdarksub and make it corrected frame (only perform flat division)
def biasdarksub2science():
    for i in range(len(filelists)):
        print('Making corrected frame from: \t'+filelists[i].split('/')[-1])

        filt = fits.getheader(filelists[i])['filter'].lower()
        month = fits.getheader(filelists[i])['date'].split('-')[1]
        if month == '02':
            month_dir = '/twilight_flats/Feb/'
            month_string = 'feb'
        elif month == '03':
            month_dir = '/twilight_flats/Mar/'
            month_string = 'mar'
        flatfile = bias_dark_flats_dir+month_dir+'flat_'+filt+'_'+month_string+'.fits'

        with fits.open(filelists[i]) as im0:
            hlisti = []
            for j in range(1,17):
                im = im0[j].data
                flat = fits.getdata(flatfile, j)
                dat = im/flat

                hduIi = fits.ImageHDU()
                hduIi.data = dat
                hduIi.header = im0[j].header
                hduIi.header['BZERO'] = 0.0
                hlisti.append(hduIi)

            hdu0i = fits.PrimaryHDU()
            hdu0i.header = im0[0].header
            hlisti.insert(0,hdu0i)

        hduAi = fits.HDUList(hlisti)
        hduAi.writeto(filelists[i].replace('biasdarksub_','corrected_'))



if reduction_type == 'bias':
    make_bias()
if reduction_type == 'dark':
    make_dark(300)
    make_dark(600)
if reduction_type == 'opticflat':
    make_flat('ASU1')
    make_flat('ASU2')
    make_flat('ASU3')
    make_flat('ASU4')
if reduction_type == 'biasdarksub':
    make_biasdarksub()
if reduction_type == 'science':
    make_science()
if reduction_type == 'biasdarksub2science':
    biasdarksub2science()
 


