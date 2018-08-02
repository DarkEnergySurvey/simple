#!/usr/bin/env python
"""
Simple binning search algorithm utilities
"""
__author__ = "Keith Bechtol, Sid Mau"

# Python libraries
import os
import glob
import yaml
import numpy as np
#from matplotlib import mlab
import numpy.lib.recfunctions
import healpy as hp
import astropy.io.fits as pyfits # migrate to fitsio
import fitsio as fits
import sys
import pylab as plt

import numpy as np
from operator import add
from scipy import interpolate
from scipy.signal import argrelextrema
import scipy.ndimage

import pylab as plt
import pyfits
import matplotlib
from matplotlib import mlab
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec

# Ugali libraries
#import ugali.utils.mlab
import ugali.utils.healpix

import simple.filters
import simple.simple_utils
import simple.diagnostic_plots

########################################################################

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']
    isoname = cfg[survey]['isoname']
    isosurvey = cfg[survey]['isosurvey']
    mag_max = cfg[survey]['mag_max']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']

    mode = cfg[survey]['mode']
    sim_population = cfg[survey]['sim_population']
    sim_dir = cfg[survey]['sim_dir']

    fracdet_map = cfg[survey]['fracdet']
    
    mag_g = cfg[survey]['mag_g']
    mag_r = cfg[survey]['mag_r']
    mag_g_err = cfg[survey]['mag_g_err']
    mag_r_err = cfg[survey]['mag_r_err']
    
########################################################################

try:
    ra_select, dec_select = float(sys.argv[1]), float(sys.argv[2])
except:
    sys.exit('ERROR! Coordinates not given in correct format.')

# Now cut for a single pixel
pix_nside_select = ugali.utils.healpix.angToPix(nside, ra_select, dec_select)
pix_nside_neighbors = np.concatenate([[pix_nside_select], hp.get_all_neighbours(nside, pix_nside_select)])

# Construct data
file_array = []
for pix_nside in pix_nside_neighbors:
    inlist = glob.glob('{}/*_{:05d}.fits'.format(datadir, pix_nside))
    for infile in inlist:
        if not os.path.exists(infile):
            continue
        file_array.append(infile)

data_array = []
for infile in file_array:
    data_array.append(fits.read(infile))

data = np.concatenate(data_array)
#data = filters.dered_mag(survey, data)

#filter = filters.star_filter(survey, data)
#data = data[filter]

###
#data = data[(data['MAG_PSF_R'] < 90) & (data['MAG_PSF_G'] < 90) & (data['MAG_PSF_I'] < 90)]
data = data[(data['MAG_PSF_G'] < 90) & (data['MAG_PSF_I'] < 90)]

#print("EXPNUM_G: {}".format(np.unique(data['EXPNUM_G'])))
#print("EXPNUM_R: {}".format(np.unique(data['EXPNUM_R'])))
#print("CCDNUM_G: {}".format(np.unique(data['CCDNUM_G'])))
#print("CCDNUM_R: {}".format(np.unique(data['CCDNUM_R'])))

## Scatter Plot
##proj = ugali.utils.projector.Projector(ra_select, dec_select)
##for expnum in np.unique(data['EXPNUM_R']):
##    x, y = proj.sphereToImage(data[basis_1][data['EXPNUM_R'] == expnum], data[basis_2][data['EXPNUM_R'] == expnum])
##    plt.scatter(x, y, edgecolor='none', s=3, label='EXPNUM_R = {}'.format(expnum))
##    plt.xlim(0.5, -0.5)
##    plt.ylim(-0.5, 0.5)
##    plt.gca().set_aspect('equal')
##    plt.xlabel(r'$\Delta \alpha$ (deg)')
##    plt.ylabel(r'$\Delta \delta$ (deg)')
##    plt.legend(loc='upper left')
##    plt.title('Stars at (RA, Dec) = ({}, {})'.format(ra_select, dec_select))
##    plt.savefig('v2_scatter_test_expnum-r_{}.png'.format(expnum))
##    plt.close()
#proj = ugali.utils.projector.Projector(ra_select, dec_select)
#x, y = proj.sphereToImage(data[basis_1], data[basis_2])
#plt.scatter(x, y, edgecolor='none', s=3)
#g_ra = [178.216, 178.233, 177.809]
#g_dec = [-41.8225, -41.806, -41.841]
#r_ra = [176.82, 176.832, 176.844, 176.856, 177.809, 177.038]
#r_dec = [-41.738, -41.7227, -41.7074, -41.6922, -41.841, -41.294]
#g_x, g_y = proj.sphereToImage(g_ra, g_dec)
#r_x, r_y = proj.sphereToImage(r_ra, r_dec)
#plt.scatter(g_x, g_y, edgecolor='none', c='g', s=10, label='g')
#plt.scatter(r_x, r_y, edgecolor='none', c='r', s=10, label='r')
#plt.xlim(1.0, -1.0)
#plt.ylim(-1.0, 1.0)
#plt.gca().set_aspect('equal')
#plt.xlabel(r'$\Delta \alpha$ (deg)')
#plt.ylabel(r'$\Delta \delta$ (deg)')
#plt.legend(loc='upper left')
#plt.title('Stars at (RA, Dec) = ({}, {})'.format(ra_select, dec_select))
#plt.savefig('v2_scatter_test.png')
#plt.close()

## CMD Plot
#angsep_1 = ugali.utils.projector.angsep(ra_select-0.5, dec_select, data[basis_1], data[basis_2])
#angsep_2 = ugali.utils.projector.angsep(ra_select+0.5, dec_select, data[basis_1], data[basis_2])
#g_radius = 0
#annulus_1 = (angsep_1 > g_radius) & (angsep_1 < 1.)
#annulus_2 = (angsep_2 > g_radius) & (angsep_2 < 1.)
#
## Plot background objects
#plt.scatter(data[mag_g][annulus_1] - data[mag_r][annulus_1], data[mag_g][annulus_1], c='r', alpha=0.5, edgecolor='none', s=1, label='RA-0.5')
#plt.scatter(data[mag_g][annulus_2] - data[mag_r][annulus_2], data[mag_g][annulus_2], c='g', alpha=0.5, edgecolor='none', s=1, label='RA+0.5')
#plt.axvline(x=np.mean(data[mag_g][annulus_1] - data[mag_r][annulus_1]), c='r', label='mean = {}'.format(np.mean(data[mag_g][annulus_1] - data[mag_r][annulus_1])))
#plt.axvline(x=np.mean(data[mag_g][annulus_2] - data[mag_r][annulus_2]), c='g', label='mean = {}'.format(np.mean(data[mag_g][annulus_2] - data[mag_r][annulus_2])))
#
#plt.axis([-0.5, 1, 16, mag_max])
#plt.gca().invert_yaxis()
#plt.gca().set_aspect(1./4.)
#plt.xlabel('g-r (mag)')
#plt.ylabel('g (mag)')
#plt.legend(loc='upper left')
#plt.title('CMD at (RA, Dec) = ({}, {})'.format(ra_select, dec_select))
#plt.savefig('cmd_test.png')
#plt.close()

# Density plot
proj = ugali.utils.projector.Projector(ra_select, dec_select)
x, y = proj.sphereToImage(data[basis_1], data[basis_2])

bound = 0.5 #1.
steps = 100.
bins = np.linspace(-bound, bound, steps)

signal = np.histogram2d(x, y, bins=[bins, bins])[0]

g_radius = 0.5
sigma = 0.01 * (0.25 * np.arctan(0.25*g_radius*60. - 1.5) + 1.3)

convolution = scipy.ndimage.filters.gaussian_filter(signal, sigma/(bound/steps))
plt.pcolormesh(bins, bins, convolution.T, cmap='Greys')

plt.xlim(bound, -bound)
plt.ylim(-bound, bound)
plt.gca().set_aspect('equal')
plt.xlabel(r'$\Delta \alpha$ (deg)')
plt.ylabel(r'$\Delta \delta$ (deg)')
plt.colorbar()
plt.title('Density (MAG < 90 in g, i) at (RA, Dec) = ({}, {})'.format(ra_select, dec_select))
plt.savefig('density_g-i.png')
plt.close()
