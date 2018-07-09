#!/usr/bin/env python

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

import filters
import simple_utils
import diagnostic_plots

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

    candidate_list = cfg[survey]['candidate_list']
    mode = cfg[survey]['mode']
    sim_catalog = cfg[survey]['sim_catalog']
    sim_population = cfg[survey]['sim_population']
    sim_dir = cfg[survey]['sim_dir']

    fracdet_map = cfg[survey]['fracdet']
    
    mag_g = cfg[survey]['mag_g']
    mag_r = cfg[survey]['mag_r']
    mag_g_err = cfg[survey]['mag_g_err']
    mag_r_err = cfg[survey]['mag_r_err']
    
########################################################################

data = fits.read(candidate_list)
#sim_pop = fits.read(sim_population)[:-1] # DES
sim_pop = fits.read(sim_population)[:] # PS

#data = data[(sim_pop['DIFFICULTY'] == 0) & (sim_pop['N_CATALOG'] > 0)]
#sim_pop = sim_pop[(sim_pop['DIFFICULTY'] == 0) & (sim_pop['N_CATALOG'] > 0)]

# Sort by MC_SOURCE_ID
data.sort(order='MC_SOURCE_ID')
sim_pop.sort(order='MC_SOURCE_ID')

# STELLAR_MASS, R_PHYSICAL, DISTANCE
f, (ax1, ax2, ax3) = plt.subplots(1,3)
f.suptitle('PS1 Simulations Spot Check')

ax1.scatter(sim_pop['STELLAR_MASS'], data['SIG'], c=sim_pop['DIFFICULTY'], s=1)
ax1.set(xlabel='STELLAR_MASS')
ax1.set(ylabel='SIG')

ax2.scatter(sim_pop['R_PHYSICAL'], data['SIG'], c=sim_pop['DIFFICULTY'], s=1)
ax2.set(xlabel='R_PHYSICAL')
#ax2.set(ylabel='SIG')

ax3.scatter(sim_pop['DISTANCE'], data['SIG'], c=sim_pop['DIFFICULTY'], s=1)
ax3.set(xlabel='DISTANCE')
#ax3.set(ylabel='SIG')

#f.colorbar()
f.savefig('{}_sim_spot_checks.png'.format(survey))
plt.close()

# N_G22
plt.scatter(sim_pop['N_G22'], data['SIG'], c=sim_pop['DIFFICULTY'], s=1)

plt.xlabel('N_G22')
plt.ylabel('SIG')
plt.colorbar()

plt.savefig('{}_sim_N_G22.png'.format(survey))
plt.close()
