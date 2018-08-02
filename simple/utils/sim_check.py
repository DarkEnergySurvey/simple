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

    candidate_list = cfg[survey]['candidate_list']
    mode = cfg[survey]['mode']
    sim_population = cfg[survey]['sim_population']
    sim_dir = cfg[survey]['sim_dir']

    fracdet_map = cfg[survey]['fracdet']
    
    mag_g = cfg[survey]['mag_g']
    mag_r = cfg[survey]['mag_r']
    mag_g_err = cfg[survey]['mag_g_err']
    mag_r_err = cfg[survey]['mag_r_err']
    
########################################################################

data = fits.read(candidate_list)
sim_pop = fits.read(sim_population)[:]

#data = data[(sim_pop['DIFFICULTY'] == 0) & (sim_pop['N_CATALOG'] > 0)]
#sim_pop = sim_pop[(sim_pop['DIFFICULTY'] == 0) & (sim_pop['N_CATALOG'] > 0)]

# Sort by MC_SOURCE_ID
data.sort(order='MC_SOURCE_ID')
sim_pop.sort(order='MC_SOURCE_ID')

# STELLAR MASS
plt.scatter(sim_pop['STELLAR_MASS'], data['SIG'], c=sim_pop['DIFFICULTY'], s=1)
plt.xlabel('STELLAR_MASS')
plt.ylabel('SIG')
plt.title('{} STELLAR_MASS'.format(survey))
plt.colorbar().set_label('DIFFICULTY')
plt.savefig('{}_sim_STELLAR_MASS.png'.format(survey))
plt.close()

# R_Physical
plt.scatter(sim_pop['R_PHYSICAL'], data['SIG'], c=sim_pop['DIFFICULTY'], s=1)
plt.xlabel('R_PHYSICAL')
plt.ylabel('SIG')
plt.title('{} R_PHYSICAL'.format(survey))
plt.colorbar().set_label('DIFFICULTY')
plt.savefig('{}_sim_R_PHYSICAL.png'.format(survey))
plt.close()

# DISTANCE
plt.scatter(sim_pop['DISTANCE'], data['SIG'], c=sim_pop['DIFFICULTY'], s=1)
plt.xlabel('DISTANCE')
plt.ylabel('SIG')
plt.title('{} DISTANCE'.format(survey))
plt.colorbar().set_label('DIFFICULTY')
plt.savefig('{}_sim_DISTANCE.png'.format(survey))
plt.close()


# N_G22
plt.scatter(sim_pop['N_G22'], data['SIG'], c=sim_pop['DIFFICULTY'], s=1)
plt.xscale('log')
plt.xlabel('log(N_G22)')
plt.ylabel('SIG')
plt.title('{} N_G22'.format(survey))
plt.colorbar().set_label('DIFFICULTY')
plt.savefig('{}_sim_N_G22.png'.format(survey))
plt.close()

# N_G22 vs SURFACE_BRIGHTNESS
plt.scatter(sim_pop['N_G22'], sim_pop['SURFACE_BRIGHTNESS'], c=data['SIG'], s=1)
plt.xscale('log')
plt.xlim(1, np.max(sim_pop['N_G22']))
plt.xlabel('log(N_G22)')
plt.ylabel('SURFACE_BRIGHTNESS')
plt.title('{} Surface Brightness vs N_G22'.format(survey))
plt.colorbar().set_label('SIG')
plt.savefig('{}_sim_N_G22-SURFACE_BRIGHTNESS.png'.format(survey))
plt.close()
