#!/usr/bin/env python
"""
Classy data
"""
__author__ = "Sid Mau"

# Python libraries
import os
import glob
import numpy as np
import numpy.lib.recfunctions
import healpy as hp
import fitsio as fits
import scipy.ndimage
import itertools

# Ugali libraries
import ugali.utils.healpix
import ugali.utils.projector
import ugali.isochrone

# Simple libraries
from simple.objects.data import Data

DES = Data(survey = "DES",
           nside  = 32,
           dirname = "/home/s1/kadrlica/projects/y3a2/data/gold/v2.0/healpix",
           fracdet = "y3a2_griz_o.4096_t.32768_coverfoot_EQU_decompressed.fits",
           band_1 = "G",
           band_2 = "R",
           mag = "SOF_PSF_MAG_{}",
           mag_err = "SOF_PSF_MAG_ERR_{}",
           mag_dered = "SOF_PSF_MAG_CORRECTED_{}",
           basis_1 = "RA",
           basis_2 = "DEC",
           mag_max = 24.5)

import pdb;pdb.set_trace()

