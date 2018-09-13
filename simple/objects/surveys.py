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

# Currently unused; likely to be removed.

class DES(Data):
    """
    Class for working with DES data.
    """
    def __init__(self):
        self.survey    = 'DES'
        self.nside     = 32
        self.datadir   = '/home/s1/kadrlica/projects/y3a2/data/gold/v2.0/healpix'
        self.fracdet   = 'y3a2_griz_o.4096_t.32768_coverfoot_EQU_decompressed.fits'
        self.band_1    = 'G'
        self.band_2    = 'R"'
        self.mag       = 'SOF_PSF_MAG_{}'
        self.mag_err   = 'SOF_PSF_MAG_ERR_{}'
        self.mag_dered = 'SOF_PSF_MAG_CORRECTED_{}'
        self.basis_1   = 'RA'
        self.basis_2   = 'DEC'
        self.mag_max   = 24.5

    def quality_mask(self, data):
        """
        Return cut on data quality.
        """
        return (data[self.mag_1] < self.mag_max)

    def star_mask(self, data):
        """
        Return cut on stellar-like objects.
        """
        return ((data['EXTENDED_CLASS_MASH_SOF'] >= 0)\
               &(data['EXTENDED_CLASS_MASH_SOF'] <= 2))

    def galaxy_mask(self, data):
        """
        Return cut on galaxy-like objects.
        """
        return (data['EXTENDED_CLASS_MASH_SOF'] > 2)

    def color_mask(self, data):
        """
        Return cut on blue objects.
        """
        return ((data[self.mag_1] - data[self.mag_2]) < 0.4) # 0.2

import pdb;pdb.set_trace()

