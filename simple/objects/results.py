#!/usr/bin/env python
"""
Classy results
"""
__author__ = "Sid Mau"

# Python libraries
import os
import glob
import yaml
import numpy as np
import numpy.lib.recfunctions
import fitsio as fits
import itertools

# Ugali libraries
import ugali.utils.healpix
import ugali.utils.projector
import ugali.isochrone

########################################################################

# not usable yet
class Result:
    """
    Class object to write, store, and read results.
    """
    def __init__(self, rdir, rfile):
        self.rdir  = rdir
        self.rfile = rfile

    def write_output(self, ra_peak_array, dec_peak_array, r_peak_array, distance_modulus_array, 
                    n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, 
                    sig_peak_array, mc_source_id_array,  outfile):
        writer = open(outfile, 'a') # don't overwrite; append
        for ii in range(0, len(sig_peak_array)):
            # SIG, RA, DEC, MODULUS, r, n_obs, n_model, mc_source_id
            writer.write('{:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}\n'.format(sig_peak_array[ii], 
                                                                                                                             ra_peak_array[ii], 
                                                                                                                             dec_peak_array[ii], 
                                                                                                                             distance_modulus_array[ii], 
                                                                                                                             r_peak_array[ii],
                                                                                                                             n_obs_peak_array[ii],
                                                                                                                             n_obs_half_peak_array[ii],
                                                                                                                             n_model_peak_array[ii],
                                                                                                                             mc_source_id_array[ii]))
    
    def make_list(self):
       #

    def read_output(self, results_dir):
        infiles = sorted(glob.glob(results_dir + '/*.csv'))
        results = np.concatenate([np.genfromtxt(infile, delimiter=',', names=['SIG', 'RA', 'DEC', 'DISTANCE_MODULUS', 'R_PEAK', 'N_OBS_PEAK', 'N_OBS_HALF_PEAK', 'N_MODEL_PEAK', 'MC_SOURCE_ID']) for infile in infiles])
        return results
