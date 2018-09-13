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
    def __init__(self, result_dir, result_file):
        self.result_dir  = result_dir
        self.result_file = result_file

        # Initialize result arrays
        self.ra_peak_array          = []
        self.dec_peak_array         = [] 
        self.r_peak_array           = []
        self.sig_peak_array         = []
        self.distance_modulus_array = []
        self.mc_source_id_array     = []
        self.n_obs_peak_array       = []
        self.n_obs_half_peak_array  = []
        self.n_model_peak_array     = []
        self.mc_source_id_array     = []

    def append_results(self, ra_peaks, dec_peaks, r_peaks, sig_peaks, distance_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks):
        """
        Append results from Search.search_by_distance().
        """
        # Append results to result arrays
        self.ra_peak_array.append(ra_peaks)
        self.dec_peak_array.append(dec_peaks)
        self.r_peak_array.append(r_peaks)
        self.sig_peak_array.append(sig_peaks)
        self.distance_modulus_array.append(dist_moduli)
        self.n_obs_peak_array.append(n_obs_peaks)
        self.n_obs_half_peak_array.append(n_obs_half_peaks)
        self.n_model_peak_array.append(n_model_peaks)
        self.mc_source_id_array.append(np.tile(0, len(sig_peaks)))

        return

    def concatenate_results(self):
        """
        Concatenate result arrays.
        """
        # Concatenate result arrays
        self.ra_peak_array          = np.concatenate(self.ra_peak_array)
        self.dec_peak_array         = np.concatenate(self.dec_peak_array)
        self.r_peak_array           = np.concatenate(self.r_peak_array)
        self.sig_peak_array         = np.concatenate(self.sig_peak_array)
        self.distance_modulus_array = np.concatenate(self.distance_modulus_array)
        self.n_obs_peak_array       = np.concatenate(self.n_obs_peak_array)
        self.n_obs_half_peak_array  = np.concatenate(self.n_obs_half_peak_array)
        self.n_model_peak_array     = np.concatenate(self.n_model_peak_array)
        self.mc_source_id_array     = np.concatenate(self.mc_source_id_array)

        return

    def sort_results(self):
        """
        Sort result arrays.
        """
        # Sort peaks according to significance
        index_sort                  = np.argsort(self.sig_peak_array)[::-1]
        self.ra_peak_array          = self.ra_peak_array[index_sort]
        self.dec_peak_array         = self.dec_peak_array[index_sort]
        self.r_peak_array           = self.r_peak_array[index_sort]
        self.sig_peak_array         = self.sig_peak_array[index_sort]
        self.distance_modulus_array = self.distance_modulus_array[index_sort]
        self.n_obs_peak_array       = self.n_obs_peak_array[index_sort]
        self.n_obs_half_peak_array  = self.n_obs_half_peak_array[index_sort]
        self.n_model_peak_array     = self.n_model_peak_array[index_sort]
        self.mc_source_id_array     = self.mc_source_id_array[index_sort]

        return

    #def write_output(self, ra_peak_array, dec_peak_array, r_peak_array, distance_modulus_array, 
    #                n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, 
    #                sig_peak_array, mc_source_id_array,  outfile):
    #    writer = open(outfile, 'a') # don't overwrite; append
    #    for ii in range(0, len(sig_peak_array)):
    #        # SIG, RA, DEC, MODULUS, r, n_obs, n_model, mc_source_id
    #        writer.write('{:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}\n'.format(sig_peak_array[ii], 
    #                                                                                                                         ra_peak_array[ii], 
    #                                                                                                                         dec_peak_array[ii], 
    #                                                                                                                         distance_modulus_array[ii], 
    #                                                                                                                         r_peak_array[ii],
    #                                                                                                                         n_obs_peak_array[ii],
    #                                                                                                                         n_obs_half_peak_array[ii],
    #                                                                                                                         n_model_peak_array[ii],
    #                                                                                                                         mc_source_id_array[ii]))
    
    #def make_list(self):
    #   #

    #def read_output(self, results_dir):
    #    infiles = sorted(glob.glob(results_dir + '/*.csv'))
    #    results = np.concatenate([np.genfromtxt(infile, delimiter=',', names=['SIG', 'RA', 'DEC', 'DISTANCE_MODULUS', 'R_PEAK', 'N_OBS_PEAK', 'N_OBS_HALF_PEAK', 'N_MODEL_PEAK', 'MC_SOURCE_ID']) for infile in infiles])
    #    return results
