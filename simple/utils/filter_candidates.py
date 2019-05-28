#!/usr/bin/env python

import numpy as np
import healpy
import astropy.io.fits as pyfits
import pylab
import argparse
import yaml

import ugali.utils.projector
import ugali.utils.healpix

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    simple_dir = cfg['setup']['simple_dir']
    save_dir = cfg['output']['save_dir']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']

def make_cut(candidate_list):
    survey = 'ps1'
    infile = candidate_list
    r = pyfits.open(infile)
    d_original = r[1].data
    r.close()
    
    # Minimum significance threshold
    d = d_original#[d_original['SIG'] > 5.] [d_original['TS'] > 25?]
    
    SIG = 'SIG'
    
    ### Define and apply cuts
    
    # Consolidate nearby peaks, iterave approach
    done = False
    while not done:
        match_1, match_2, angsep = ugali.utils.projector.match(d['ra'], d['dec'], d['ra'], d['dec'], tol=0.5, nnearest=2)
        index_exclude = np.where(d[SIG][match_1] > d[SIG][match_2], match_2, match_1)
        if len(index_exclude) == 0:
            done = True
        cut_consolidate = np.tile(True, len(d))
        cut_consolidate[index_exclude] = False
        d = d[cut_consolidate]
    
    # Geometric cuts
    pix = ugali.utils.healpix.angToPix(4096, d['ra'], d['dec'], nest=True)
    mask = ugali.utils.healpix.read_map('healpix_mask_{}.fits.gz'.format(survey), nest=True)
    
    cut_footprint = np.where(mask[pix] & 0b10000, False, True)
    cut_ebv = np.where(mask[pix] & 0b0001, False, True)
    cut_associate = np.where(mask[pix] & 0b0010, False, True)
    cut_dwarfs = np.where(mask[pix] & 0b0100, False, True)
    cut_bsc = np.where(mask[pix] & 0b1000, False, True)
    
    # Other cuts (modulus, size, shape)
    if survey == 'ps1':
        cut_modulus = (d['MODULUS'] < 21.75)
    elif survey == 'des':
        cut_modulus = (d['MODULUS'] < 23.5)
    
    cut_size = (d['r'] >= 0.020)
    model = d['N_MODEL'] / (np.pi * d['r']**2)
    core = d['N_OBS_HALF'] / (np.pi * (0.5 * d['r'])**2)
    full = d['N_OBS'] / (np.pi * d['r']**2)
    ratio = (core - model) / (full - model)
    cut_core = (ratio > 1.)
    
    # Significance cut
    if survey == 'ps1':
        min_sig = 6.
    elif survey == 'des':
        min_sig = 7.
    cut_sig = d[SIG] > min_sig
    
    cut_bulk = cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc
    cut_bulk = cut_bulk & cut_core & cut_size
    
    cut_final = cut_bulk & cut_dwarfs & cut_sig & cut_consolidate

    return(d[cut_final])
