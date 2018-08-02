#!/usr/bin/env python
"""
Pseudo-fracdet map
"""
__author__ = "Sidney Mau"

import os
import glob
import yaml
import numpy as np
import healpy as hp
import fitsio as fits


with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']


############################################################

infiles = glob.glob ('{}/*.fits'.format(datadir))

nside = 2048
pix = []
for infile in infiles:
    print('loading {}'.format(infile))
    data = fits.read(infile, columns=[basis_1,basis_2])
    p = hp.ang2pix(nside, data[basis_1], data[basis_2], lonlat=True)
    pix.append(np.unique(p))

print('Constructing map')
pix = np.concatenate(pix)
pix = np.unique(pix)
coverage_map = np.tile(hp.UNSEEN, hp.nside2npix(nside))
coverage_map[pix] = 1

print('Writing output')
result = '{}_pseudo_fracdet.fits.gz'.format(survey)
hp.write_map(result, coverage_map)
