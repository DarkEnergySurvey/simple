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

    survey = cfg['data']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']


############################################################

infiles = glob.glob ('{}/*.fits'.format(datadir))

nside = 2048
pix = []
for infile in infiles:
    data = fits.read(infile, columns=['RA','DEC'])
    p = hp.ang2pix(nside, data['RA'], data['DEC'], lonlat=True)
    pix.append(p)
pix = np.unique(pix)
max = hp.zeros(hp.nside2npix(nside))
map[pix] = 1

result = '{}_pseudo_fracdet.csv'.format(survey)
hp.write_map(result, map)
