#!/usr/bin/env python
"""
Compile candidate list from results_dir
"""
__author__ = "Sidney Mau"

import glob
import yaml

from astropy.io import fits
# import fitsio
import numpy as np
import csv

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

candidate_list = cfg[cfg['data']]['candidate_list']


# Parse results from results_dir into a list of values
results = []
for file in glob.glob('{}/*.csv'.format(cfg['output']['results_dir'])):
    with open(file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            results.append([float(val) for val in row])

data = np.asarray(results)

# Create fits columns
c1 = fits.Column(name='SIG',     format='E', array=data[:,0])
c2 = fits.Column(name='RA',      format='E', array=data[:,1])
c3 = fits.Column(name='DEC',     format='E', array=data[:,2])
c4 = fits.Column(name='MODULUS', format='E', array=data[:,3])
c5 = fits.Column(name='r',       format='E', array=data[:,4])

# Write fits output
t = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5])
t.writeto(candidate_list)
