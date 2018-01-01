#!/usr/bin/env python
"""
Compile candidate list from results_dir
"""
__author__ = "Sidney Mau"

import glob
import yaml

from astropy.io import fits
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
    csvfile.close()

data = np.asarray(results)
dt = np.dtype([('SIG', '>f4'), ('RA', '>f4'), ('DEC', '>f4'), ('MODULUS', '>f4'), ('r', '>f4')])
data = data.astype(dt)

# Create fits columns
c0 = fits.Column(name='SIG',     format='E', array=data[:,0])
c1 = fits.Column(name='RA',      format='E', array=data[:,1])
c2 = fits.Column(name='DEC',     format='E', array=data[:,2])
c3 = fits.Column(name='MODULUS', format='E', array=data[:,3])
c4 = fits.Column(name='r',       format='E', array=data[:,4])
#c5 = fits.Column(name='ASSOC',   format='E', array=data[:,5])
#c6 = fits.Column(name='ANGSEP',  format='E', array=data[:,6])

# Write fits output
t = fits.BinTableHDU.from_columns([c0, c1, c2, c3, c4])
t.writeto(candidate_list, overwrite=True)

#from fitsio import FITS
#
#fits = FITS(candidate_list, 'rw')
##names = ['SIG', 'RA', 'DEC', 'MODULUS', 'r']
##fits.write(data, names=names)
#fits.write(data)
#fits.close()
