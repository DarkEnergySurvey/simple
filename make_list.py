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

candidate_list_csv = 'candidate_list.csv'
candidate_list = cfg[cfg['data']]['candidate_list']


# intermediate csv write before producing fits
# this step should be removed
results_file = open(candidate_list_csv, 'w')
results_file.write('SIG, RA, DEC, MODULUS, r\n')
for file in glob.glob('{}/*.csv'.format(cfg['output']['results_dir'])):
    writer = open(file, 'r')
    for line in writer:
        results_file.write(line)
    writer.close()
results_file.close()


# Write candidate_list as .fits file
# Currently following procedure from http://docs.astropy.org/en/stable/io/fits/
# Should move away from astropy.io.fits to fitsio

data = []
with open(candidate_list_csv, 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        data.append([val for val in row])

data = data[1:] # remove header
data = np.asarray(data) # convert to numpy.ndarray
data = data.astype(float) # convert elements from strings to floats

c1 = fits.Column(name='SIG',     format='E', array=data[:,0])
c2 = fits.Column(name='RA',      format='E', array=data[:,1])
c3 = fits.Column(name='DEC',     format='E', array=data[:,2])
c4 = fits.Column(name='MODULUS', format='E', array=data[:,3])
c5 = fits.Column(name='r',       format='E', array=data[:,4])

t = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5])
t.writeto(candidate_list)

