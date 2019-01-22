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
import fitsio

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']

    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']
    candidate_list = cfg[survey]['candidate_list']


# Parse results from results_dir into a list of values
results = []
for file in glob.glob('{}/*.csv'.format(cfg['output']['results_dir'])):
    with open(file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            results.append([float(val) for val in row])
    csvfile.close()

data = np.asarray(results)
#data = data[np.unique([data[i][-1] for i in range(len(data))], return_index=True)[1]]

# Create fits columns
c0 = fits.Column(name='SIG',          format='E', array=data[:,0])
c1 = fits.Column(name=basis_1,        format='E', array=data[:,1])
c2 = fits.Column(name=basis_2,        format='E', array=data[:,2])
c3 = fits.Column(name='MODULUS',      format='E', array=data[:,3])
c4 = fits.Column(name='r',            format='E', array=data[:,4])
c5 = fits.Column(name='N_OBS',        format='E', array=data[:,5])
c6 = fits.Column(name='N_OBS_HALF',   format='E', array=data[:,6])
c7 = fits.Column(name='N_MODEL',      format='E', array=data[:,7])
c8 = fits.Column(name='MC_SOURCE_ID', format='E', array=data[:,8])

# Write fits output
t = fits.BinTableHDU.from_columns([c0, c1, c2, c3, c4, c5, c6, c7, c8])
t.writeto(candidate_list, overwrite=True)

#from fitsio import FITS
#
#fits = FITS(candidate_list, 'rw')
##names = ['SIG', basis_1, basis_2, 'MODULUS', 'r']
##fits.write(data, names=names)
#fits.write(data)
#fits.close()

# Diagnostic output
data = fitsio.read(candidate_list)
print("{} hotspots found.").format(len(data))
cut_0 = (data['SIG'] > 5.5)
print("{} hotspots found with SIG > 5.5.").format(len(data[cut_0]))
cut_1 = (data['SIG'] > 10)
print("{} hotspots found with SIG > 10.").format(len(data[cut_1]))
cut_2 = (data['SIG'] > 15)
print("{} hotspots found with SIG > 15.").format(len(data[cut_2]))
cut_3 = (data['SIG'] > 20)
print("{} hotspots found with SIG > 20.").format(len(data[cut_3]))
cut_4 = (data['SIG'] > 25)
print("{} hotspots found with SIG > 25.").format(len(data[cut_4]))
cut_5 = (data['SIG'] >= 37.5)
print("{} hotspots found with SIG >= 37.55").format(len(data[cut_5]))
