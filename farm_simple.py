#!/usr/bin/env python
"""
Perform simple binning search
"""
__author__ = "Sidney Mau"

import os
import time
import subprocess
import glob
import pyfits
import healpy
import numpy

import ugali.utils.healpix

import yaml

############################################################

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

simple_dir = cfg['setup']['simple_dir']

jobs = cfg['batch']['jobs']

nside = cfg[cfg['data']]['nside']
datadir = cfg[cfg['data']]['datadir']

results_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'])
if not os.path.exists(results_dir):
    os.mkdir(results_dir)

log_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'], cfg['output']['log_dir'])
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

infiles = glob.glob ('{}/*.fits'.format(datadir))

############################################################

print('Pixelizing...')
pix_nside = [] # Equatorial coordinates, RING ordering scheme
for infile in infiles:
    pix_nside.append(int(infile.split('.fits')[0].split('_')[-1]))

############################################################

for ii in range(0, len(pix_nside)):
    ra, dec = ugali.utils.healpix.pixToAng(nside, pix_nside[ii])

    print('({}/{})').format(ii, len(pix_nside))

    logfile = '{}/results_nside_{}_{}.log'.format(log_dir, nside, pix_nside[ii])
    batch = 'csub -n {} -o {} '.format(jobs, logfile)
    command = 'python {}/search_algorithm.py {:0.2} {:0.2f}'.format(simple_dir, ra, dec)
    command_queue = batch + command
    print(command_queue)
    #os.system('./' + command) # Run locally
    os.system(command_queue) # Submit to queue
