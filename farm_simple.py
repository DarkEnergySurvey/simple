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
import healpy as hp
import numpy as np
import fitsio as fits

import ugali.utils.healpix

import yaml

############################################################

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']
    simple_dir = cfg['setup']['simple_dir']
    jobs = cfg['batch']['jobs']
    nside = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']
    mode = cfg[survey]['mode']
    sim_catalog = cfg[survey]['sim_catalog']
    sim_population = cfg[survey]['sim_population']

results_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'])
if not os.path.exists(results_dir):
    os.mkdir(results_dir)

log_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'], cfg['output']['log_dir'])
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

############################################################

def submitJob(ra, dec, pix, mc_source_id, mode):
    if (mode == 0):
        logfile = '{}/results_nside_{}_{}.log'.format(log_dir, nside, pix)
    elif (mode == 1):
        logfile = '{}/results_mc_source_id_{}.log'.format(log_dir, mc_source_id) # all values in mc_source_id_array should be the same
    batch = 'csub -n {} -o {} '.format(jobs, logfile)
    command = 'python {}/search_algorithm.py {:0.2f} {:0.2f} {:0.2f}'.format(simple_dir, ra, dec, mc_source_id)
    command_queue = batch + command
    print(command_queue)
    #os.system('./' + command) # Run locally
    os.system(command_queue) # Submit to queue

    return

def executeJob(ra, dec, pix, mc_source_id, mode):
    if (mode == 0):
        logfile = 'results_nside_{}_{}.log'.format(nside, pix)
    elif (mode == 1):
        logfile = 'results_mc_source_id_{}.log'.format(mc_source_id) # all values in mc_source_id_array should be the same
    command = 'python {}/search_algorithm.py {:0.2f} {:0.2f} {:0.2f} >> {}/{}'.format(simple_dir, ra, dec, mc_source_id, log_dir, logfile)
    print(command)
    os.system(command) # Submit to queue # TODO use subprocess instead

    return

############################################################

print('mode: {}...'.format(mode))

if (mode == 0): # real
    infiles = glob.glob ('{}/*.fits'.format(datadir))
    
    print('Pixelizing...')
    pix_nside = [] # Equatorial coordinates, RING ordering scheme
    for infile in infiles:
        pix_nside.append(int(infile.split('.fits')[0].split('_')[-1]))

    for ii in range(0, len(pix_nside)):
        ra, dec = ugali.utils.healpix.pixToAng(nside, pix_nside[ii])
    
        submitJob(ra, dec, pix_nside[ii], 0) # TODO: mc_source_id (0 for real)
        print('({}/{})').format(ii, len(pix_nside))
    
elif (mode == 1): # real+sim
    sim_pop = fits.read(sim_population)
    for sim in sim_pop:
        ra, dec, mc_source_id = sim[basis_1], sim[basis_2], sim['MC_SOURCE_ID']
        pix = hp.ang2pix(nside, ra, dec, lonlat=True)

        submitJob(ra, dec, pix, mc_source_id, mode)

        #batch = 'csub -n {} -o {} '.format(jobs, logfile)
