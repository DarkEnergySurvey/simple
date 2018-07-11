#!/usr/bin/env python
"""
Perform simple binning search
"""
__author__ = "Sidney Mau"

import os
import time
import subprocess
import glob
import healpy as hp
import numpy as np
import fitsio as fits

import ugali.utils.healpix

import yaml

from multiprocessing import Pool

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
    sim_population = cfg[survey]['sim_population']
    sim_dir = cfg[survey]['sim_dir']
    object_list = cfg[survey]['object_list']

results_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'])
if not os.path.exists(results_dir):
    os.mkdir(results_dir)

log_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'], cfg['output']['log_dir'])
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

############################################################

def submitJob(ra, dec, pix, mc_source_id, mode, **population_file):
    if (mode == 0):
        outfile = '{}/results_nside_{}_{}.csv'.format(results_dir, nside, pix)
        logfile = '{}/results_nside_{}_{}.log'.format(log_dir, nside, pix)
        #command = 'python {}/search_algorithm.py {:0.2f} {:0.2f} {:0.2f} {} {}'.format(simple_dir, ra, dec, mc_source_id, outfile, logfile)
    elif (mode == 1):
        outfile = '{}/results_mc_source_id_{}.csv'.format(results_dir, mc_source_id) # all values in mc_source_id_array should be the same
        logfile = '{}/results_mc_source_id_{}.log'.format(log_dir, mc_source_id) # all values in mc_source_id_array should be the same
        #command = 'python {}/search_algorithm.py {:0.2f} {:0.2f} {:0.2f} {} {} {}'.format(simple_dir, ra, dec, mc_source_id, outfile, logfile, population_file)
    elif (mode == 2):
        outfile = '{}/results_mc_source_id_{}.csv'.format(results_dir, mc_source_id) # all values in mc_source_id_array should be the same
        logfile = '{}/results_mc_source_id_{}.log'.format(log_dir, mc_source_id) # all values in mc_source_id_array should be the same
    batch = 'csub -n {} -o {} '.format(jobs, logfile)
    command = 'python {}/search_algorithm.py {:0.2f} {:0.2f} {:0.2f} {} {}'.format(simple_dir, ra, dec, mc_source_id, outfile, logfile)
    command_queue = batch + command

    print(command_queue)
    os.system(command_queue) # Submit to queue

    return

#############################################################

print('mode: {}...'.format(mode))

if (mode == 0): # real
    infiles = glob.glob ('{}/*.fits'.format(datadir))
    
    print('Pixelizing...')
    pix_nside = [] # Equatorial coordinates, RING ordering scheme
    for infile in infiles:
        pix_nside.append(int(infile.split('.fits')[0].split('_')[-1]))

    #for ii in range(0, len(pix_nside)):
    for ii in range(0, 19):
        ra, dec = ugali.utils.healpix.pixToAng(nside, pix_nside[ii])
    
        submitJob(ra, dec, pix_nside[ii], 0, mode) # TODO: mc_source_id (0 for real)
        print('({}/{})').format(ii, len(pix_nside))
    
elif (mode == 1): # real+sim
    sim_pop = fits.read(sim_population)
    for sim in sim_pop[:]:
        #if (sim['DIFFICULTY'] == 0): # check to make sure simulated object has stars
        ra, dec, mc_source_id = sim[basis_1], sim[basis_2], sim['MC_SOURCE_ID']
        pix = hp.ang2pix(nside, ra, dec, lonlat=True)

        submitJob(ra, dec, pix, mc_source_id, mode)

elif (mode == 2): # real objects
    #sim_pop = fits.read(object_list)
    sim_pop = np.genfromtxt(object_list, delimiter=',', names=['RA', 'DEC', 'DISTANCE_MODULUS', 'MC_SOURCE_ID', 'NAME'])[1:]
    for sim in sim_pop[:]:
        #if (sim['DIFFICULTY'] == 0): # check to make sure simulated object has stars
        ra, dec, mc_source_id = sim[basis_1], sim[basis_2], sim['MC_SOURCE_ID']
        pix = hp.ang2pix(nside, ra, dec, lonlat=True)

        submitJob(ra, dec, pix, mc_source_id, mode)
