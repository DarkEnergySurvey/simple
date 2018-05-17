#!/usr/bin/env python
"""
Perform simple binning search
"""
__author__ = "Sidney Mau"

import os
import time
import subprocess
import healpy as hp
import numpy as np
import fitsio as fits

import yaml

import batch_submit

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

def execute_batch_jobs(sub_sim_list, outfile, logfile):
    for sim in sub_sim_list:
        if (sim['DIFFICULTY'] == 0):
            ra, dec, mc_source_id = sim[basis_1], sim[basis_2], sim['MC_SOURCE_ID']
            pix = hp.ang2pix(nside, ra, dec, lonlat=True)
    
            command = 'python {}/search_algorithm.py {:0.2f} {:0.2f} {:0.2f} {} >> {}'.format(simple_dir, ra, dec, mc_source_id, outfile, logfile)
            print(command)
            os.system(command) # Submit to queue

    return

############################################################

# Main

sim_pop = fits.read(sim_population)
sub_sim_pop_list = [sim_pop[x:x+100] for x in range(0, len(sim_pop), 100)]

for sub_sim_list in sub_sim_pop_list:
    outfile = '{}/results_mc_source_id_{}-{}.csv'.format(results_dir, sub_sim_list[0]['MC_SOURCE_ID'], sub_sim_list[-1]['MC_SOURCE_ID'])
    logfile = '{}/results_mc_source_id_{}-{}.log'.format(log_dir, sub_sim_list[0]['MC_SOURCE_ID'], sub_sim_list[-1]['MC_SOURCE_ID'])

    batch = 'csub -n {} -o {} '.format(jobs, logfile)

    command = 'python {}/batch_submit.py {:0.2f} {:0.2f} {:0.2f}'.format(simple_dir, outfile, logfile)
    command_queue = batch + command
    print(command_queue)
    #os.system('./' + command) # Run locally
    os.system(command_queue) # Submit to queue
