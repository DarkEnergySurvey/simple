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

results_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'])
if not os.path.exists(results_dir):
    os.mkdir(results_dir)

log_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'], cfg['output']['log_dir'])
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

############################################################

# Main

sim_pop = fits.read(sim_population)
sub_sim_pop_list = [sim_pop[x:x+100] for x in range(0, len(sim_pop), 100)]

for sub_sim_list in sub_sim_pop_list:
    mc_first = sub_sim_list[0]['MC_SOURCE_ID']
    mc_last = sub_sim_list[-1]['MC_SOURCE_ID']
    outfile = '{}/results_mc_source_id_{}-{}.csv'.format(results_dir, mc_first, mc_last)
    logfile = '{}/results_mc_source_id_{}-{}.log'.format(log_dir, mc_first, mc_last)

    batch = 'csub -n {} -o {} '.format(jobs, logfile)

    command = 'python {}/batch_submit.py {:0.2f} {:0.2f} {:0.2f} {:0.2f}'.format(simple_dir, mc_first, mc_last, outfile, logfile)
    command_queue = batch + command

    print(command_queue)
    os.system(command_queue) # Submit to queue
