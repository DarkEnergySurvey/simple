#!/usr/bin/env python
"""
Batch submit
"""
__author__ = "Sidney Mau"

import sys
import os
import subprocess
import healpy as hp
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
    sim_catalog = cfg[survey]['sim_catalog']
    sim_population = cfg[survey]['sim_population']

results_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'])
if not os.path.exists(results_dir):
    os.mkdir(results_dir)

log_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'], cfg['output']['log_dir'])
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

############################################################

try:
    mc_first, mc_last, outfile, logfile = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3], sys.argv[4]
except:
    sys.exit('ERROR!')

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

# Main:

sim_pop = fits.read(sim_population)
sub_sim_list = sim_pop[(sim_pop['MC_SOURCE_ID'] >= mc_first) & (sim_pop['MC_SOURCE_ID'] <= mc_last)]

execute_batch_jobs(sub_sim_list, outfile, logfile)
