#!/usr/bin/env python
"""
Create simple binner style plots for ugali or simple candidate lists
"""
__author__ = "Sidney Mau"

import os
import sys
import time
import subprocess
import glob
import pyfits
import healpy
import numpy
import numpy as np

import ugali.utils.healpix
import fitsio as fits

import yaml

############################################################

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    simple_dir = cfg['setup']['simple_dir']
    jobs = cfg['batch']['jobs']
    candidate_list = cfg[survey]['candidate_list']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']

save_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'])
if not os.path.exists(save_dir):
    os.mkdir(save_dir)

log_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'], cfg['output']['log_dir'])
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

try:
    sig_cut = float(sys.argv[1])
except:
    sig_cut = 5.5

print('Plotting hotspots with sig > {}'.format(sig_cut))

candidate_list = fits.read(candidate_list)
try: # simple
    candidate_list = candidate_list[candidate_list['SIG'] > sig_cut]
except: # ugali
    candidate_list = candidate_list[candidate_list['TS'] > 25]

print('{} candidates found...').format(len(candidate_list))

############################################################

for candidate in candidate_list:
    try: # simple
        sig = round(candidate['SIG'], 2)
    except: # ugali
        sig = round(candidate['TS'], 2)
    ra      = round(candidate[basis_1], 2)
    dec     = round(candidate[basis_2], 2)
    mod     = round(candidate['MODULUS'], 2)
    mc_source_id = round(candidate['MC_SOURCE_ID'], 2)
    if 'N_MODEL' in candidate_list.dtype.names:
        field_density = round(candidate['N_MODEL'] / (np.pi * (candidate['r'] * 60.)**2), 4) # field density (arcmin^-2)
    
    logfile = '{}/candidate_{}_{}.log'.format(log_dir, ra, dec)
    batch = 'csub -n {} -o {} '.format(jobs, logfile)
    if 'N_MODEL' in candidate_list.dtype.names:
        command = 'python {}/make_plot.py {:0.2f} {:0.2f} {:0.2f} {:0.2f} {:0.2f} {:0.4f}'.format(simple_dir, ra, dec, mod, sig, mc_source_id, field_density)
    else:
        command = 'python {}/make_plot.py {:0.2f} {:0.2f} {:0.2f} {:0.2f} {:0.2f}'.format(simple_dir, ra, dec, mod, sig, mc_source_id)
    command_queue = batch + command

    #print(command)
    #os.system(command)
    print(command_queue)
    os.system(command_queue) # Submit to queue
