#!/usr/bin/env python
"""
Arrange and produce plots
"""
__author__ = "Sidney Mau"

# Set the backend first!
import matplotlib
matplotlib.use('Agg')

import sys
import os
import yaml
import numpy as np
import matplotlib.pylab as plt

import ugali.utils.projector
import ugali.candidate.associate

import simple.plotting.diagnostic_plots

print(matplotlib.get_backend())

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']

save_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'])
if not os.path.exists(save_dir):
    os.mkdir(save_dir)

#################################################################

try:
    targ_ra, targ_dec, mod, sig, mc_source_id, field_density, = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]) 
except:
    try:
        targ_ra, targ_dec, mod, sig, mc_source_id, = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
        field_density = None
    except:
        sys.exit('ERROR! Coordinates not given in correct format.')

data, iso, g_radius, nbhd = simple.plotting.diagnostic_plots.analysis(targ_ra, targ_dec, mod, mc_source_id)

print('Making diagnostic plots for ({}, {}) = ({}, {})...'.format(basis_1, basis_2, targ_ra, targ_dec))

file_name = 'candidate_{:0.2f}_{:0.2f}'.format(targ_ra, targ_dec)

fig, axs = plt.subplots(1, 5, figsize=(24, 4))
fig.subplots_adjust(wspace=0.5, hspace=0.5)

simple.plotting.diagnostic_plots.star_plot(axs[0], targ_ra, targ_dec, data, iso, g_radius, nbhd)
simple.plotting.diagnostic_plots.density_plot(axs[1], targ_ra, targ_dec, data, iso, g_radius, nbhd, 'stars')
simple.plotting.diagnostic_plots.density_plot(axs[2], targ_ra, targ_dec, data, iso, g_radius, nbhd, 'galaxies')
simple.plotting.diagnostic_plots.cm_plot(axs[3], targ_ra, targ_dec, data, iso, g_radius, nbhd, 'stars')
simple.plotting.diagnostic_plots.hess_plot(axs[4], targ_ra, targ_dec, data, iso, g_radius, nbhd)

plt.savefig(save_dir+'/'+file_name+'.png', bbox_inches='tight')
plt.close()
