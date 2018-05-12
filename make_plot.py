#!/usr/bin/env python
"""
Arrange and produce plots
"""
__author__ = "Sidney Mau"

# Set the backend first!
import matplotlib
#matplotlib.use('Agg')

import sys
import os
import yaml
import numpy as np
import numpy
import pylab as plt
from matplotlib import gridspec

import ugali.utils.projector
import ugali.candidate.associate
import diagnostic_plots

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
    targ_ra, targ_dec, mod, sig, mc_source_id, = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
except:
    sys.exit('ERROR! Coordinates not given in correct format.')

fig = plt.figure(figsize=(20, 17))
fig.subplots_adjust(wspace=0.3, hspace=0.3)
gs = gridspec.GridSpec(3, 3)

data, iso, g_radius, nbhd = diagnostic_plots.analysis(targ_ra, targ_dec, mod, mc_source_id)

print('Making diagnostic plots for ({}, {}) = ({}, {})...'.format(basis_1, basis_2, targ_ra, targ_dec))

fig.add_subplot(gs[0,0])
diagnostic_plots.densityPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd, 'stars')

fig.add_subplot(gs[1,0])
diagnostic_plots.densityPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd, 'galaxies')

fig.add_subplot(gs[2,0])
diagnostic_plots.densityPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd, 'blue_stars')

fig.add_subplot(gs[0,1])
diagnostic_plots.cmPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd, 'stars')

fig.add_subplot(gs[1,1])
diagnostic_plots.cmPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd, 'galaxies')

fig.add_subplot(gs[0,2])
diagnostic_plots.hessPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd)

fig.add_subplot(gs[1,2])
diagnostic_plots.starPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd)

fig.add_subplot(gs[2,1:3])
diagnostic_plots.radialPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd)

# Name
try: # ugali
    association_string = candidate_list[candidate]['NAME']
except: # simple
    # Check for possible associations
    glon_peak, glat_peak = ugali.utils.projector.celToGal(targ_ra, targ_dec)
    catalog_array = ['McConnachie15', 'Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'WEBDA14','ExtraDwarfs','ExtraClusters']
    catalog = ugali.candidate.associate.SourceCatalog()
    for catalog_name in catalog_array:
        catalog += ugali.candidate.associate.catalogFactory(catalog_name)

    idx1, idx2, sep = catalog.match(glon_peak, glat_peak, tol=0.5, nnearest=1)
    match = catalog[idx2]
    if len(match) > 0:
        association_string = '{} at {:0.3f} deg'.format(match[0]['name'], float(sep))
    else:
        association_string = 'No association within 0.5 deg'

# Currently in development -- not passing candidate into this program anymore
#try: # ugali
#    association_string = candidate_list[candidate]['NAME']
#except: # simple
#    assoc  = candidate_list[candidate]['ASSOC']
#    angsep = candidate_list[candidate]['ANGSEP']
#    if assoc == 'nan':
#        association_string = 'No association within 0.5 deg'
#    else:
#        association_string = '{} at {:0.3f} deg'.format(assoc, float(angsep))

association_string = str(np.char.strip(association_string))
#association_string = association_string.encode('utf-8')
association_string = repr(association_string)

plt.suptitle('{}\n'.format(association_string) + r'($\alpha$, $\delta$, $\mu$, $\sigma$) = ({}, {}, {}, {})'.format(targ_ra, targ_dec, mod, sig), fontsize=24)

file_name = 'candidate_{}_{}'.format(targ_ra, targ_dec)
plt.savefig(save_dir+'/'+file_name+'.png',  bbox_inches='tight')
plt.close()
