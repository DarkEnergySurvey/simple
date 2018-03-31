#!/usr/bin/env python
"""
Quick plot of hotspot distribution
"""
__author__ = "Sidney Mau"

# Set the backend first!
import matplotlib
matplotlib.use('Agg')

import fitsio as fits
import numpy as np
import pylab as plt
import glob
import yaml
import ugali
import ugali.candidate.associate

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']
    candidate_list = cfg[survey]['candidate_list']

data = fits.read(candidate_list)
#data = data[data['SIG'] > 15]

plt.figure()
plt.xlabel("RA")
plt.ylabel("Dec")
plt.title("{} Hotspots".format(survey))

# Plot hotspots
plt.scatter(data['RA'], data['DEC'], edgecolor='none', s=1, c='black', alpha=0.5)

# Overplot known objects
catalog_array = ['McConnachie15', 'ExtraDwarfs']
catalog = ugali.candidate.associate.SourceCatalog()
for catalog_name in catalog_array:
    catalog += ugali.candidate.associate.catalogFactory(catalog_name)
plt.scatter(catalog['ra'], catalog['dec'], edgecolor='red', s=20, linewidth=0.5, c='none', label='known objects')

plt.xlim([min(data['RA']), max(data['RA'])])
plt.ylim([min(data['DEC']), max(data['DEC'])])
plt.gca().invert_xaxis()
plt.legend(loc='upper right')
plt.savefig("{}_hotspots.png".format(survey), bbox_inches='tight', dpi=300)
plt.close()
