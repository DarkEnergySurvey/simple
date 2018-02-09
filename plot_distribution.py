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

    survey = cfg['data']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']
    candidate_list = cfg[survey]['candidate_list']

data = fits.read(candidate_list)

ra = data['RA']
dec = data['DEC']

plt.figure()
plt.scatter(ra, dec, edgecolor='none', s=1, c='black', alpha=0.5)
plt.xlabel("RA")
plt.ylabel("Dec")
plt.title("{} Hotspots".format(survey))

# Overplot known dwarves
catalog_array = ['McConnachie15', 'ExtraDwarfs']
catalog = ugali.candidate.associate.SourceCatalog()
for catalog_name in catalog_array:
    catalog += ugali.candidate.associate.catalogFactory(catalog_name)
cut = ((catalog['ra'] >= min(ra)) & (catalog['ra'] <= max(ra)) & (catalog['dec'] >= min(dec)) & (catalog['dec'] <= max(dec)))
catalog = catalog[cut]
known_ra = catalog['ra']
known_dec = catalog['dec']
plt.scatter(known_ra, known_dec, edgecolor='none', s=5, c='red', label='known dwarfs')

plt.legend(loc='upper right')
plt.savefig("{}_hotspots.png".format(survey), bbox_inches='tight', dpi=300)
plt.close()
