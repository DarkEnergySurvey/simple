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

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['data']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']
    candidate_list = cfg[survey]['candidate_list']

file = "candidate_list.fits"
data = fits.read(file)

ra = data['RA']
dec = data['DEC']

plt.figure()
plt.scatter(ra, dec, edgecolor='none', s=1, c='black')
plt.xlabel("RA")
plt.ylabel("Dec")
plt.title("{} Hotspots".format(survey))

plt.savefig("{}_hotspots.png".format(survey), bbox_inches='tight')
plt.close()
