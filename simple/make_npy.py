#!/usr/bin/env python
"""
Compile candidate list from results_dir
"""
__author__ = "Sidney Mau"

import glob
import yaml

import numpy as np

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']

    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']
    candidate_list = cfg[survey]['candidate_list']


fs = glob.glob('{}/*.npy'.format(cfg['output']['results_dir']))
data = np.concatenate([np.load(f) for f in fs])
np.save(candidate_list, data)

# Diagnostic output
data = np.load(candidate_list)
print("{} hotspots found.").format(len(data))
cut_0 = (data['SIG'] > 5.5)
print("{} hotspots found with SIG > 5.5.").format(len(data[cut_0]))
cut_1 = (data['SIG'] > 10)
print("{} hotspots found with SIG > 10.").format(len(data[cut_1]))
cut_2 = (data['SIG'] > 15)
print("{} hotspots found with SIG > 15.").format(len(data[cut_2]))
cut_3 = (data['SIG'] > 20)
print("{} hotspots found with SIG > 20.").format(len(data[cut_3]))
cut_4 = (data['SIG'] > 25)
print("{} hotspots found with SIG > 25.").format(len(data[cut_4]))
cut_5 = (data['SIG'] >= 37.5)
print("{} hotspots found with SIG >= 37.55").format(len(data[cut_5]))
