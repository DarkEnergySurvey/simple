#!/usr/bin/env python
"""
Generating html
"""
__author__ = "Sidney Mau"

import yaml
import fits
import www

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    simple_dir = cfg['setup']['simple_dir']
    candidate_list = cfg[survey]['candidate_list']
    save_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'])
    log_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'], cfg['output']['log_dir'])

candidate_list = fits.read(candidate_list)

outfile = 'test.html'

www.create_index_html(outfile, candidates)
