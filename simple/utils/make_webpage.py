#!/usr/bin/env python
"""
Generating html
"""
__author__ = "Sidney Mau"

import glob
import numpy as np
import pandas as pd
import yaml
import fitsio as fits
import os
import ugali.utils.projector
import ugali.candidate.associate
import simple.utils.filter_candidates
import simple.utils.query_image

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    simple_dir = cfg['setup']['simple_dir']
    candidate_list = cfg[survey]['candidate_list']
    jobs = cfg['batch']['jobs']
    save_dir = cfg['output']['save_dir']
    log_dir = cfg['output']['log_dir']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']
    candidate_list = cfg[survey]['candidate_list']

log_dir = os.path.join(cfg['output']['save_dir'], cfg['output']['log_dir'])
if not os.path.exists(save_dir):
    os.mkdir(save_dir)
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

outfile = 'index.html'
if survey == 'panstarrs':
    survey = 'ps1'
elif survey == 'y3_gold_2_0':
    survey = 'des'

plots = os.listdir(save_dir)

################################################################################

#candidate_list = fits.read(candidate_list)
#candidate_list = candidate_list[(candidate_list['SIG'] < 7) & (candidate_list['SIG'] > 6)]
#candidate_list[::-1].sort(order='SIG')

#candidate_list = simple.utils.filter_candidates.make_cut(candidate_list)
#candidate_list[::-1].sort(order='SIG')

#candidate_list = np.genfromtxt('remains_des_simple.txt', names=True)
#candidate_list[::-1].sort(order='SIG')

candidate_list = fits.read(candidate_list)
candidate_list[::-1].sort(order='SIG')

################################################################################


INDEX = """
<html>
  <head>
    <title>Candidates</title>
</head>
<body>
%(table)s
</body>
</html>
"""

TABLE = """
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: center;">
      <th>Candidate</th>
      <th>Image</th>
      <th>Diagnostic Plots</th>
    </tr>
  </thead>
  <tbody>
%(rows)s
  </tbody>
"""

ROW = """
    <tr>
      <th>%(name)s</th>
      <td> <a id="%(fname)s"></a><a href="%(fname)s.png"><img src="%(fname)s.png" width="400"></a></td>
      <td> <a id="%(fname)s"></a><a href="%(fname)s.png"><img src="%(fname)s.png"></a></td>
    </tr>  
"""

def find_name(candidate):
    # Name
    try: # ugali
        association_string = candidate['NAME']
    except: # simple
        # Check for possible associations
        glon_peak, glat_peak = ugali.utils.projector.celToGal(candidate['RA'], candidate['DEC'])
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
    
    association_string = str(np.char.strip(association_string))
    association_string = repr(association_string)

    return(association_string)


def create_entry(candidate):
    """Create a diagnostic row for candidate"""

    print('Creating entry for {}'.format(candidate))
    name    = find_name(candidate)
    sig     = round(candidate['SIG'], 2)
    ra      = round(candidate['RA'], 2)
    #ra      = round(candidate['ra'], 2)
    dec     = round(candidate['DEC'], 2)
    #dec     = round(candidate['dec'], 2)
    #lon     =
    #lat     =
    mod     = round(candidate['MODULUS'], 2)
    #mod     = round(candidate['modulus'], 2)
    #mc_source_id = round(candidate['MC_SOURCE_ID'], 2)
    mc_source_id = 0
    #dist    =
    #use astropy to convert
    #how to get distance?
    
    candidate_file = 'candidate_{:0.2f}_{:0.2f}.png'.format(ra, dec)
    plot = '{}/{}'.format(save_dir, candidate_file)
    fname = plot.strip('.png')
    if candidate_file not in plots:
        print('Plot not found; making plot for {}'.format(name))
        if 'N_MODEL' in candidate_list.dtype.names:
            field_density = round(candidate['N_MODEL'] / (np.pi * (candidate['r'] * 60.)**2), 4) # field density (arcmin^-2)
        logfile = '{}/candidate_{:0.2f}_{:0.2f}.log'.format(log_dir, ra, dec)
        batch = 'csub -n {} -o {} --host all '.format(jobs, logfile) # testing condor updates
        if 'N_MODEL' in candidate_list.dtype.names:
            command = 'python {}/plotting/web_plot.py {:0.2f} {:0.2f} {:0.2f} {:0.2f} {:0.2f} {:0.4f}'.format(simple_dir, ra, dec, mod, sig, mc_source_id, field_density)
        else:
            command = 'python {}/plotting/web_plot.py {:0.2f} {:0.2f} {:0.2f} {:0.2f} {:0.2f}'.format(simple_dir, ra, dec, mod, sig, mc_source_id)
        command_queue = batch + command

        print(command_queue)
        os.system(command_queue) # Submit to queue

    image_file = 'image_{:0.2f}_{:0.2f}.png'.format(ra, dec)
    #if image_file not in plots:
    #    print('Image not found; retrieving image for {}'.format(name))
    #    image_url = simple.utils.query_image.retrieve_image(image_file, ra, dec, survey)
    image_url = simple.utils.query_image.retrieve_image(image_file, ra, dec, survey)
    image = '{}/{}'.format(save_dir, image_file)
    image_name = image.strip('.png')

    tablerow = ROW%dict(name=plot.strip('.png'), fname=save_dir+'/'+plot.strip('.png'))
    tablerow = """
        <tr>
          <th>{}<br>sig = {}<br>(RA, Dec) = ({}, {})<br>mod = {}</th>
          <td> <a id={}></a><a href={}><img src={} width="200"></a></td>
          <td> <a id={}></a><a href={}><img src={} width="1000"></a></td>
        </tr>  
    """.format(name, sig, ra, dec, mod, image_name, image_url, image, fname, plot, plot)

    return(tablerow)

def create_index_html(filename,candidate_list):
    """Create the index.html"""
    #tablerows = [ROW%dict(name=plot.strip('.png'), fname=save_dir+'/'+plot.strip('.png')) for plot in plots if plot.endswith('.png')]

    entries = [create_entry(candidate) for candidate in candidate_list]
    table = TABLE%dict(rows='\n'.join(entries))
    index = INDEX%dict(table=table)
        
    with open(filename,'w') as out:
        out.write(index)

################################################################################

create_index_html(outfile, candidate_list)
