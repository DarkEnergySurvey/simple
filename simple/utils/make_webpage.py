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
#import www

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    simple_dir = cfg['setup']['simple_dir']
    candidate_list = cfg[survey]['candidate_list']
    save_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'])

outfile = 'test.html'

plots = os.listdir(save_dir)

#www.create_index_html(outfile, candidates)

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
      <th></th>
      <th></th>
      <th></th>
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
    </tr>  
"""

def create_index_html(filename,candidates):
    """Create the index.html"""
    tablerows = [ROW%dict(name=plot, fname=plot) for plot in plots if plot.endswith('.png')]

    table = TABLE%dict(rows='\n'.join(tablerows))
    index = INDEX%dict(table=table)
        
    with open(filename,'w') as out:
        out.write(index)

create_index_html(outfile, plots)
