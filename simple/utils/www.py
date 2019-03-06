#!/usr/bin/env python
"""
Generating html
"""
__author__ = "Alex Drlica-Wagner"
import glob

import os
import glob
from astropy.io import fits 
import numpy as np
import pandas as pd

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
      <td> <a id="%(fname)s_dist"></a><a href="%(fname)s_dist.png"><img src="%(fname)s_dist.png" width="400"></a></td>
      <td> <a id="%(fname)s_scat"></a><a href="%(fname)s_scat.png"><img src="%(fname)s_scat.png" width="400"></a></td>
    </tr>  
"""

def create_index_html(filename,candidates):
    """Create the index.html"""
    tablerows = []
    for idx,c in enumerate(candidates):
        tablerows.append(ROW%dict(name=c['name'],fname=c['name'].lower().replace(' ','_')))

    table = TABLE%dict(rows='\n'.join(tablerows))
    index = INDEX%dict(table=table)
        
    with open(filename,'w') as out:
        out.write(index)
