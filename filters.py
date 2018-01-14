#!/usr/bin/env python
"""
Filters for different surveys
"""
__author__ = "Sidney Mau"

import yaml
import numpy as np
from matplotlib import mlab

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

mag_g = cfg[cfg['data']]['mag_g']
mag_r = cfg[cfg['data']]['mag_r']

def quality_filter(survey, data):
    """Return data above a quality threshold"""
    if survey == 'des':
        cut = (data['PSF_MAG_SFD_G'] < 25)
    elif survey == 'bliss':
        cut = (data['PSF_MAG_SFD_G'] < 25)
    elif survey == 'maglites':
        cut = (data['PSF_MAG_SFD_G'] < 25)
    elif survey == 'panstarrs':
        cut = (np.bitwise_and(data['QUALITYFLAG'], 16) > 0) \
            & (data['NSTACKDETECTIONS'] > 1) \
            & (np.bitwise_and(data['GINFOFLAG'], 8) == 0) \
            & (np.bitwise_and(data['RINFOFLAG'], 8) == 0) \
            & (np.bitwise_and(data['IINFOFLAG'], 8) == 0) \
            & (np.bitwise_and(data['GINFOFLAG'], 2048) == 0) \
            & (np.bitwise_and(data['RINFOFLAG'], 2048) == 0) \
            & (np.bitwise_and(data['IINFOFLAG'], 2048) == 0) \
            & (np.bitwise_and(data['GINFOFLAG2'], 4194304) == 0) \
            & (np.bitwise_and(data['RINFOFLAG2'], 4194304) == 0) \
            & (np.bitwise_and(data['IINFOFLAG2'], 4194304) == 0) \
            & (np.bitwise_and(data['GINFOFLAG2'], 8192) == 0) \
            & (np.bitwise_and(data['RINFOFLAG2'], 8192) == 0) \
            & (np.bitwise_and(data['IINFOFLAG2'], 8192) == 0)
    return cut

def star_filter(survey, data):
    """Return stellar-like objects"""
    if survey == 'des':
        cut = (data['EXTENDED_CLASS_MASH'] >= 0) \
            & (data['EXTENDED_CLASS_MASH'] <= 2)
    elif survey == 'bliss':
        cut = (data['CM_T'] < 0.003) # use CM_T_ERR?
    elif survey == 'maglites':
        cut = (data['CM_T'] < 0.003) # use CM_T_ERR?
    elif survey == 'panstarrs':
        cut = ((data['IFPSFMAG'] - data['IFKRONMAG']) < 0.05) 
    return cut

def galaxy_filter(survey, data):
    """Return stellar-like objects"""
    if survey == 'des':
        cut = (data['EXTENDED_CLASS_MASH'] > 2)
    elif survey == 'bliss':
        cut = (data['CM_T'] > 0.005) # use CM_T_ERR?
    elif survey == 'maglites':
        cut = (data['CM_T'] > 0.005) # use CM_T_ERR?
    elif survey == 'panstarrs':
        cut = ((data['IFPSFMAG'] - data['IFKRONMAG']) > 0.05) # just a guess
    return cut

def color_filter(survey, data):
    """Return blue objects"""
    if survey == 'des':
        cut = ((data['PSF_MAG_SFD_G'] - data['PSF_MAG_SFD_R']) < 0.4) # 0.2
    #elif survey == 'bliss':
    #elif survey == 'maglites':
    elif survey == 'panstarrs':
        cut = ((data['GFPSFMAG'] - data['IFPSFMAG']) > -0.5) \
            & ((data['GFPSFMAG'] - data['IFPSFMAG']) < 1.0)
    return cut

def dered_mag(survey, data):
    """Return the data with an added flag for dereddened (extinction
       corrected) magnitude"""
    if survey == 'des':
        data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['PSF_MAG_SFD_G'], data['PSF_MAG_SFD_R']])
    elif survey == 'bliss':
        #data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['CM_MAG_G'] - data['EXINCTION_G'], data['CM_MAG_R'] - data['EXTINCTION_R']])
        data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['PSF_MAG_SFD_G'], data['PSF_MAG_SFD_R']])
    elif survey == 'maglites':
        data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['WAVG_MAG_PSF_G'] - data['EXINCTION_G'], data['WAVG_MAG_PSF_R'] - data['EXTINCTION_R']])
    elif survey == 'panstarrs':
        data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['EXTSFD_G'], data['EXTSFD_R']]) # is this right?
    return data
