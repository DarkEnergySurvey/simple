#!/usr/bin/env python
"""
Filters for different surveys
"""
__author__ = "Sidney Mau"

#import yaml
import numpy as np

#with open('config.yaml', 'r') as ymlfile:
#    cfg = yaml.load(ymlfile)


## De-redden magnitudes
#try:
#	data = mlab.rec_append_fields(data, ['WAVG_MAG_PSF_DRED_G', 'WAVG_MAG_PSF_DRED_R'], [data[mag_g_dred_flag], data[mag_r_dred_flag]])
##except:
##    data = mlab.rec_append_fields(data, ['WAVG_MAG_PSF_DRED_G', 'WAVG_MAG_PSF_DRED_R'], [data[mag_g_flag] - data[extinction_g_flag], data[mag_r_flag] - data[extinction_r_flag]])
#except:
#	data = mlab.rec_append_fields(data, ['WAVG_MAG_PSF_DRED_G', 'WAVG_MAG_PSF_DRED_R'], [data[mag_g_flag], data[mag_r_flag]])
#
#mag_g = data['WAVG_MAG_PSF_DRED_G']
#mag_r = data['WAVG_MAG_PSF_DRED_R']


def quality_filter(survey, data):
    """Return data above a quality threshold"""
    if survey == 'des':
        cut = (data['PSF_MAG_SFD_G'] < 24)
    elif survey == 'bliss':
        cut = (data['PSF_MAG_SFD_G'] < 24)
    elif survey == 'maglites':
        cut = (data['PSF_MAG_SFD_G'] < 24)
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
    data = data[cut]

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
    data = data[cut]

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
    data = data[cut]

def color_filter(survey, data):
    """Return blue objects"""
    if survey == 'des':
        cut = ((data['PSF_MAG_SFD_G'] - data['PSF_MAG_SFD_R']) < 0.4) # 0.2
    elif survey == 'bliss':
    elif survey == 'maglites':
    elif survey == 'panstarrs':
        cut = ((data['GFPSFMAG'] - data['IFPSFMAG']) > -0.5) \
            & ((data['GFPSFMAG'] - data['IFPSFMAG']) < 1.0)
    data = data[cut]
