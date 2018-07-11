#!/usr/bin/env python
"""
Filters for different surveys
"""
__author__ = "Sidney Mau"

import yaml
import numpy as np
from matplotlib import mlab
import numpy.lib.recfunctions
#import ugali.utils.mlab

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']

    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']
    mag_g = cfg[survey]['mag_g']
    mag_r = cfg[survey]['mag_r']
    mag_max = cfg[survey]['mag_max']

    band_1 = cfg[survey]['band_1']
    band_2 = cfg[survey]['band_2']
    mag = cfg[survey]['mag']
    mag_err = cfg[survey]['mag_err']
    mag_dered = cfg[survey]['mag_dered']

# construct mags
mag_1 = mag + band_1
mag_2 = mag + band_2
mag_err_1 = mag_err + band_1
mag_err_2 = mag_err + band_2
mag_dered_1 = mag_dered + band_1
mag_dered_2 = mag_dered + band_2

###

def quality_filter(survey, data):
    """Return data above a quality threshold"""
    if survey == 'y3_gold_2_0':
        cut = (data[mag_1] < mag_max)
    elif survey == 'y3a2':
        cut = (data['PSF_MAG_SFD_G'] < mag_max)
    elif survey == 'bliss':
        #cut = (data['PSF_MAG_SFD_G'] < 25)
        #cut = (data[mag_g] < 25)
        cut = (data['WAVG_MAG_PSF_G'] < mag_max) \
            & (data['MAG_PSF_G'] < 90) \
            & (data['MAG_PSF_R'] < 90)
            #& (data['SEXTRACTOR_FLAGS_G'] < 4) \
            #& (data['SEXTRACTOR_FLAGS_R'] < 4)
            #& ((data['PSF_MAG_SFD_G'] - data['PSF_MAG_SFD_R']) < 1.)
    elif survey == 'maglites':
        cut = (data['PSF_MAG_SFD_G'] < mag_max) \
            & (data['SEXTRACTOR_FLAGS_G'] < 4) \
            & (data['SEXTRACTOR_FLAGS_R'] < 4)
            #& ((data['PSF_MAG_SFD_G'] - data['PSF_MAG_SFD_R']) < 1.)
    elif survey == 'panstarrs':
        #cut = (np.bitwise_and(data['QUALITYFLAG'], 16) > 0) \
        #    & (data['NSTACKDETECTIONS'] > 1) \
        #    & (data['NDETECTIONS'] > 0) \
        #    & (np.bitwise_and(data['GINFOFLAG'], 8) == 0) \
        #    & (np.bitwise_and(data['RINFOFLAG'], 8) == 0) \
        #    & (np.bitwise_and(data['IINFOFLAG'], 8) == 0) \
        #    & (np.bitwise_and(data['GINFOFLAG'], 2048) == 0) \
        #    & (np.bitwise_and(data['RINFOFLAG'], 2048) == 0) \
        #    & (np.bitwise_and(data['IINFOFLAG'], 2048) == 0) \
        #    & (np.bitwise_and(data['GINFOFLAG2'], 4194304) == 0) \
        #    & (np.bitwise_and(data['RINFOFLAG2'], 4194304) == 0) \
        #    & (np.bitwise_and(data['IINFOFLAG2'], 4194304) == 0) \
        #    & (np.bitwise_and(data['GINFOFLAG2'], 8192) == 0) \
        #    & (np.bitwise_and(data['RINFOFLAG2'], 8192) == 0) \
        #    & (np.bitwise_and(data['IINFOFLAG2'], 8192) == 0) \
        #    & (data['RFPSFMAGERR'] < 0.1) # replacing (data['GFPSFMAG'] < 22.5) after Keith's investigations
        #    #& (data['GFPSFMAG'] < 22.5) # observed magnitude - not extinction corrected
        #    #& (data['GINFOFLAG'] >= 0) # recommended by Alex; untested yet
        cut = (data['RFPSFMAGERR'] < 0.1) # replacing (data['GFPSFMAG'] < 22.5) after Keith's investigations
    return cut

def star_filter(survey, data):
    """Return stellar-like objects"""
    if survey == 'y3_gold_2_0':
        cut = (data['EXTENDED_CLASS_MASH_SOF'] >= 0) \
            & (data['EXTENDED_CLASS_MASH_SOF'] <= 2)
    elif survey == 'y3a2':
        cut = (data['EXTENDED_CLASS_MASH'] >= 0) \
            & (data['EXTENDED_CLASS_MASH'] <= 2)
    elif survey == 'bliss':
        #cut = (data['CM_T'] < 0.003 + data['CM_T_ERR'])
        cut = (data['WAVG_SPREAD_MODEL_G'] < 0.003 + data['SPREADERR_MODEL_G'])
    elif survey == 'maglites':
        cut = (data['CM_T'] < 0.003 + data['CM_T_ERR'])
    elif survey == 'panstarrs':
        #cut = ((data['IFPSFMAG'] - data['IFKRONMAG']) < 0.05)
        ## Label objects that fail fit stars (this happens in dense regions)
        #cut = cut | (data['IFPSFMAG'] == -999.) 
        ## Label objects that have exceptionally bright Kron magnitudes relative to PSF magnitudes stars (this happens in dense regions)
        #cut = cut | ((data['IFPSFMAG'] - data['IFKRONMAG']) > 4.0)
        cut = (data['EXTENDED_CLASS'] == 0)
    return cut

def galaxy_filter(survey, data):
    """Return stellar-like objects"""
    if survey == 'y3_gold_2_0':
        cut = (data['EXTENDED_CLASS_MASH_SOF'] > 2)
    elif survey == 'y3a2':
        cut = (data['EXTENDED_CLASS_MASH'] > 2)
    elif survey == 'bliss':
        #cut = (data['CM_T'] > 0.003 + data['CM_T_ERR']) # 0.005?
        cut = (data['WAVG_SPREAD_MODEL_G'] > 0.003 + data['SPREADERR_MODEL_G'])
    elif survey == 'maglites':
        cut = (data['CM_T'] > 0.003 + data['CM_T_ERR']) # 0.005?
    elif survey == 'panstarrs':
        #cut = ((data['IFPSFMAG'] - data['IFKRONMAG']) > 0.05) # just a guess
        cut = (data['EXTENDED_CLASS'] == 1)
    return cut

def color_filter(survey, data):
    """Return blue objects"""
    if survey == 'y3_gold_2_0':
        #cut = ((data['SOF_PSF_MAG_CORRECTED_G'] - data['SOF_PSF_MAG_CORRECTED_R']) < 0.4) # 0.2
        cut = ((data[mag_1] - data[mag_2]) < 0.4) # 0.2
    elif survey == 'y3a2':
        cut = ((data['PSF_MAG_SFD_G'] - data['PSF_MAG_SFD_R']) < 0.4) # 0.2
    elif survey == 'bliss':
        cut = ((data[mag_g] - data[mag_r]) < 0.4) # 0.2
    elif survey == 'maglites':
        cut = ((data[mag_g] - data[mag_r]) < 0.4) # 0.2
    elif survey == 'panstarrs':
        cut = ((data[mag_g] - data[mag_r]) < 0.4) # 0.2
        #cut = ((data['GFPSFMAG'] - data['IFPSFMAG']) > -0.5) \
        #    & ((data['GFPSFMAG'] - data['IFPSFMAG']) < 1.0)
    return cut

def dered_mag(survey, data):
    """Return the data with an added flag for dereddened (extinction
       corrected) magnitude"""
    if survey == 'y3_gold_2_0':
        #data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['SOF_PSF_MAG_CORRECTED_G'], data['SOF_PSF_MAG_CORRECTED_R']])
        data = numpy.lib.recfunctions.append_fields(data, [mag_dered_1, mag_dered_2], [data[mag_1], data[mag_2]], usemask=False, asrecarray=True)
        #data = ugali.utils.mlab.rec_append_fields(data, [mag_g, mag_r], [data['SOF_PSF_MAG_CORRECTED_G'], data['SOF_PSF_MAG_CORRECTED_R']])
    elif survey == 'y3a2':
        #data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['PSF_MAG_SFD_G'], data['PSF_MAG_SFD_R']])
        data = numpy.lib.recfunctions.append_fields(data, [mag_g, mag_r], [data['PSF_MAG_SFD_G'], data['PSF_MAG_SFD_R']], usemask=False, asrecarray=True)
        #data = ugali.utils.mlab.rec_append_fields(data, [mag_g, mag_r], [data['PSF_MAG_SFD_G'], data['PSF_MAG_SFD_R']])
    elif survey == 'bliss':
        #data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['CM_MAG_G'] - data['EXINCTION_G'], data['CM_MAG_R'] - data['EXTINCTION_R']])
        #data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['WAVG_MAG_PSF_G'], data['WAVG_MAG_PSF_R']])
        data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['MAG_PSF_SFD_G'], data['MAG_PSF_SFD_R']])
        #data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['PSF_MAG_SFD_G'], data['PSF_MAG_SFD_R']])
        #data = numpy.lib.recfunctions.append_fields(data, [mag_g, mag_r], [data['PSF_MAG_SFD_G'], data['PSF_MAG_SFD_R']], 
        #                                            usemask=False, asrecarray=True)
        #data = ugali.utils.mlab.rec_append_fields(data, [mag_g, mag_r], [data['PSF_MAG_SFD_G'], data['PSF_MAG_SFD_R']])
    elif survey == 'maglites':
        #data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['WAVG_MAG_PSF_G'] - data['EXINCTION_G'], data['WAVG_MAG_PSF_R'] - data['EXTINCTION_R']])
        data = numpy.lib.recfunctions.append_fields(data, [mag_g, mag_r], [data['WAVG_MAG_PSF_G'] - data['EXINCTION_G'], data['WAVG_MAG_PSF_R'] - data['EXTINCTION_R']], usemask=False, asrecarray=True)
        #data = ugali.uitls.mlab.rec_append_fields(data, [mag_g, mag_r], [data['WAVG_MAG_PSF_G'] - data['EXINCTION_G'], data['WAVG_MAG_PSF_R'] - data['EXTINCTION_R']])
    elif survey == 'panstarrs':
        #data = mlab.rec_append_fields(data, [mag_g, mag_r], [data['GFPSFMAG'] - data['EXTSFD_G'], data['RFPSFMAG'] - data['EXTSFD_R']])
        #data = numpy.lib.recfunctions.append_fields(data, [mag_g, mag_r], [data['GFPSFMAG'] - data['EXTSFD_G'], data['RFPSFMAG'] - data['EXTSFD_R']], 
        #                                            usemask=False, asrecarray=True)
        #data = ugali.utils.mlab.rec_append_fields(data, [mag_g, mag_r], [data['GFPSFMAG'] - data['EXTSFD_G'], data['RFPSFMAG'] - data['EXTSFD_R']])
        data = numpy.lib.recfunctions.append_fields(data, [mag_g, mag_r], [data['GFPSFMAG_SFD'], data['RFPSFMAG_SFD']], usemask=False, asrecarray=True)
    return data
