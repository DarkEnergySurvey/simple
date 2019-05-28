#!/usr/bin/env python
"""
Simple binning search algorithm
"""
__author__ = "Keith Bechtol"

# Python libraries
import sys
import os
import glob
import yaml
from matplotlib import mlab
import numpy as np
import healpy as hp
import astropy.io.fits as pyfits
import fitsio as fits

# Ugali libraries
import ugali.utils.healpix
import ugali.utils.projector

# Simple binner modules
import simple.filters
import simple.simple_utils

###########################################################

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    nside = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']
    mag_max = cfg[survey]['mag_max']

    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']

    band_1 = cfg[survey]['band_1']
    band_2 = cfg[survey]['band_2']
    mag = cfg[survey]['mag']
    mag_err = cfg[survey]['mag_err']
    mag_dered = cfg[survey]['mag_dered']

    mode = cfg[survey]['mode']
    sim_population = cfg[survey]['sim_population']
    sim_dir = cfg[survey]['sim_dir']
    #object_list = cfg[survey]['object_list']

    fracdet_map = cfg[survey]['fracdet']

    #mag_g = cfg[survey]['mag_g']
    #mag_r = cfg[survey]['mag_r']
    #mag_g_err = cfg[survey]['mag_g_err']
    #mag_r_err = cfg[survey]['mag_r_err']

    results_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

############################################################

# construct mags
mag_1 = mag.format(band_1.upper())
mag_2 = mag.format(band_2.upper())
mag_err_1 = mag_err.format(band_1.upper())
mag_err_2 = mag_err.format(band_2.upper())
mag_dered_1 = mag_dered.format(band_1.upper())
mag_dered_2 = mag_dered.format(band_2.upper())

# main

try:
    ra_select, dec_select = float(sys.argv[1]), float(sys.argv[2])
except:
    sys.exit('ERROR! Coordinates not given in correct format.')

try:
    mc_source_id = float(sys.argv[3])
except:
    mc_source_id = 0

try:
    outfile = sys.argv[4]
except:
    sys.exit('ERROR: no outfile given.')
#    if (mode == 0):
#        outfile = '{}/results_nside_{}_{}.csv'.format(results_dir, nside, pix_nside_select)
#    elif (mode == 1):
#        outfile = '{}/results_mc_source_id_{}.csv'.format(results_dir, mc_source_id_array[0]) # all values in mc_source_id_array should be the same
    

print('Search coordinates: (RA, Dec) = ({:0.2f}, {:0.2f})').format(ra_select, dec_select)

# Now cut for a single pixel
pix_nside_select = ugali.utils.healpix.angToPix(nside, ra_select, dec_select)
#ra_select, dec_select = ugali.utils.healpix.pixToAng(nside, pix_nside_select)
pix_nside_neighbors = np.concatenate([[pix_nside_select], hp.get_all_neighbours(nside, pix_nside_select)])
print('Center healpixel: {}'.format(pix_nside_select))
print('Healpixels: {}'.format(pix_nside_neighbors))

# Construct data
#data = simple_utils.construct_modal_data(mode, pix_nside_neighbors, mc_source_id)
data = simple.simple_utils.construct_real_data(pix_nside_neighbors)

print('MC_SOURCE_ID = {}'.format(mc_source_id))
if (mode == 0):
    print('mode = 0: running only on real data')
elif (mode == 1):
    print('mode = 1: running on real data and simulated data')
    
    # inject objects for simulated object of mc_source_id
    sim_data = simple.simple_utils.construct_sim_data(pix_nside_neighbors, mc_source_id)
    data = simple.simple_utils.inject_sim(data, sim_data, mc_source_id)
elif (mode == 2):
    print('mode = 2: running only on real data')
else:
    print('No/unsupported mode specified; running only on real data')

# Quality cut
quality = simple.filters.quality_filter(survey, data)
data = data[quality]

# Deredden magnitudes
data = simple.filters.dered_mag(survey, data)

print('Found {} objects...').format(len(data))
if (len(data) == 0):
    print('Ending search prematurely. Look at data for debugging.')
    nan_array = [np.nan]
    simple.simple_utils.write_output(results_dir, nside, pix_nside_select,
                             nan_array, nan_array, nan_array, nan_array, 
                             nan_array, nan_array, nan_array, nan_array,
                             [mc_source_id], mode, outfile)
    exit()

print('Applying cuts...')
cut = simple.filters.star_filter(survey, data)
cut_gal = simple.filters.galaxy_filter(survey, data)

data_gal = data[cut_gal] # this isn't used at all other than for noting number of galaxy-like objects in ROI
data = data[cut]

print('{} star-like objects in ROI...'.format(len(data)))
print('{} galaxy-like objects in ROI...'.format(len(data_gal)))
if (mode == 1):
    print('{} simulated objects in ROI...'.format(np.sum(data['MC_SOURCE_ID'] != 0)))

# Read in fracdet map
if (fracdet_map is not None) and (fracdet_map.lower().strip() != 'none') and (fracdet_map != ''):
    print('Reading fracdet map {} ...').format(fracdet_map)
    fracdet = ugali.utils.healpix.read_map(fracdet_map)
else:
    print('No fracdet map specified ...')
    fracdet = None

distance_modulus_search_array = np.arange(16., mag_max, 0.5)

ra_peak_array = []
dec_peak_array = [] 
r_peak_array = []
sig_peak_array = []
distance_modulus_array = []
mc_source_id_array = []
n_obs_peak_array = []
n_obs_half_peak_array = []
n_model_peak_array = []

if (mode == 0):
    for distance_modulus in distance_modulus_search_array:
        ra_peaks, dec_peaks, r_peaks, sig_peaks, dist_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = simple.simple_utils.search_by_distance(nside, data, distance_modulus, pix_nside_select, ra_select, dec_select, mag_max, fracdet)
        ra_peak_array.append(ra_peaks)
        dec_peak_array.append(dec_peaks)
        r_peak_array.append(r_peaks)
        sig_peak_array.append(sig_peaks)
        distance_modulus_array.append(dist_moduli)
        n_obs_peak_array.append(n_obs_peaks)
        n_obs_half_peak_array.append(n_obs_half_peaks)
        n_model_peak_array.append(n_model_peaks)
        mc_source_id_array.append(np.tile(0, len(sig_peaks)))
elif (mode == 1):
    # grab distance_modulus from population
    sim_pop = fits.read(sim_population)
    distance_modulus_select = sim_pop[sim_pop['MC_SOURCE_ID'] == mc_source_id]['DISTANCE_MODULUS'][0]

    distance_modulus = distance_modulus_search_array[np.argmin(np.fabs(distance_modulus_search_array - distance_modulus_select))]
    ra_peaks, dec_peaks, r_peaks, sig_peaks, dist_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = simple.simple_utils.search_by_simulation(nside, data, distance_modulus, pix_nside_select, ra_select, dec_select, mag_max, fracdet)
    ra_peak_array.append(ra_peaks)
    dec_peak_array.append(dec_peaks)
    r_peak_array.append(r_peaks)
    sig_peak_array.append(sig_peaks)
    distance_modulus_array.append(dist_moduli)
    n_obs_peak_array.append(n_obs_peaks)
    n_obs_half_peak_array.append(n_obs_half_peaks)
    n_model_peak_array.append(n_model_peaks)
    mc_source_id_array.append(np.tile(mc_source_id, len(sig_peaks)))
elif (mode == 2):
    # grab distance_modulus from population
    #sim_pop = fits.read(object_list)
    sim_pop = np.genfromtxt(object_list, delimiter=',', names=['RA', 'DEC', 'DISTANCE_MODULUS', 'MC_SOURCE_ID', 'NAME'])[1:]
    distance_modulus_select = sim_pop[sim_pop['MC_SOURCE_ID'] == mc_source_id]['DISTANCE_MODULUS'][0]

    distance_modulus = distance_modulus_search_array[np.argmin(np.fabs(distance_modulus_search_array - distance_modulus_select))]
    ra_peaks, dec_peaks, r_peaks, sig_peaks, dist_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = simple.simple_utils.search_by_simulation(nside, data, distance_modulus, pix_nside_select, ra_select, dec_select, mag_max, fracdet)
    ra_peak_array.append(ra_peaks)
    dec_peak_array.append(dec_peaks)
    r_peak_array.append(r_peaks)
    sig_peak_array.append(sig_peaks)
    distance_modulus_array.append(dist_moduli)
    n_obs_peak_array.append(n_obs_peaks)
    n_obs_half_peak_array.append(n_obs_half_peaks)
    n_model_peak_array.append(n_model_peaks)
    mc_source_id_array.append(np.tile(mc_source_id, len(sig_peaks)))

ra_peak_array = np.concatenate(ra_peak_array)
dec_peak_array = np.concatenate(dec_peak_array)
r_peak_array = np.concatenate(r_peak_array)
sig_peak_array = np.concatenate(sig_peak_array)
distance_modulus_array = np.concatenate(distance_modulus_array)
n_obs_peak_array = np.concatenate(n_obs_peak_array)
n_obs_half_peak_array = np.concatenate(n_obs_half_peak_array)
n_model_peak_array = np.concatenate(n_model_peak_array)
mc_source_id_array = np.concatenate(mc_source_id_array)

# Sort peaks according to significance
index_sort = np.argsort(sig_peak_array)[::-1]
ra_peak_array = ra_peak_array[index_sort]
dec_peak_array = dec_peak_array[index_sort]
r_peak_array = r_peak_array[index_sort]
sig_peak_array = sig_peak_array[index_sort]
distance_modulus_array = distance_modulus_array[index_sort]
n_obs_peak_array = n_obs_peak_array[index_sort]
n_obs_half_peak_array = n_obs_half_peak_array[index_sort]
n_model_peak_array = n_model_peak_array[index_sort]
mc_source_id_array = mc_source_id_array[index_sort]

# Collect overlapping peaks
for ii in range(0, len(sig_peak_array)):
    if sig_peak_array[ii] < 0:
        continue
    angsep = ugali.utils.projector.angsep(ra_peak_array[ii], dec_peak_array[ii], ra_peak_array, dec_peak_array)
    sig_peak_array[(angsep < r_peak_array[ii]) & (np.arange(len(sig_peak_array)) > ii)] = -1.
    #sig_peak_array[(angsep < 0.5) & (np.arange(len(sig_peak_array)) > ii)] = -1. # 0.5 deg radius

if (mode == 0):
    # Prune the list of peaks
    ra_peak_array = ra_peak_array[sig_peak_array > 0.]
    dec_peak_array = dec_peak_array[sig_peak_array > 0.]
    r_peak_array = r_peak_array[sig_peak_array > 0.]
    distance_modulus_array = distance_modulus_array[sig_peak_array > 0.]
    n_obs_peak_array = n_obs_peak_array[sig_peak_array > 0.]
    n_obs_half_peak_array = n_obs_half_peak_array[sig_peak_array > 0.]
    n_model_peak_array = n_model_peak_array[sig_peak_array > 0.]
    mc_source_id_array = mc_source_id_array[sig_peak_array > 0.]
    sig_peak_array = sig_peak_array[sig_peak_array > 0.] # Update the sig_peak_array last!

for ii in range(0, len(sig_peak_array)):
    print('{:0.2f} sigma; (RA, Dec, d) = ({:0.2f}, {:0.2f}); r = {:0.2f} deg; d = {:0.1f}, mu = {:0.2f} mag), mc_source_id: {:0.2f}'.format(sig_peak_array[ii], 
                 ra_peak_array[ii], 
                 dec_peak_array[ii], 
                 r_peak_array[ii],
                 ugali.utils.projector.distanceModulusToDistance(distance_modulus_array[ii]),
                 distance_modulus_array[ii],
                 mc_source_id_array[ii]))

# Write output
if (len(sig_peak_array) > 0):
    simple.simple_utils.write_output(results_dir, nside, pix_nside_select, ra_peak_array, dec_peak_array, r_peak_array, distance_modulus_array, 
                             n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, 
                             sig_peak_array, mc_source_id_array, mode, outfile)
else:
    print('No significant hotspots found.')
    nan_array = [np.nan]
    simple.simple_utils.write_output(results_dir, nside, pix_nside_select,
                             nan_array, nan_array, nan_array, nan_array, 
                             nan_array, nan_array, nan_array, nan_array,
                             [mc_source_id], mode, outfile)
