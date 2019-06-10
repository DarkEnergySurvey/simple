#!/usr/bin/env python
"""
Simple binning search algorithm utilities
"""
__author__ = "Keith Bechtol, Sid Mau"

# Python libraries
import os
import glob
import yaml
import numpy as np
#from matplotlib import mlab
import numpy.lib.recfunctions
import healpy as hp
import astropy.io.fits as pyfits # migrate to fitsio
import fitsio as fits
import scipy.interpolate
import scipy.ndimage
import itertools

# Ugali libraries
#import ugali.utils.mlab
import ugali.utils.healpix
import ugali.utils.projector
import ugali.isochrone

########################################################################

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']
    isoname = cfg[survey]['isoname']
    isosurvey = cfg[survey]['isosurvey']
    mag_max = cfg[survey]['mag_max']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']

    mode = cfg[survey]['mode']
    sim_population = cfg[survey]['sim_population']
    sim_dir = cfg[survey]['sim_dir']

    fracdet_map = cfg[survey]['fracdet']
    
    #mag_g = cfg[survey]['mag_g']
    #mag_r = cfg[survey]['mag_r']
    #mag_g_err = cfg[survey]['mag_g_err']
    #mag_r_err = cfg[survey]['mag_r_err']

    band_1 = cfg[survey]['band_1']
    band_2 = cfg[survey]['band_2']
    mag = cfg[survey]['mag']
    mag_err = cfg[survey]['mag_err']
    mag_dered = cfg[survey]['mag_dered']

    results_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

# construct mags
mag_1 = mag.format(band_1.upper())
mag_2 = mag.format(band_2.upper())
mag_err_1 = mag_err.format(band_1.upper())
mag_err_2 = mag_err.format(band_2.upper())
mag_dered_1 = mag_dered.format(band_1.upper())
mag_dered_2 = mag_dered.format(band_2.upper())

# specify fits columns
#COLS = [basis_1, basis_2, mag_1, mag_2, mag_err_1, mag_err_2, mag_dered_1, mag_dered_2,...]

########################################################################

def get_catalog_file(catalog_dir, mc_source_id):
    """
    Inputs:
        catalog_dir = string corresponding to directory containing the stellar catalog infiles
        mc_source_id = integer corresponding the target MC_SOURCE_ID value
    Outputs:
        catalog_infile = string corresponding to filename of stellar catalog containing mc_source_id
    """
    catalog_infiles = sorted(glob.glob(catalog_dir + '/*catalog*.fits'))
    mc_source_id_array = []
    catalog_infile_index_array = []
    for ii, catalog_infile in enumerate(catalog_infiles):
        mc_source_id_min = int(os.path.basename(catalog_infile).split('.')[0].split('mc_source_id_')[-1].split('-')[0])
        mc_source_id_max = int(os.path.basename(catalog_infile).split('.')[0].split('mc_source_id_')[-1].split('-')[1])
        assert (mc_source_id_max > mc_source_id_min) & (mc_source_id_min >= 1), 'Found invalue MC_SOURCE_ID values in filenames'
        mc_source_id_array.append(np.arange(mc_source_id_min, mc_source_id_max + 1))
        catalog_infile_index_array.append(np.tile(ii, 1 + (mc_source_id_max - mc_source_id_min)))

    mc_source_id_array = np.concatenate(mc_source_id_array)
    catalog_infile_index_array = np.concatenate(catalog_infile_index_array)

    assert len(mc_source_id_array) == len(np.unique(mc_source_id_array)), 'Found non-unique MC_SOURCE_ID values in filenames'
    assert np.in1d(mc_source_id, mc_source_id_array), 'Requested MC_SOURCE_ID value not among files'
    mc_source_id_index = np.nonzero(mc_source_id == mc_source_id_array)[0][0] # second [0] added by smau 7/23/18 to fix incompatiable type bug
    return catalog_infiles[catalog_infile_index_array[mc_source_id_index]]

def construct_real_data(pix_nside_neighbors):
    data_array = []
    for pix_nside in pix_nside_neighbors:
        inlist = glob.glob('{}/*_{:05d}.fits'.format(datadir, pix_nside))
        for infile in inlist:
            if not os.path.exists(infile):
                continue
            #reader = pyfits.open(infile)
            #data_array.append(reader[1].data)
            #reader.close()
            data_array.append(fits.read(infile)) # add COLS=??? to speed up IO; use yaml
    data = np.concatenate(data_array)

    # Guarantee data has MC_SOURCE_ID
    try:
        data['MC_SOURCE_ID']
    except:
        # real data given MC_SOURCE_ID of 0
        #data = mlab.rec_append_fields(data, ['MC_SOURCE_ID'], [0])
        data = numpy.lib.recfunctions.append_fields(data, 'MC_SOURCE_ID', np.tile(0, len(data)), usemask=False, asrecarray=True)
        #data = ugali.utils.mlab.rec_append_field(data, ['MC_SOURCE_ID'], [0])
    
    return data

def construct_sim_data(pix_nside_neighbors, mc_source_id):
    data_array = []
    cat_file = get_catalog_file(sim_dir, mc_source_id)
    cat_data = fits.read(cat_file)
    pix = hp.ang2pix(nside,cat_data[basis_1],cat_data[basis_2],lonlat=True)
    pix_data = cat_data[np.in1d(pix, pix_nside_neighbors)]
    data_array.append(pix_data)
    data = np.concatenate(data_array)

    return data

def inject_sim(real_data, sim_data, mc_source_id):
    inject_data = sim_data[sim_data['MC_SOURCE_ID'] == mc_source_id]
    names = [name for name in real_data.dtype.names]
    inject_data = inject_data[names]
    # inject_data has boolean extended class
    # real_data has pseudo-boolean extended class
    real_data = real_data.astype(inject_data.dtype)
    data = np.concatenate((real_data, inject_data), axis=0)

    return data

def construct_modal_data(mode, pix_nside_neighbors, mc_source_id):
    print('Assembling data...')
    if (mode == 0): # real data
        data = construct_real_data(pix_nside_neighbors)
    elif (mode == 1): # real and simulated data
        data_real = construct_real_data(pix_nside_neighbors)
        data_sim  = construct_sim_data(pix_nside_neighbors, mc_source_id)

        #columns = data_sim.dtype
        #data_real = np.ndarray(data_real.shape, columns, data_real, 0, data_real.strides)
        #data = np.concatenate((data_real, data_sim), axis=0)
        data = inject_sim(data_real, data_sim, mc_source_id)
    else: # assume mode = 0
        data = construct_real_data(pix_nside_neighbors)

    return data

########################################################################

def cut_isochrone_path(g, r, g_err, r_err, isochrone, radius=0.1, return_all=False):
    """
    Cut to identify objects within isochrone cookie-cutter.
    """
    if np.all(isochrone.stage == 'Main'):
        # Dotter case
        index_transition = len(isochrone.stage)
    else:
        # Other cases
        index_transition = np.nonzero(isochrone.stage >= isochrone.hb_stage)[0][0] + 1    

    mag_1_rgb = isochrone.mag_1[0: index_transition] + isochrone.distance_modulus
    mag_2_rgb = isochrone.mag_2[0: index_transition] + isochrone.distance_modulus
    
    mag_1_rgb = mag_1_rgb[::-1]
    mag_2_rgb = mag_2_rgb[::-1]

    # Cut one way...
    f_isochrone = scipy.interpolate.interp1d(mag_2_rgb, mag_1_rgb - mag_2_rgb, bounds_error=False, fill_value = 999.)
    color_diff = np.fabs((g - r) - f_isochrone(r))
    cut_2 = (color_diff < np.sqrt(0.1**2 + r_err**2 + g_err**2))

     # ...and now the other
    f_isochrone = scipy.interpolate.interp1d(mag_1_rgb, mag_1_rgb - mag_2_rgb, bounds_error=False, fill_value = 999.)
    color_diff = np.fabs((g - r) - f_isochrone(g))
    cut_1 = (color_diff < np.sqrt(0.1**2 + r_err**2 + g_err**2))

    cut = np.logical_or(cut_1, cut_2)

    #mag_bins = np.arange(17., 24.1, 0.1)
    mag_bins = np.arange(17., mag_max+0.1, 0.1)
    mag_centers = 0.5 * (mag_bins[1:] + mag_bins[0:-1])
    magerr = np.tile(0., len(mag_centers))
    for ii in range(0, len(mag_bins) - 1):
        cut_mag_bin = (g > mag_bins[ii]) & (g < mag_bins[ii + 1])
        magerr[ii] = np.median(np.sqrt(0.1**2 + r_err[cut_mag_bin]**2 + g_err[cut_mag_bin]**2))

    if return_all:
        return cut, mag_centers[f_isochrone(mag_centers) < 100], (f_isochrone(mag_centers) + magerr)[f_isochrone(mag_centers) < 100], (f_isochrone(mag_centers) - magerr)[f_isochrone(mag_centers) < 100]
    else:
        return cut

########################################################################

def compute_char_density(nside, data, ra_select, dec_select, magnitude_threshold=mag_max, fracdet=None):
    """
    Compute the characteristic density of a region
    Convlve the field and find overdensity peaks
    """

    cut_magnitude_threshold = (data[mag_1] < magnitude_threshold)

    proj = ugali.utils.projector.Projector(ra_select, dec_select)
    x, y = proj.sphereToImage(data[basis_1][cut_magnitude_threshold], data[basis_2][cut_magnitude_threshold]) # Trimmed magnitude range for hotspot finding
    #x_full, y_full = proj.sphereToImage(data[basis_1], data[basis_2]) # If we want to use full magnitude range for significance evaluation
    delta_x = 0.01
    area = delta_x**2
    smoothing = 2. / 60. # Was 3 arcmin
    bins = np.arange(-8., 8. + 1.e-10, delta_x)
    centers = 0.5 * (bins[0: -1] + bins[1:])
    yy, xx = np.meshgrid(centers, centers)

    h = np.histogram2d(x, y, bins=[bins, bins])[0]

    h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)

    #cut_goodcoverage = (data['NEPOCHS_G'][cut_magnitude_threshold] >= 2) & (data['NEPOCHS_R'][cut_magnitude_threshold] >= 2)
    # expect NEPOCHS to be good in DES data

    delta_x_coverage = 0.1
    area_coverage = (delta_x_coverage)**2
    bins_coverage = np.arange(-5., 5. + 1.e-10, delta_x_coverage)
    h_coverage = np.histogram2d(x, y, bins=[bins_coverage, bins_coverage])[0]
    #h_goodcoverage = np.histogram2d(x[cut_goodcoverage], y[cut_goodcoverage], bins=[bins_coverage, bins_coverage])[0]
    h_goodcoverage = np.histogram2d(x, y, bins=[bins_coverage, bins_coverage])[0]

    n_goodcoverage = h_coverage[h_goodcoverage > 0].flatten()

    #characteristic_density = np.mean(n_goodcoverage) / area_coverage # per square degree
    characteristic_density = np.median(n_goodcoverage) / area_coverage # per square degree
    print('Characteristic density = {:0.1f} deg^-2').format(characteristic_density)

    # Use pixels with fracdet ~1.0 to estimate the characteristic density
    if fracdet is not None:
        fracdet_zero = np.tile(0., len(fracdet))
        cut = (fracdet != hp.UNSEEN)
        fracdet_zero[cut] = fracdet[cut]

        nside_fracdet = hp.npix2nside(len(fracdet))
        
        subpix_region_array = []
        for pix in np.unique(ugali.utils.healpix.angToPix(nside, data[basis_1], data[basis_2])):
            subpix_region_array.append(ugali.utils.healpix.subpixel(pix, nside, nside_fracdet))
        subpix_region_array = np.concatenate(subpix_region_array)

        # Compute mean fracdet in the region so that this is available as a correction factor
        cut = (fracdet[subpix_region_array] != hp.UNSEEN)
        mean_fracdet = np.mean(fracdet[subpix_region_array[cut]])

        # smau: this doesn't seem to be used in the non-local density estimation
        subpix_region_array = subpix_region_array[fracdet[subpix_region_array] > 0.99]
        subpix = ugali.utils.healpix.angToPix(nside_fracdet, 
                                              data[basis_1][cut_magnitude_threshold], 
                                              data[basis_2][cut_magnitude_threshold]) # Remember to apply mag threshold to objects
        characteristic_density_fracdet = float(np.sum(np.in1d(subpix, subpix_region_array))) \
                                         / (hp.nside2pixarea(nside_fracdet, degrees=True) * len(subpix_region_array)) # deg^-2
        print('Characteristic density fracdet = {:0.1f} deg^-2').format(characteristic_density_fracdet)
        
        # Correct the characteristic density by the mean fracdet value
        characteristic_density_raw = 1. * characteristic_density
        characteristic_density /= mean_fracdet 
        print('Characteristic density (fracdet corrected) = {:0.1f} deg^-2').format(characteristic_density)

    return characteristic_density

def compute_local_char_density(nside, data, characteristic_density, ra_select, dec_select, x_peak, y_peak, angsep_peak, magnitude_threshold=mag_max, fracdet=None):
    """
    Compute the local characteristic density of a region
    """

    # The following is all computed elsewhere but needed in here... should either
    # abstract into its own function or somehow else circumvent the need to copy
    magnitude_threshold = mag_max # make this function argument?
    cut_magnitude_threshold = (data[mag_dered_1] < magnitude_threshold)

    proj = ugali.utils.projector.Projector(ra_select, dec_select)
    x, y = proj.sphereToImage(data[basis_1][cut_magnitude_threshold], data[basis_2][cut_magnitude_threshold]) # Trimmed magnitude range for hotspot finding
    #x_full, y_full = proj.sphereToImage(data[basis_1], data[basis_2]) # If we want to use full magnitude range for significance evaluation

    # If fracdet map is available, use that information to either compute local density,
    # or in regions of spotty coverage, use the typical density of the region
    if fracdet is not None:
        # The following is copied from how it's used in compute_char_density
        fracdet_zero = np.tile(0., len(fracdet))
        cut = (fracdet != hp.UNSEEN)
        fracdet_zero[cut] = fracdet[cut]

        nside_fracdet = hp.npix2nside(len(fracdet))
        
        subpix_region_array = []
        for pix in np.unique(ugali.utils.healpix.angToPix(nside, data[basis_1], data[basis_2])):
            subpix_region_array.append(ugali.utils.healpix.subpixel(pix, nside, nside_fracdet))
        subpix_region_array = np.concatenate(subpix_region_array)

        # Compute mean fracdet in the region so that this is available as a correction factor
        cut = (fracdet[subpix_region_array] != hp.UNSEEN)
        mean_fracdet = np.mean(fracdet[subpix_region_array[cut]])

        subpix_region_array = subpix_region_array[fracdet[subpix_region_array] > 0.99]
        subpix = ugali.utils.healpix.angToPix(nside_fracdet, 
                                              data[basis_1][cut_magnitude_threshold], 
                                              data[basis_2][cut_magnitude_threshold]) # Remember to apply mag threshold to objects

        # This is where the local computation begins
        ra_peak, dec_peak = proj.imageToSphere(x_peak, y_peak)
        subpix_all = ugali.utils.healpix.angToDisc(nside_fracdet, ra_peak, dec_peak, 0.5)
        subpix_inner = ugali.utils.healpix.angToDisc(nside_fracdet, ra_peak, dec_peak, 0.3)
        subpix_annulus = subpix_all[~np.in1d(subpix_all, subpix_inner)]
        mean_fracdet = np.mean(fracdet_zero[subpix_annulus])
        print('mean_fracdet {}'.format(mean_fracdet))
        if mean_fracdet < 0.5:
            characteristic_density_local = characteristic_density
            print('characteristic_density_local baseline {}').format(characteristic_density_local)
        else:
            # Check pixels in annulus with complete coverage
            subpix_annulus_region = np.intersect1d(subpix_region_array, subpix_annulus)
            print('{} percent pixels with complete coverage'.format(float(len(subpix_annulus_region)) / len(subpix_annulus)))
            if (float(len(subpix_annulus_region)) / len(subpix_annulus)) < 0.25:
                characteristic_density_local = characteristic_density
                print('characteristic_density_local spotty {}'.format(characteristic_density_local))
            else:
                characteristic_density_local = float(np.sum(np.in1d(subpix, subpix_annulus_region))) \
                                               / (hp.nside2pixarea(nside_fracdet, degrees=True) * len(subpix_annulus_region)) # deg^-2
                print('characteristic_density_local cleaned up {}'.format(characteristic_density_local))
    else:
        # Compute the local characteristic density
        area_field = np.pi * (0.5**2 - 0.3**2)
        n_field = np.sum((angsep_peak > 0.3) & (angsep_peak < 0.5))
        characteristic_density_local = n_field / area_field

        # If not good azimuthal coverage, revert
        cut_annulus = (angsep_peak > 0.3) & (angsep_peak < 0.5) 
        #phi = np.degrees(np.arctan2(y_full[cut_annulus] - y_peak, x_full[cut_annulus] - x_peak)) # Use full magnitude range, NOT TESTED!!!
        phi = np.degrees(np.arctan2(y[cut_annulus] - y_peak, x[cut_annulus] - x_peak)) # Impose magnitude threshold
        h = np.histogram(phi, bins=np.linspace(-180., 180., 13))[0]
        if np.sum(h > 0) < 10 or np.sum(h > 0.5 * np.median(h)) < 10:
            #angsep_peak = np.sqrt((x - x_peak)**2 + (y - y_peak)**2)
            characteristic_density_local = characteristic_density

    print('Characteristic density local = {:0.1f} deg^-2 = {:0.3f} arcmin^-2'.format(characteristic_density_local, characteristic_density_local / 60.**2))

    return characteristic_density_local

########################################################################

def find_peaks(nside, data, characteristic_density, distance_modulus, pix_nside_select, ra_select, dec_select, magnitude_threshold=mag_max, fracdet=None):
    """
    Convolve field to find characteristic density and peaks within the selected pixel
    """

    # convolve field and find peaks
    cut_magnitude_threshold = (data[mag_dered_1] < magnitude_threshold)

    proj = ugali.utils.projector.Projector(ra_select, dec_select)
    x, y = proj.sphereToImage(data[basis_1][cut_magnitude_threshold], data[basis_2][cut_magnitude_threshold]) # Trimmed magnitude range for hotspot finding
    #x_full, y_full = proj.sphereToImage(data[basis_1], data[basis_2]) # If we want to use full magnitude range for significance evaluation
    delta_x = 0.01
    area = delta_x**2
    smoothing = 2. / 60. # Was 3 arcmin
    bins = np.arange(-8., 8. + 1.e-10, delta_x)
    centers = 0.5 * (bins[0: -1] + bins[1:])
    yy, xx = np.meshgrid(centers, centers)

    h = np.histogram2d(x, y, bins=[bins, bins])[0]
    
    h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)

    factor_array = np.arange(1., 5., 0.05)
    rara, decdec = proj.imageToSphere(xx.flatten(), yy.flatten())
    cutcut = (ugali.utils.healpix.angToPix(nside, rara, decdec) == pix_nside_select).reshape(xx.shape)
    threshold_density = 5 * characteristic_density * area
    for factor in factor_array:
        h_region, n_region = scipy.ndimage.measurements.label((h_g * cutcut) > (area * characteristic_density * factor))
        #print 'factor', factor, n_region, n_region < 10
        if n_region < 10:
            threshold_density = area * characteristic_density * factor
            break

    h_region, n_region = scipy.ndimage.measurements.label((h_g * cutcut) > threshold_density)
    h_region = np.ma.array(h_region, mask=(h_region < 1))

    x_peak_array = []
    y_peak_array = []
    angsep_peak_array = []

    for index in range(1, n_region + 1): # loop over peaks
        index_peak = np.argmax(h_g * (h_region == index))
        x_peak, y_peak = xx.flatten()[index_peak], yy.flatten()[index_peak]
        #print index, np.max(h_g * (h_region == index))
        
        #angsep_peak = np.sqrt((x_full - x_peak)**2 + (y_full - y_peak)**2) # Use full magnitude range, NOT TESTED!!!
        angsep_peak = np.sqrt((x - x_peak)**2 + (y - y_peak)**2) # Impose magnitude threshold

        x_peak_array.append(x_peak)
        y_peak_array.append(y_peak)
        angsep_peak_array.append(angsep_peak)
    
    return x_peak_array, y_peak_array, angsep_peak_array

########################################################################

def fit_aperture(proj, distance_modulus, characteristic_density_local, x_peak, y_peak, angsep_peak):
    """
    Fit aperture by varing radius and computing the significance
    """

    ra_peak_array = []
    dec_peak_array = []
    r_peak_array = []
    sig_peak_array = []
    distance_modulus_array = []
    n_obs_peak_array = []
    n_obs_half_peak_array = []
    n_model_peak_array = []

    size_array = np.arange(0.01, 0.3, 0.01)
    sig_array = np.tile(0., len(size_array))
    
    size_array_zero = np.concatenate([[0.], size_array])
    area_array = np.pi * (size_array_zero[1:]**2 - size_array_zero[0:-1]**2)

    n_obs_array = np.tile(0, len(size_array))
    n_model_array = np.tile(0., len(size_array))
    for ii in range(0, len(size_array)):
        n_obs = np.sum(angsep_peak < size_array[ii])
        n_model = characteristic_density_local * (np.pi * size_array[ii]**2)
        sig_array[ii] = np.clip(scipy.stats.norm.isf(scipy.stats.poisson.sf(n_obs, n_model)), 0., 37.5) # Clip at 37.5
        n_obs_array[ii] = n_obs
        n_model_array[ii] = n_model

    ra_peak, dec_peak = proj.imageToSphere(x_peak, y_peak)

    index_peak = np.argmax(sig_array)
    r_peak = size_array[index_peak]
    #if np.max(sig_array) >= 37.5:
    #    r_peak = 0.5
    n_obs_peak = n_obs_array[index_peak]
    n_model_peak = n_model_array[index_peak]
    n_obs_half_peak = np.sum(angsep_peak < (0.5 * r_peak))

    # Compile resilts
    print('Candidate: x_peak: {:12.3f}, y_peak: {:12.3f}, r_peak: {:12.3f}, sig: {:12.3f}, ra_peak: {:12.3f}, dec_peak: {:12.3f}'.format(x_peak, y_peak, r_peak, np.max(sig_array), ra_peak, dec_peak))
    ra_peak_array.append(ra_peak)
    dec_peak_array.append(dec_peak)
    r_peak_array.append(r_peak)
    #sig_peak_array.append(np.max(sig_array))
    sig_peak_array.append(sig_array[index_peak])
    distance_modulus_array.append(distance_modulus)
    n_obs_peak_array.append(n_obs_peak)
    n_obs_half_peak_array.append(n_obs_half_peak)
    n_model_peak_array.append(n_model_peak)

    return ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array

########################################################################

# mode = 0
def search_by_distance(nside, data, distance_modulus, pix_nside_select, ra_select, dec_select, magnitude_threshold=mag_max, fracdet=None):
    """
    Idea: 
    Send a data extension that goes to faint magnitudes, e.g., g < 24.
    Use the whole region to identify hotspots using a slightly brighter 
    magnitude threshold, e.g., g < 23, so not susceptible to variations 
    in depth. Then compute the local field density using a small annulus 
    around each individual hotspot, e.g., radius 0.3 to 0.5 deg.

    fracdet corresponds to a fracdet map (numpy array, assumed to be EQUATORIAL and RING)
    """

    print('Distance = {:0.1f} kpc (m-M = {:0.1f})').format(ugali.utils.projector.distanceModulusToDistance(distance_modulus), distance_modulus)

    iso = ugali.isochrone.factory(name=isoname, survey=isosurvey, band_1=band_1.lower(), band_2=band_2.lower())
    iso.age = 12.
    iso.metallicity = 0.0001
    iso.distance_modulus = distance_modulus

    cut = cut_isochrone_path(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1)
    data = data[cut]

    print('{} objects left after isochrone cut...').format(len(data))

    if (len(data) == 0):
        return [], [], [], [], [], [], [], []

    # Compute characteristic density at this distance
    characteristic_density = compute_char_density(nside, data, ra_select, dec_select, mag_max, fracdet)

    ra_peak_array = []
    dec_peak_array = []
    r_peak_array = []
    sig_peak_array = []
    distance_modulus_array = []
    n_obs_peak_array = []
    n_obs_half_peak_array = []
    n_model_peak_array = []

    proj = ugali.utils.projector.Projector(ra_select, dec_select)

    x_peak_array, y_peak_array, angsep_peak_array = find_peaks(nside, data, characteristic_density, distance_modulus, pix_nside_select, ra_select, dec_select, magnitude_threshold, fracdet)

    for x_peak, y_peak, angsep_peak in itertools.izip(x_peak_array, y_peak_array, angsep_peak_array):
        characteristic_density_local = compute_local_char_density(nside, data, characteristic_density, ra_select, dec_select, x_peak, y_peak, angsep_peak, mag_max, fracdet)
        # Aperture fitting
        print('Fitting aperture to hotspot...')
        ra_peaks, dec_peaks, r_peaks, sig_peaks, distance_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = fit_aperture(proj, distance_modulus, characteristic_density_local, x_peak, y_peak, angsep_peak)
        
        ra_peak_array.append(ra_peaks)
        dec_peak_array.append(dec_peaks)
        r_peak_array.append(r_peaks)
        sig_peak_array.append(sig_peaks)
        distance_modulus_array.append(distance_moduli)
        n_obs_peak_array.append(n_obs_peaks)
        n_obs_half_peak_array.append(n_obs_half_peaks)
        n_model_peak_array.append(n_model_peaks)

    ra_peak_array = np.concatenate(ra_peak_array)
    dec_peak_array = np.concatenate(dec_peak_array)
    r_peak_array = np.concatenate(r_peak_array)
    sig_peak_array = np.concatenate(sig_peak_array)
    distance_modulus_array = np.concatenate(distance_modulus_array)
    n_obs_peak_array = np.concatenate(n_obs_peak_array)
    n_obs_half_peak_array = np.concatenate(n_obs_half_peak_array)
    n_model_peak_array = np.concatenate(n_model_peak_array)

    return ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array

########################################################################

# mode = 1
def search_by_simulation(nside, data, distance_modulus, pix_nside_select, ra_select, dec_select, magnitude_threshold=mag_max, fracdet=None):
    """
    Idea: 
    Send a data extension that goes to faint magnitudes, e.g., g < 24.
    Use the whole region to identify hotspots using a slightly brighter 
    magnitude threshold, e.g., g < 23, so not susceptible to variations 
    in depth. Then compute the local field density using a small annulus 
    around each individual hotspot, e.g., radius 0.3 to 0.5 deg.

    fracdet corresponds to a fracdet map (numpy array, assumed to be EQUATORIAL and RING)
    """

    print('Distance = {:0.1f} kpc (m-M = {:0.1f})').format(ugali.utils.projector.distanceModulusToDistance(distance_modulus), distance_modulus)

    iso = ugali.isochrone.factory(name=isoname, survey=isosurvey, band_1=band_1.lower(), band_2=band_2.lower())
    iso.age = 12.
    iso.metallicity = 0.0001
    iso.distance_modulus = distance_modulus

    cut = cut_isochrone_path(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1)
    data = data[cut]

    print('{} objects left after isochrone cut...'.format(len(data)))
    print('{} simulated objects left after isochrone cut...'.format(np.sum(data['MC_SOURCE_ID'] != 0)))

    if (len(data) == 0):
        return [], [], [], [], [], [], [], []

    # Compute characteristic density at this distance
    characteristic_density = compute_char_density(nside, data, ra_select, dec_select, mag_max, fracdet)

    ra_peak_array = []
    dec_peak_array = []
    r_peak_array = []
    sig_peak_array = []
    distance_modulus_array = []
    n_obs_peak_array = []
    n_obs_half_peak_array = []
    n_model_peak_array = []

    proj = ugali.utils.projector.Projector(ra_select, dec_select)

    cut_magnitude_threshold = (data[mag_dered_1] < magnitude_threshold)
    x_peak, y_peak = proj.sphereToImage(ra_select, dec_select) # = 0, 0
    x, y = proj.sphereToImage(data[basis_1][cut_magnitude_threshold], data[basis_2][cut_magnitude_threshold])
    angsep_peak = np.sqrt((x - x_peak)**2 + (y - y_peak)**2)

    characteristic_density_local = compute_local_char_density(nside, data, characteristic_density, ra_select, dec_select, x_peak, y_peak, angsep_peak, mag_max, fracdet)

    # Aperture fitting
    print('Fitting aperture to hotspot...')
    ra_peaks, dec_peaks, r_peaks, sig_peaks, distance_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = fit_aperture(proj, distance_modulus, characteristic_density_local, x_peak, y_peak, angsep_peak)
    
    ra_peak_array.append(ra_peaks)
    dec_peak_array.append(dec_peaks)
    r_peak_array.append(r_peaks)
    sig_peak_array.append(sig_peaks)
    distance_modulus_array.append(distance_moduli)
    n_obs_peak_array.append(n_obs_peaks)
    n_obs_half_peak_array.append(n_obs_half_peaks)
    n_model_peak_array.append(n_model_peaks)

    ra_peak_array = np.concatenate(ra_peak_array)
    dec_peak_array = np.concatenate(dec_peak_array)
    r_peak_array = np.concatenate(r_peak_array)
    sig_peak_array = np.concatenate(sig_peak_array)
    distance_modulus_array = np.concatenate(distance_modulus_array)
    n_obs_peak_array = np.concatenate(n_obs_peak_array)
    n_obs_half_peak_array = np.concatenate(n_obs_half_peak_array)
    n_model_peak_array = np.concatenate(n_model_peak_array)

    return ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array

########################################################################

# mode = 2
def search_by_object(nside, data, distance_modulus, pix_nside_select, ra_select, dec_select, magnitude_threshold=mag_max, fracdet=None):
    """
    Idea: 
    Send a data extension that goes to faint magnitudes, e.g., g < 24.
    Use the whole region to identify hotspots using a slightly brighter 
    magnitude threshold, e.g., g < 23, so not susceptible to variations 
    in depth. Then compute the local field density using a small annulus 
    around each individual hotspot, e.g., radius 0.3 to 0.5 deg.

    fracdet corresponds to a fracdet map (numpy array, assumed to be EQUATORIAL and RING)
    """

    print('Distance = {:0.1f} kpc (m-M = {:0.1f})').format(ugali.utils.projector.distanceModulusToDistance(distance_modulus), distance_modulus)

    iso = ugali.isochrone.factory(name=isoname, survey=isosurvey, band_1=band_1.lower(), band_2=band_2.lower())
    iso.age = 12.
    iso.metallicity = 0.0001
    iso.distance_modulus = distance_modulus

    cut = cut_isochrone_path(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1)
    data = data[cut]

    print('{} objects left after isochrone cut...'.format(len(data)))
    print('{} simulated objects left after isochrone cut...'.format(np.sum(data['MC_SOURCE_ID'] != 0)))

    if (len(data) == 0):
        return [], [], [], [], [], [], [], []

    # Compute characteristic density at this distance
    characteristic_density = compute_char_density(nside, data, ra_select, dec_select, mag_max, fracdet)

    ra_peak_array = []
    dec_peak_array = []
    r_peak_array = []
    sig_peak_array = []
    distance_modulus_array = []
    n_obs_peak_array = []
    n_obs_half_peak_array = []
    n_model_peak_array = []

    proj = ugali.utils.projector.Projector(ra_select, dec_select)

    cut_magnitude_threshold = (data[mag_dered_1] < magnitude_threshold)
    x_peak, y_peak = proj.sphereToImage(ra_select, dec_select) # = 0, 0
    x, y = proj.sphereToImage(data[basis_1][cut_magnitude_threshold], data[basis_2][cut_magnitude_threshold])
    angsep_peak = np.sqrt((x - x_peak)**2 + (y - y_peak)**2)

    characteristic_density_local = compute_local_char_density(nside, data, characteristic_density, ra_select, dec_select, x_peak, y_peak, angsep_peak, mag_max, fracdet)

    # Aperture fitting
    print('Fitting aperture to hotspot...')
    ra_peaks, dec_peaks, r_peaks, sig_peaks, distance_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = fit_aperture(proj, distance_modulus, characteristic_density_local, x_peak, y_peak, angsep_peak)
    
    ra_peak_array.append(ra_peaks)
    dec_peak_array.append(dec_peaks)
    r_peak_array.append(r_peaks)
    sig_peak_array.append(sig_peaks)
    distance_modulus_array.append(distance_moduli)
    n_obs_peak_array.append(n_obs_peaks)
    n_obs_half_peak_array.append(n_obs_half_peaks)
    n_model_peak_array.append(n_model_peaks)

    ra_peak_array = np.concatenate(ra_peak_array)
    dec_peak_array = np.concatenate(dec_peak_array)
    r_peak_array = np.concatenate(r_peak_array)
    sig_peak_array = np.concatenate(sig_peak_array)
    distance_modulus_array = np.concatenate(distance_modulus_array)
    n_obs_peak_array = np.concatenate(n_obs_peak_array)
    n_obs_half_peak_array = np.concatenate(n_obs_half_peak_array)
    n_model_peak_array = np.concatenate(n_model_peak_array)

    return ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array

########################################################################

def write_output(results_dir, nside, pix_nside_select, ra_peak_array, dec_peak_array, r_peak_array, distance_modulus_array, 
                n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, 
                sig_peak_array, mc_source_id_array, mode, outfile):
    #writer = open(outfile, 'a') # append if exists
    #for ii in range(0, len(sig_peak_array)):
    #    # SIG, RA, DEC, MODULUS, r, n_obs, n_model, mc_source_id
    #    writer.write('{:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}\n'.format(sig_peak_array[ii], 
    #                                                                                                                     ra_peak_array[ii], 
    #                                                                                                                     dec_peak_array[ii], 
    #                                                                                                                     distance_modulus_array[ii], 
    #                                                                                                                     r_peak_array[ii],
    #                                                                                                                     n_obs_peak_array[ii],
    #                                                                                                                     n_obs_half_peak_array[ii],
    #                                                                                                                     n_model_peak_array[ii],
    #                                                                                                                     mc_source_id_array[ii]))
    data = [tuple(row) for row in np.stack([sig_peak_array, ra_peak_array, dec_peak_array, distance_modulus_array, r_peak_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, mc_source_id_array], axis=-1)]
    arr = np.array(data, dtype=[('SIG', float), ('RA', float), ('DEC', float), ('MODULUS', float), ('R', float), ('N_OBS', float), ('N_OBS_HALF', float), ('N_MODEL', float), ('MC_SOURCE_ID', float)])
    np.save(outfile, arr)

########################################################################

def read_output(results_dir, pix):
    #infiles = sorted(glob.glob(results_dir + '/*{}.csv'.format(pix)))
    #if (len(infiles) > 0):
    #    results = [np.genfromtxt(infile, delimiter=',', names=['SIG', 'RA', 'DEC', 'DISTANCE_MODULUS', 'R_PEAK', 'N_OBS_PEAK', 'N_OBS_HALF_PEAK', 'N_MODEL_PEAK', 'MC_SOURCE_ID']) for infile in infiles]
    #    results = [result for result in results if (result.size != 0)]
    #    #results = np.concatenate(np.vstack(results))
    #    results = np.concatenate(results)
    infile = glob.glob(results_dir + '/*{}.csv'.format(pix))
    if (len(infile) > 0):
        results = np.genfromtxt(infile[0], delimiter=',', names=['SIG', 'RA', 'DEC', 'MODULUS', 'R', 'N_OBS', 'N_OBS_HALF', 'N_MODEL', 'MC_SOURCE_ID'])
    else:
        results = []

    return results
