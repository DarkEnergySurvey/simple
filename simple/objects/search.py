#!/usr/bin/env python
"""
Search algorithm
"""
__author__ = "Sid Mau"

# Python libraries
import numpy as np
import itertools

# Ugali libraries
import ugali.utils.projector
import ugali.isochrone

# Simple libraries
import simple.simple_utils
import simple.objects.point
import simple.objects.result

# TODO:
# - use point
# - make search functions

class Search:
    """
    Class for searching.
    """
    def __init__(self, data):
        self.data    = data
        self.nside   = data.nside
        self.band_1  = data.band_1
        self.band_2  = data.band_2
        self.basis_1 = data.basis_1
        self.basis_2 = data.basis_2

    def search_by_distance(self, nside, data, distance_modulus, ra, dec):
        """
        Idea: 
        Send a data extension that goes to faint magnitudes, e.g., g < 24.
        Use the whole region to identify hotspots using a slightly brighter 
        magnitude threshold, e.g., g < 23, so not susceptible to variations 
        in depth. Then compute the local field density using a small annulus 
        around each individual hotspot, e.g., radius 0.3 to 0.5 deg.
        """
    
        print('Distance = {:0.1f} kpc (m-M = {:0.1f})').format(ugali.utils.projector.distanceModulusToDistance(distance_modulus),
                                                               distance_modulus)

        point = simple.objects.point.Point(ra, dec)
        pix_nside_select = point.healpixel(nside)
        proj = point.projector
    
        iso = ugali.isochrone.factory(name=isoname,
                                      survey=isosurvey,
                                      age=12.,
                                      metallicity=0.0001,
                                      distance_modulus = distance_modulus,
                                      band_1=self.band_1.lower(),
                                      band_2=self.band_2.lower())
    
        cut = simple.simple_utils.cut_isochrone_path(data[data.mag_dered_1], data[data.mag_dered_2], data[data.mag_err_1], data[data.mag_err_2], iso, radius=0.1)
        data = data[cut]
    
        print('{} objects left after isochrone cut...').format(len(data))
    
        characteristic_density = data.compute_char_density(ra, dec)
    
        search_result = simple.objects.result.Result()
    
        fracdet = data.load_fracdet
        x_peak_array, y_peak_array, angsep_peak_array = data.find_peaks(ra, dec, distance_modulus)
    
        for x_peak, y_peak, angsep_peak in itertools.izip(x_peak_array, y_peak_array, angsep_peak_array):
            characteristic_density_local = data.compute_local_char_density(ra, dec, x_peak, y_peak, angsep_peak)
            # Aperture fitting
            print('Fitting aperture to hotspot...')
            results = data.fit_aperture(ra, dec, proj, distance_modulus, x_peak, y_peak, angsep_peak)
            search_result.append_results(*results)
    
        search_result.concatenate_results()
    
        return search_result

    ## mode = 1
    #def search_by_simulation(nside, data, distance_modulus, pix_nside_select, ra_select, dec_select, magnitude_threshold=mag_max, fracdet=None):
    #    """
    #    Idea: 
    #    Send a data extension that goes to faint magnitudes, e.g., g < 24.
    #    Use the whole region to identify hotspots using a slightly brighter 
    #    magnitude threshold, e.g., g < 23, so not susceptible to variations 
    #    in depth. Then compute the local field density using a small annulus 
    #    around each individual hotspot, e.g., radius 0.3 to 0.5 deg.
    #
    #    fracdet corresponds to a fracdet map (numpy array, assumed to be EQUATORIAL and RING)
    #    """
    #
    #    print('Distance = {:0.1f} kpc (m-M = {:0.1f})').format(ugali.utils.projector.distanceModulusToDistance(distance_modulus), distance_modulus)
    #
    #    iso = ugali.isochrone.factory(name=isoname, survey=isosurvey, band_1=band_1.lower(), band_2=band_2.lower())
    #    iso.age = 12.
    #    iso.metallicity = 0.0001
    #    iso.distance_modulus = distance_modulus
    #
    #    cut = simple.simple_utils.cut_isochrone_path(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1)
    #    data = data[cut]
    #
    #    print('{} objects left after isochrone cut...'.format(len(data)))
    #    print('{} simulated objects left after isochrone cut...'.format(np.sum(data['MC_SOURCE_ID'] != 0)))
    #
    #    # Compute characteristic density at this distance
    #    characteristic_density = compute_char_density(nside, data, ra_select, dec_select, mag_max, fracdet)
    #
    #    ra_peak_array = []
    #    dec_peak_array = []
    #    r_peak_array = []
    #    sig_peak_array = []
    #    distance_modulus_array = []
    #    n_obs_peak_array = []
    #    n_obs_half_peak_array = []
    #    n_model_peak_array = []
    #
    #    proj = ugali.utils.projector.Projector(ra_select, dec_select)
    #
    #    cut_magnitude_threshold = (data[mag_dered_1] < magnitude_threshold)
    #    x_peak, y_peak = proj.sphereToImage(ra_select, dec_select) # = 0, 0
    #    x, y = proj.sphereToImage(data[basis_1][cut_magnitude_threshold], data[basis_2][cut_magnitude_threshold])
    #    angsep_peak = np.sqrt((x - x_peak)**2 + (y - y_peak)**2)
    #
    #    characteristic_density_local = compute_local_char_density(nside, data, characteristic_density, ra_select, dec_select, x_peak, y_peak, angsep_peak, mag_max, fracdet)
    #
    #    # Aperture fitting
    #    print('Fitting aperture to hotspot...')
    #    ra_peaks, dec_peaks, r_peaks, sig_peaks, distance_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = fit_aperture(proj, distance_modulus, characteristic_density_local, x_peak, y_peak, angsep_peak)
    #    
    #    ra_peak_array.append(ra_peaks)
    #    dec_peak_array.append(dec_peaks)
    #    r_peak_array.append(r_peaks)
    #    sig_peak_array.append(sig_peaks)
    #    distance_modulus_array.append(distance_moduli)
    #    n_obs_peak_array.append(n_obs_peaks)
    #    n_obs_half_peak_array.append(n_obs_half_peaks)
    #    n_model_peak_array.append(n_model_peaks)
    #
    #    ra_peak_array = np.concatenate(ra_peak_array)
    #    dec_peak_array = np.concatenate(dec_peak_array)
    #    r_peak_array = np.concatenate(r_peak_array)
    #    sig_peak_array = np.concatenate(sig_peak_array)
    #    distance_modulus_array = np.concatenate(distance_modulus_array)
    #    n_obs_peak_array = np.concatenate(n_obs_peak_array)
    #    n_obs_half_peak_array = np.concatenate(n_obs_half_peak_array)
    #    n_model_peak_array = np.concatenate(n_model_peak_array)
    #
    #    return ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array
    #
    #########################################################################
    #
    ## mode = 2
    #def search_by_object(nside, data, distance_modulus, pix_nside_select, ra_select, dec_select, magnitude_threshold=mag_max, fracdet=None):
    #    """
    #    Idea: 
    #    Send a data extension that goes to faint magnitudes, e.g., g < 24.
    #    Use the whole region to identify hotspots using a slightly brighter 
    #    magnitude threshold, e.g., g < 23, so not susceptible to variations 
    #    in depth. Then compute the local field density using a small annulus 
    #    around each individual hotspot, e.g., radius 0.3 to 0.5 deg.
    #
    #    fracdet corresponds to a fracdet map (numpy array, assumed to be EQUATORIAL and RING)
    #    """
    #
    #    print('Distance = {:0.1f} kpc (m-M = {:0.1f})').format(ugali.utils.projector.distanceModulusToDistance(distance_modulus), distance_modulus)
    #
    #    iso = ugali.isochrone.factory(name=isoname, survey=isosurvey, band_1=band_1.lower(), band_2=band_2.lower())
    #    iso.age = 12.
    #    iso.metallicity = 0.0001
    #    iso.distance_modulus = distance_modulus
    #
    #    cut = simple.simple_utils.cut_isochrone_path(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1)
    #    data = data[cut]
    #
    #    print('{} objects left after isochrone cut...'.format(len(data)))
    #    print('{} simulated objects left after isochrone cut...'.format(np.sum(data['MC_SOURCE_ID'] != 0)))
    #
    #    # Compute characteristic density at this distance
    #    characteristic_density = compute_char_density(nside, data, ra_select, dec_select, mag_max, fracdet)
    #
    #    ra_peak_array = []
    #    dec_peak_array = []
    #    r_peak_array = []
    #    sig_peak_array = []
    #    distance_modulus_array = []
    #    n_obs_peak_array = []
    #    n_obs_half_peak_array = []
    #    n_model_peak_array = []
    #
    #    proj = ugali.utils.projector.Projector(ra_select, dec_select)
    #
    #    cut_magnitude_threshold = (data[mag_dered_1] < magnitude_threshold)
    #    x_peak, y_peak = proj.sphereToImage(ra_select, dec_select) # = 0, 0
    #    x, y = proj.sphereToImage(data[basis_1][cut_magnitude_threshold], data[basis_2][cut_magnitude_threshold])
    #    angsep_peak = np.sqrt((x - x_peak)**2 + (y - y_peak)**2)
    #
    #    characteristic_density_local = compute_local_char_density(nside, data, characteristic_density, ra_select, dec_select, x_peak, y_peak, angsep_peak, mag_max, fracdet)
    #
    #    # Aperture fitting
    #    print('Fitting aperture to hotspot...')
    #    ra_peaks, dec_peaks, r_peaks, sig_peaks, distance_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = fit_aperture(proj, distance_modulus, characteristic_density_local, x_peak, y_peak, angsep_peak)
    #    
    #    ra_peak_array.append(ra_peaks)
    #    dec_peak_array.append(dec_peaks)
    #    r_peak_array.append(r_peaks)
    #    sig_peak_array.append(sig_peaks)
    #    distance_modulus_array.append(distance_moduli)
    #    n_obs_peak_array.append(n_obs_peaks)
    #    n_obs_half_peak_array.append(n_obs_half_peaks)
    #    n_model_peak_array.append(n_model_peaks)
    #
    #    ra_peak_array = np.concatenate(ra_peak_array)
    #    dec_peak_array = np.concatenate(dec_peak_array)
    #    r_peak_array = np.concatenate(r_peak_array)
    #    sig_peak_array = np.concatenate(sig_peak_array)
    #    distance_modulus_array = np.concatenate(distance_modulus_array)
    #    n_obs_peak_array = np.concatenate(n_obs_peak_array)
    #    n_obs_half_peak_array = np.concatenate(n_obs_half_peak_array)
    #    n_model_peak_array = np.concatenate(n_model_peak_array)
    #
    #    return ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array
    #
