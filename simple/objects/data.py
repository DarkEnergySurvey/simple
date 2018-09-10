#!/usr/bin/env python
"""
Classy data
"""
__author__ = "Sid Mau"

# Python libraries
import os
import glob
import yaml
import numpy as np
import numpy.lib.recfunctions
import healpy as hp
import fitsio as fits
import scipy.interpolate
import scipy.ndimage
import itertools

# Ugali libraries
import ugali.utils.healpix
import ugali.utils.projector
import ugali.isochrone

########################################################################

class Data:
    """
    Class object for analyzing photometric data.
    """
    def __init__(self, survey, nside, dirname, fracdet, band_1, band_2, mag, mag_err, mag_dered, basis_1, basis_2, mag_max):
        self.survey      = survey
        self.nside       = nside
        self.dirname     = dirname
        self.fracdet     = fracdet
        self.band_1      = band_1
        self.band_2      = band_2
        self.mag_1       = mag.format(band_1.upper())
        self.mag_2       = mag.format(band_2.upper())
        self.mag_err_1   = mag_err.format(band_1.upper())
        self.mag_err_2   = mag_err.format(band_2.upper())
        self.mag_dered_1 = mag_dered.format(band_1.upper())
        self.mag_dered_2 = mag_dered.format(band_2.upper())
        self.basis_1     = basis_1
        self.basis_2     = basis_2
        self.mag_max     = mag_max

    def load_local_data(self, ra, dec):
        """
        Load data corresponding to the 8 nearest neighbors of the
        centroid at (ra, dec) into memory.
        """
        hpixel = ugali.utils.healpix.angToPix(self.nside, ra, dec)
        hpixels = np.concatenate([[hpixel], hp.get_all_neighbours(self.nside, hpixel)])

        data_array = []
        for hpix in hpixels:
            inlist = glob.glob('{}/*_{:05d}.fits'.format(self.dirname, hpix))
            for infile in inlist:
                if not os.path.exists(infile):
                    continue
                data_array.append(fits.read(infile))
        data = np.concatenate(data_array)

        return data

    def compute_char_density(self, ra, dec):
        """
        Compute the characteristic density of a region.
        """
    
        data = self.load_local_data(ra, dec)

        mag_cut = (data[self.mag_1] < self.mag_max)
    
        proj = ugali.utils.projector.Projector(ra, dec)
        x, y = proj.sphereToImage(data[self.basis_1][mag_cut], data[self.basis_2][mag_cut]) # Trimmed magnitude range for hotspot finding
        delta_x = 0.01
        area = delta_x**2
        smoothing = 2. / 60. # Was 3 arcmin
        bins = np.arange(-8., 8. + 1.e-10, delta_x)
        centers = 0.5 * (bins[0: -1] + bins[1:])
        yy, xx = np.meshgrid(centers, centers)
    
        h = np.histogram2d(x, y, bins=[bins, bins])[0]
    
        h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)
    
        delta_x_coverage = 0.1
        area_coverage = (delta_x_coverage)**2
        bins_coverage = np.arange(-5., 5. + 1.e-10, delta_x_coverage)
        h_coverage = np.histogram2d(x, y, bins=[bins_coverage, bins_coverage])[0]
        h_goodcoverage = np.histogram2d(x, y, bins=[bins_coverage, bins_coverage])[0]
    
        n_goodcoverage = h_coverage[h_goodcoverage > 0].flatten()
    
        characteristic_density = np.median(n_goodcoverage) / area_coverage # per square degree
        print('Characteristic density = {:0.1f} deg^-2').format(characteristic_density)
    
        # Use pixels with fracdet ~1.0 to estimate the characteristic density
        if self.fracdet is not None:
            fracdet = fits.read(self.fracdet)
            fracdet_zero = np.tile(0., len(fracdet))
            cut = (fracdet != hp.UNSEEN)
            fracdet_zero[cut] = fracdet[cut]
    
            nside_fracdet = hp.npix2nside(len(fracdet))
            
            subpix_region_array = []
            for pix in np.unique(ugali.utils.healpix.angToPix(self.nside, data[self.basis_1], data[self.basis_2])):
                subpix_region_array.append(ugali.utils.healpix.subpixel(pix, self.nside, nside_fracdet))
            subpix_region_array = np.concatenate(subpix_region_array)
    
            # Compute mean fracdet in the region so that this is available as a correction factor
            cut = (fracdet[subpix_region_array] != hp.UNSEEN)
            mean_fracdet = np.mean(fracdet[subpix_region_array[cut]])
    
            # smau: this doesn't seem to be used in the non-local density estimation
            subpix_region_array = subpix_region_array[fracdet[subpix_region_array] > 0.99]
            subpix = ugali.utils.healpix.angToPix(nside_fracdet, 
                                                  data[self.basis_1][mag_cut], 
                                                  data[self.basis_2][mag_cut]) # Remember to apply mag threshold to objects
            characteristic_density_fracdet = float(np.sum(np.in1d(subpix, subpix_region_array))) \
                                             / (hp.nside2pixarea(nside_fracdet, degrees=True) * len(subpix_region_array)) # deg^-2
            print('Characteristic density fracdet = {:0.1f} deg^-2').format(characteristic_density_fracdet)
            
            # Correct the characteristic density by the mean fracdet value
            characteristic_density_raw = 1. * characteristic_density
            characteristic_density /= mean_fracdet 
            print('Characteristic density (fracdet corrected) = {:0.1f} deg^-2').format(characteristic_density)
    
        return characteristic_density
    
    def compute_local_char_density(self, ra, dec, x_peak, y_peak, angsep_peak):
        """
        Compute the local characteristic density of a region.
        """
    
        data = self.load_local_data(ra, dec)

        characteristic_density = self.compute_char_density(ra, dec)

        mag_cut = (data[mag_dered_1] < self.mag_max)
    
        proj = ugali.utils.projector.Projector(ra, dec)
        x, y = proj.sphereToImage(data[basis_1][mag_cut], data[basis_2][mag_cut]) # Trimmed magnitude range for hotspot finding
    
        # If fracdet map is available, use that information to either compute local density,
        # or in regions of spotty coverage, use the typical density of the region
        if self.fracdet is not None:
            fracdet = fits.read(self.fracdet)
            fracdet_zero = np.tile(0., len(fracdet))
            cut = (fracdet != hp.UNSEEN)
            fracdet_zero[cut] = fracdet[cut]
    
            nside_fracdet = hp.npix2nside(len(fracdet))
            
            subpix_region_array = []
            for pix in np.unique(ugali.utils.healpix.angToPix(nside, data[basis_1], data[basis_2])):
                subpix_region_array.append(ugali.utils.healpix.subpixel(pix, self.nside, nside_fracdet))
            subpix_region_array = np.concatenate(subpix_region_array)
    
            # Compute mean fracdet in the region so that this is available as a correction factor
            cut = (fracdet[subpix_region_array] != hp.UNSEEN)
            mean_fracdet = np.mean(fracdet[subpix_region_array[cut]])
    
            subpix_region_array = subpix_region_array[fracdet[subpix_region_array] > 0.99]
            subpix = ugali.utils.healpix.angToPix(nside_fracdet, 
                                                  data[self.basis_1][mag_cut], 
                                                  data[self.basis_2][mag_cut]) # Remember to apply mag threshold to objects
    
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


    def find_peaks(self, ra, dec, distance_modulus):
        """
        Convolve field to find characteristic density and peaks within the selected pixel.
        """

        pix_nside_select = ugali.utils.healpix.angToPix(self.nside, ra, dec)
        data = self.load_local_data(ra, dec)
        characteristic_density = self.compute_char_density(ra, dec)
        mag_cut = (data[mag_dered_1] < self.mag_max)
    
        # convolve field and find peaks
        proj = ugali.utils.projector.Projector(ra, dec)
        x, y = proj.sphereToImage(data[self.basis_1][mag_cut], data[self.basis_2][mag_cut]) # Trimmed magnitude range for hotspot finding
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
        cutcut = (ugali.utils.healpix.angToPix(self.nside, rara, decdec) == pix_nside_select).reshape(xx.shape)
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

    def fit_aperture(self, proj, distance_modulus, characteristic_density_local, x_peak, y_peak, angsep_peak):
        """
        Fit aperture by varing radius and computing the significance.
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
