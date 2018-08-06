#!/usr/bin/env python
"""
Diagnostic plot functions
"""
__author__ = "Sidney Mau"

import os
import glob
import yaml

import fitsio as fits
from astropy.coordinates import SkyCoord
from ugali.utils import healpix
from ugali.isochrone import factory as isochrone_factory
import ugali.utils.projector
import ugali.utils.plotting
import healpy

import numpy as np
from operator import add
from scipy import interpolate
from scipy.signal import argrelextrema
import scipy.ndimage

import pylab as plt
import matplotlib
from matplotlib import mlab
from mpl_toolkits.axes_grid1 import make_axes_locatable

import simple.filters
import simple.simple_utils

################################################################################

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
    
    mag_g = cfg[survey]['mag_g']
    mag_r = cfg[survey]['mag_r']
    mag_g_err = cfg[survey]['mag_g_err']
    mag_r_err = cfg[survey]['mag_r_err']

    band_1 = cfg[survey]['band_1']
    band_2 = cfg[survey]['band_2']
    mag = cfg[survey]['mag']
    mag_err = cfg[survey]['mag_err']
    mag_dered = cfg[survey]['mag_dered']

# construct mags
mag_1 = mag.format(band_1.upper())
mag_2 = mag.format(band_2.upper())
mag_err_1 = mag_err.format(band_1.upper())
mag_err_2 = mag_err.format(band_2.upper())
mag_dered_1 = mag_dered.format(band_1.upper())
mag_dered_2 = mag_dered.format(band_2.upper())
    
################################################################################

def analysis(targ_ra, targ_dec, mod, mc_source_id):
    """Analyze a candidate"""

    pix_nside_select = ugali.utils.healpix.angToPix(nside, targ_ra, targ_dec)
    ra_select, dec_select = ugali.utils.healpix.pixToAng(nside, pix_nside_select)
    pix_nside_neighbors = np.concatenate([[pix_nside_select], healpy.get_all_neighbours(nside, pix_nside_select)])

    # Construct data
    #data = simple_utils.construct_modal_data(mode, pix_nside_neighbors)
    data = simple_utils.construct_real_data(pix_nside_neighbors)
    if (mode == 0):
        print('mode = 0: running only on real data')
    elif (mode == 1):
        print('mode = 1: running on real data and simulated data')
        
        # inject objects for simulated object of mc_source_id
        sim_data = simple_utils.construct_sim_data(pix_nside_neighbors, mc_source_id)
        data = simple_utils.inject_sim(data, sim_data, mc_source_id)
    else:
        print('No mode specified; running only on real data')

    print('Loading data...')
    data = simple_utils.construct_modal_data(mode, pix_nside_neighbors, mc_source_id)
    quality_cut = filters.quality_filter(survey, data)
    data = data[quality_cut]
    print('Found {} objects...').format(len(data))

    data = filters.dered_mag(survey, data)

    # This should be generalized to also take the survey
    iso = isochrone_factory(name=isoname, survey=isosurvey, age=12, z=0.0001, distance_modulus=mod, band_1=band_1.lower(), band_2=band_2.lower())

    # g_radius estimate
    filter = filters.star_filter(survey, data)

    iso_filter = simple_utils.cutIsochronePath(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1, return_all=False)

    angsep = ugali.utils.projector.angsep(targ_ra, targ_dec, data[basis_1], data[basis_2])

    bins = np.linspace(0, 0.4, 21) # deg
    centers = 0.5*(bins[1:] + bins[0:-1])
    area = np.pi*(bins[1:]**2 - bins[0:-1]**2) * 60**2
    hist = np.histogram(angsep[(angsep < 0.4) & filter & iso_filter], bins=bins)[0] # counts

    f_interp = interpolate.interp1d(np.linspace(centers[0], centers[-1], len(hist)), hist/area, 'cubic')
    f_range = np.linspace(centers[0], centers[-1], 1000)
    f_val = f_interp(f_range)

    pairs = zip(f_range, f_val)

    peak = max(pairs[:len(pairs)/4], key=lambda x: x[1]) # find peak within first quarter

    def peak_index(pairs, peak):
        for i in range(len(pairs)):
            if pairs[i] == peak:
                return i

    osc = int(0.04/0.4*1000) # +/- 0.04 (rounded down) deg oscillation about local extremum
    relmin = argrelextrema(f_val, np.less, order=osc)[0]

    try:
        if len(relmin) > 0:
            #half_point = f_range[relmin[0]]
            i = 0
            while ((f_range[relmin[i]] <= f_range[peak_index(pairs,peak)]) & (i <= len(relmin)-1)):
                i+=1
            half_point = f_range[relmin[i]]
        elif len(relmin) == 0:
            half_peak = (peak[1] + np.mean(f_val[len(f_val)/4:]))/2. # normalized to background (after first quarter)
            #half_peak = np.mean(f_val[len(f_val)/4:])
            half_pairs = []
            for i in pairs[peak_index(pairs, peak):len(pairs)/2]: # start after peak, stay within first quarter
                if i != peak:
                    half_pairs.append((i[0], abs(i[1]-half_peak)))
            half_point = min(half_pairs, key=lambda x: x[1])[0] # deg
    except:
        half_point = 0.1 # fixed value to catch errors

    g_min = 0.5/60. # deg
    g_max = 12./60. # deg

    if half_point < g_min:
        g_radius = g_min
    elif half_point > g_max:
        g_radius = g_max
    else:
        g_radius = half_point # deg

    angsep = ugali.utils.projector.angsep(targ_ra, targ_dec, data[basis_1], data[basis_2])
    nbhd = (angsep < g_radius)

    return(data, iso, g_radius, nbhd)

def densityPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd, type):
    """Stellar density plot"""

    if type == 'stars':
        filter = filters.star_filter(survey, data)
        plt.title('Stellar Density')
    elif type == 'galaxies':
        filter = filters.galaxy_filter(survey, data)
        plt.title('Galactic Density')
    elif type == 'blue_stars':
        filter = filters.color_filter(survey, data) \
               & filters.star_filter(survey, data)
        plt.title('Blue Stellar Density')
    
    iso_filter = simple_utils.cutIsochronePath(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1, return_all=False)

    # projection of image
    proj = ugali.utils.projector.Projector(targ_ra, targ_dec)
    x, y = proj.sphereToImage(data[filter & iso_filter][basis_1], data[filter & iso_filter][basis_2])

    bound = 0.5 #1.
    steps = 100.
    bins = np.linspace(-bound, bound, steps)

    signal = np.histogram2d(x, y, bins=[bins, bins])[0]

    sigma = 0.01 * (0.25 * np.arctan(0.25*g_radius*60. - 1.5) + 1.3)
    
    convolution = scipy.ndimage.filters.gaussian_filter(signal, sigma/(bound/steps))
    plt.pcolormesh(bins, bins, convolution.T, cmap='Greys')

    plt.xlim(bound, -bound)
    plt.ylim(-bound, bound)
    plt.gca().set_aspect('equal')
    plt.xlabel(r'$\Delta \alpha$ (deg)')
    plt.ylabel(r'$\Delta \delta$ (deg)')

    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '5%', pad=0)
    plt.colorbar(cax=cax)

def starPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd):
    """Star bin plot"""

    filter = filters.star_filter(survey, data)

    iso_filter = simple_utils.cutIsochronePath(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1, return_all=False)

    # projection of image
    proj = ugali.utils.projector.Projector(targ_ra, targ_dec)
    x, y = proj.sphereToImage(data[filter & iso_filter][basis_1], data[filter & iso_filter][basis_2])

    plt.scatter(x, y, edgecolor='none', s=3, c='black')
    plt.xlim(0.2, -0.2)
    plt.ylim(-0.2, 0.2)
    plt.gca().set_aspect('equal')
    plt.xlabel(r'$\Delta \alpha$ (deg)')
    plt.ylabel(r'$\Delta \delta$ (deg)')

    plt.title('Stars')

def cmPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd, type):
    """Color-magnitude plot"""

    angsep = ugali.utils.projector.angsep(targ_ra, targ_dec, data[basis_1], data[basis_2])
    annulus = (angsep > g_radius) & (angsep < 1.)

    if type == 'stars':
        filter = filters.star_filter(survey, data)
        plt.title('Stellar Color-Magnitude')
    elif type == 'galaxies':
        filter = filters.galaxy_filter(survey, data)
        plt.title('Galactic Color-Magnitude')

    iso_filter = simple_utils.cutIsochronePath(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1, return_all=False)

    # Plot background objects
    plt.scatter(data[mag_dered_1][filter & annulus] - data[mag_dered_2][filter & annulus], data[mag_dered_1][filter & annulus], c='k', alpha=0.1, edgecolor='none', s=1)

    # Plot isochrone
    ugali.utils.plotting.drawIsochrone(iso, lw=2, label='{} Gyr, z = {}'.format(iso.age, iso.metallicity))

    # Plot objects in nbhd
    plt.scatter(data[mag_dered_1][filter & nbhd] - data[mag_dered_2][filter & nbhd], data[mag_dered_2][filter & nbhd], c='g', s=5, label='r < {:.3f}$^\circ$'.format(g_radius))

    # Plot objects in nbhd and near isochrone
    plt.scatter(data[mag_dered_1][filter & nbhd & iso_filter] - data[mag_dered_2][filter & nbhd & iso_filter], data[mag_dered_1][filter & nbhd & iso_filter], c='r', s=5, label='$\Delta$CM < 0.1')

    plt.axis([-0.5, 1, 16, mag_max])
    plt.gca().invert_yaxis()
    plt.gca().set_aspect(1./4.)
    plt.legend(loc='upper left')
    plt.xlabel('{} - {} (mag)'.format(band_1.lower(), band_2.lower()))
    plt.ylabel('{} (mag)'.format(band_1.lower()))

def hessPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd):
    """Hess plot"""

    filter = filters.star_filter(survey, data)

    plt.title('Hess')

    c1 = SkyCoord(targ_ra, targ_dec, unit='deg')

    r_near = 2.*g_radius # annulus begins at 2*g_radius away from centroid
    r_far = np.sqrt(5.)*g_radius # annulus has same area as inner area

    #inner = (c1.separation(SkyCoord(data[basis_1], data[basis_2], unit='deg')).deg < g_radius)
    #outer = (c1.separation(SkyCoord(data[basis_1], data[basis_2], unit='deg')).deg > r_near) & (c1.separation(SkyCoord(data[basis_1], data[basis_2], unit='deg')).deg < r_far)
    angsep = ugali.utils.projector.angsep(targ_ra, targ_dec, data[basis_1], data[basis_2])
    inner = (angsep < g_radius)
    outer = ((angsep > r_near) & (angsep < r_far))

    xbins = np.arange(-0.5, 1.1, 0.1)
    ybins = np.arange(16., mag_max + 0.5, 0.5)

    foreground = np.histogram2d(data[mag_dered_1][inner & filter] - data[mag_dered_2][inner & filter], data[mag_dered_1][inner & filter], bins=[xbins, ybins])
    background = np.histogram2d(data[mag_dered_1][outer & filter] - data[mag_dered_2][outer & filter], data[mag_dered_1][outer & filter], bins=[xbins, ybins])

    fg = foreground[0].T
    bg = background[0].T

    fg_abs = np.absolute(fg)
    bg_abs = np.absolute(bg)

    mask_abs = fg_abs + bg_abs
    mask_abs[mask_abs == 0.] = np.nan # mask signficiant zeroes

    signal = fg - bg
    signal = np.ma.array(signal, mask=np.isnan(mask_abs)) # mask nan

    cmap = matplotlib.cm.viridis
    cmap.set_bad('w', 1.)
    plt.pcolormesh(xbins, ybins, signal, cmap=cmap)

    plt.colorbar()

    ugali.utils.plotting.drawIsochrone(iso, lw=2, c='k', zorder=10, label='Isocrhone')

    plt.axis([-0.5, 1.0, 16, mag_max])
    plt.gca().invert_yaxis()
    plt.gca().set_aspect(1./4.)
    plt.xlabel('{} - {} (mag)'.format(band_1.lower(), band_2.lower()))
    plt.ylabel('{} (mag)'.format(band_1.lower()))

def radialPlot(targ_ra, targ_dec, data, iso, g_radius, nbhd, field_density=None):
    """Radial distribution plot"""

    ## Deprecated?
    #filter_s = filters.star_filter(survey, data)
    #filter_g = filters.galaxy_filter(survey, data)

    plt.title('Radial Distribution')

    angsep = ugali.utils.projector.angsep(targ_ra, targ_dec, data[basis_1], data[basis_2])

    # Isochrone filtered/unfiltered
    iso_seln_f = simple_utils.cutIsochronePath(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1, return_all=False)
    iso_seln_u = ~iso_seln_f

    bins = np.linspace(0, 0.4, 21) # deg
    centers = 0.5*(bins[1:] + bins[0:-1])
    area = np.pi*(bins[1:]**2 - bins[0:-1]**2) * 60**2

    def interp_values(type, seln):
        if type == 'stars':
            filter = filters.star_filter(survey, data)
        elif type == 'galaxies':
            filter = filters.galaxy_filter(survey, data)

        if seln == 'f':
            iso_filter = iso_seln_f
        elif seln == 'u':
            iso_filter = iso_seln_u

        hist = np.histogram(angsep[(angsep < 0.4) & filter & iso_filter], bins=bins)[0] # counts

        f_interp = interpolate.interp1d(np.linspace(centers[0], centers[-1], len(hist)), hist/area, 'cubic')
        f_range = np.linspace(centers[0], centers[-1], 1000)
        f_val = f_interp(f_range)

        return(f_range, f_val)

    def value_errors(type, seln):
        if type == 'stars':
            filter = filters.star_filter(survey, data)
        elif type == 'galaxies':
            filter = filters.galaxy_filter(survey, data)
        if seln == 'f':
            iso_filter = iso_seln_f
        elif seln == 'u':
            iso_filter = iso_seln_u

        hist = np.histogram(angsep[(angsep < 0.4) & filter & iso_filter], bins=bins)[0] # counts

        val = hist/area
        yerr = np.sqrt(hist)/area

        return(val, yerr)

    f_range, f_val = interp_values('stars', 'f')
    pairs = zip(f_range, f_val)
    peak = max(pairs[:len(pairs)/4], key=lambda x: x[1]) # find peak within first quarter
    def peak_index(pairs, peak):
        for i in range(len(pairs)):
            if pairs[i] == peak:
                return i

    plt.axvline(x=f_range[peak_index(pairs,peak)], color='m', label='peak')

    plt.axvline(x=g_radius, color='r', label='g_radius')

    f_range, f_val = interp_values('galaxies', 'f')
    plt.plot(f_range, f_val, '-g', label='Filtered Galaxies')

    f_range, f_val = interp_values('stars', 'u')
    plt.plot(f_range, f_val, '-k', alpha=0.25, label='Unfiltered Stars')

    val, y_err = value_errors('stars', 'f')
    plt.plot(centers, val, '.b')
    plt.errorbar(centers, val, yerr=y_err, fmt='none', ecolor='b', elinewidth=1, capsize=5)

    f_range, f_val = interp_values('stars', 'f')
    plt.plot(f_range, f_val, '-b', label='Filtered Stars')

    if field_density:
        plt.axhline(y=field_density, color='blue', ls='--', label='Model Field')

    ymax = plt.ylim()[1]
    plt.annotate(r'$\approx %0.1f$' + str(round(g_radius, 3)) + '$^\circ$', (g_radius*1.1, ymax/50.), color='red', bbox=dict(boxstyle='round,pad=0.0', fc='white', alpha=0.75, ec='white', lw=0))
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., ymax)
    plt.legend(loc='upper right')
    plt.xlabel('Angular Separation (deg)')
    plt.ylabel('Denisty (arcmin$^{-2})$')

#def maglim_plot(targ_ra, targ_dec, data, iso, band):
#    """Maglim plots"""
#
#    reso = 0.5
#    xsize = 2.*60./reso
#
#    if band == 'g':
#        reader = pyfits.open(maglim_g)
#        m_maglim_g = reader[1].data.field('I').flatten()
#        reader.close()
#        m_maglim_g[np.isnan(m_maglim_g)] = healpy.UNSEEN
#        #healpy.gnomview(m_maglim_g, fig='summary', rot=(targ_ra, targ_dec, 0.), reso=reso, xsize=xsize, title='maglim g (S/N =10)', sub=(3, 4, 11))
#        healpy.gnomview(m_maglim_g, rot=(targ_ra, targ_dec, 0.), reso=reso, xsize=xsize, title='maglim g (S/N =10)', sub=(3, 4, 8))
#    elif band == 'r':
#        reader = pyfits.open(maglim_r)
#        m_maglim_r = reader[1].data.field('I').flatten()
#        reader.close()
#        m_maglim_r[np.isnan(m_maglim_r)] = healpy.UNSEEN
#        healpy.gnomview(m_maglim_r, rot=(targ_ra, targ_dec, 0.), reso=reso, xsize=xsize, title='maglim r (S/N =10)', sub=(3, 4, 12))