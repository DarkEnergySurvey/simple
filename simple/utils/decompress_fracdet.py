import fitsio
import astropy.io as fits
import yaml
import healpy as hp

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    survey = cfg['survey']
    nside   = cfg[survey]['nside']
    candidate_list = cfg[survey]['candidate_list']
    fracdet_gz = cfg[survey]['fracdet_gz']
    fracdet = cfg[survey]['fracdet']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']

# Read in fracdet map
if (fracdet_gz is not None) and (fracdet_gz.lower().strip() != 'none') and (fracdet_gz != ''):
    print('Reading fracdet map {} ...').format(fracdet_gz)
    fracdet_map = hp.read_map(fracdet_gz)
    hp.write_map(fracdet, fracdet_map, nest=False) # TODO
    #fits.write(fracdet_map) # TODO
else:
    print('No fracdet map specified ...')
