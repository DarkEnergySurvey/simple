import os
import sys
import yaml
import numpy as np
import fitsio as fits
import matplotlib.pyplot as plt
import ugali.candidate.associate

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

    simple_dir = cfg['setup']['simple_dir']
    jobs = cfg['batch']['jobs']
    survey = cfg['survey']
    nside   = cfg[survey]['nside']
    datadir = cfg[survey]['datadir']
    candidate_list = cfg[survey]['candidate_list']
    basis_1 = cfg[survey]['basis_1']
    basis_2 = cfg[survey]['basis_2']

save_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'])
if not os.path.exists(save_dir):
    os.mkdir(save_dir)

log_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'], cfg['output']['log_dir'])
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

data = fits.read(candidate_list)

try:
    sig_sel, plot = float(sys.argv[1]), float(sys.argv[2])
except:
    sys.exit('ERROR! Give sig selection (0---37.5) and plot choice (0 or 1).')

data.sort(order='SIG')

#x = [360 - ra if ra > 180 else ra for ra in data[basis_1]]
#y = data[basis_2]

plt.scatter(data[basis_1], data[basis_2], c=data['SIG'], s=1)
plt.colorbar(label='SIG')

# Overplot known objects
catalog_array = ['McConnachie15', 'ExtraDwarfs']
catalog = ugali.candidate.associate.SourceCatalog()
for catalog_name in catalog_array:
    catalog += ugali.candidate.associate.catalogFactory(catalog_name)
plt.scatter(catalog['ra'], catalog['dec'], edgecolor='k', s=20, linewidth=0.5, c='none', label='known objects')

for hotspot in data[data['SIG'] > sig_sel]:
    distances = [np.sqrt((hotspot[basis_1] - c['ra'])**2 + (hotspot[basis_2] - c['dec'])**2) for c in catalog]
    if np.min(distances) > 1:
        plt.scatter(hotspot[basis_1], hotspot[basis_2], edgecolor='r', s=20, linewidth=0.5, c='none')
        if plot == 1:
            sig     = round(hotspot['SIG'], 2)
            ra      = round(hotspot[basis_1], 2)
            dec     = round(hotspot[basis_2], 2)
            mod     = round(hotspot['MODULUS'], 2)
            mc_source_id = round(hotspot['MC_SOURCE_ID'], 2)
            if 'N_MODEL' in data.dtype.names:
                field_density = round(hotspot['N_MODEL'] / (np.pi * (hotspot['r'] * 60.)**2), 4) # field density (arcmin^-2)
            
            logfile = '{}/candidate_{}_{}.log'.format(log_dir, ra, dec)
            #batch = 'csub -n {} -o {} '.format(jobs, logfile)
            batch = 'csub -n {} -o {} --host all '.format(jobs, logfile) # testing condor updates
            if 'N_MODEL' in data.dtype.names:
                command = 'python {}/make_plot.py {:0.2f} {:0.2f} {:0.2f} {:0.2f} {:0.2f} {:0.4f}'.format(simple_dir, ra, dec, mod, sig, mc_source_id, field_density)
            else:
                command = 'python {}/make_plot.py {:0.2f} {:0.2f} {:0.2f} {:0.2f} {:0.2f}'.format(simple_dir, ra, dec, mod, sig, mc_source_id)
            command_queue = batch + command

            print(command_queue)
            os.system(command_queue) # Submit to queue


plt.ylim([min(data[basis_2]), max(data[basis_2])])
plt.gca().invert_xaxis()
plt.xlabel(basis_1)
plt.ylabel(basis_2)
plt.title("{} Hotspots (SIG > {})".format(survey, sig_sel))
plt.savefig("{}_hotspots_{}.png".format(survey, sig_sel), bbox_inches='tight', dpi=300)
plt.close()

