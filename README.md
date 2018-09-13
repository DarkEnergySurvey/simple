# Simple Binning Suite 

The simple binning suite includes modules for peforming a simple binning search of stellar overdensities in DES &c data, compiling a candidate list from the results of the search, and producing diagnostic plots from a candidate list.

This code has been adapted from [Keith Bechtol](https://github.com/bechtol)'s original simple binning program to be more modular and ready-to-use for future searches.

Maintained by [Sidney Mau](https://github.com/SidneyMau).

## Dependencies

This program uses the following:

python;
numpy,
scipy,
healpy,
astropy,
matplotlib,
[fitsio](https://github.com/esheldon/fitsio)

[ugali](https://github.com/DarkEnergySurvey/ugali),
[batchtools](https://github.com/kadrlica/batchtools) (for farming)

and some others that I've forgotten to write down.

## Configuration and use

`config.yaml` handles most of the setup and can be modified according to use.
You will need to modify the value for the `simple_dir` key to be the path to your `simple` folder that contains these source files.

You will need to create a directory such as `simple_run/`.
Copy `config.yaml` into `simple_run/`; this will let the the scripts know where to source the config info from as well as where to write the output.
If you have fracdet or maglim files, those should also be copied or linked here; specify their filepaths in `config.yaml`.

To run the simple binning search, run `farm_simple.py` in `simple_run/`.
This will run `search_algorithm.py` over the given data set and write the output to `results_dir/`, logging each job in `results_dir/log_dir`.
The results can then be compiled into a candidate list by running `make_list.py` from `simple_run/` (this is saved as a `.fits` file).

To produce plots from a candidate list, run `farm_plots.py` in `simple_run/`.
The output will be written to `save_dir/` with logs in `save_dir/log_dir/`.

## Notes

By default, `farm_plots.py` will only produce plots for hotspots with statistical significance greater than 5.5 sigma.
This threshold has intentionally chosen to be low (such that investigations can be made of very low-significance hotsposts and candidates) but also to minimize junk.

This code is currently going through an architectural overhaul to capitalize on the OOP nature of python.
