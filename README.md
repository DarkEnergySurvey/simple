# Simple Binning Search 

The simple binning suite includes modules for peforming a simple binning search of stellar overdensities in DES &c data, compiling a candidate list from the resulsts of the search, and producing diagnostic plots from a candidate list.

This code has been adapted from [Keith Bechtol](https://github.com/bechtol)'s original simple binning program to be more modular and ready-to-use for future searches.

Maintained by [Sidney Mau](https://github.com/SidneyMau).

## Dependencies

python

numpy

scipy

matplotlib

healpy

[fitsio](https://github.com/esheldon/fitsio)

[ugali](https://github.com/DarkEnergySurvey/ugali)

## Configuration and use

`config.yaml` handles most of the setup and can be modified according to use. You will need to modify the value for the `simple_dir` key to be the path to your `simple` folder that contains these source files.

First, you will need to create a directory such as `simple_run/` (the name can be anything). In `simple_run/`, create a symbolic link to `config.yaml` with `ln -s /path/to/simple/config.yaml`. This will let the the scripts know where to source the config info from as well as where to write the output. If you have fracdet or maglim files, those should also be copied or linked here; otherwise, they may be specified in `config.yaml`.

To run the simple binning search, run `farm_simple.py` from `simple_run/`. This will run `search_algorithm.py` over the given data set and write the output to `results_dir/`, logging each job in `results_dir/log_dir`. The results can then be compiled into a candidate list by running `make_list.py` from `simple_run/`.

To render plots from a candidate list, run `farm_plots.py` from `simple_run/`. The output will be written to `save_dir/` with logs in `save_dir/log_dir/`.
