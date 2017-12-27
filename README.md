# Simple Binning Search Program

This code has been adapted from Keith Bechtol's original simple binning code program to be more modular and ready-to-use for future searches.

[Sidney Mau](https://github.com/SidneyMau)

## Configuration and use

<!--The `config.yaml` file will point the different parts of the code to the data set over which the search will be performed and the directories where the results will be saved. Once this has been set, then `farm_simple.py` will use `search_algorithm.py` to perform the simple binnning search over the given data set, storing the results in `results_dir/`. Then, `make_list.py` will compile the data stored in `results_dir/` into a candidate list. Finally, running `farm_plots.py` will use `diagnostic_plots.py` to create diagnostic plots for each candidate, saving these plots in `save_dir/`.-->

`config.yaml` handles most of the setup and can be modified according to use. You will need to modify the value for the `simple_dir` key to be the path to your `simple` folder that contains these source files.

First, you will need to create a directory such as `simple_run/` (the name can be anything). In `simple_run/`, create a symbolic link to `config.yaml` with `ln -s /path/to/simple/config.yaml`. This will let the the scripts know where to source the config info from as well as where to write the output.

To run the simple binning search, run `farm_simple.py` from `simple_run/`. This will run `search_algorithm.py` over the given data set and write the output to `results_dir/`, logging each job in `results_dir/log_dir`. The resulsts can then be compiled into a candidate list by running `make_list.py` from `simple_run/`.

To render plots from a candidate list, run `farm_plots.py` from `simple_run/`. The output will be written to `save_dir/` with logs in `save_dir/log_dir/`.
