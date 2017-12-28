#!/usr/bin/env bash

# Generate simple_dir from config.yaml
echo "Sourcing config..."
simple_dir=$(grep 'simple_dir' config.yaml | awk 'NF>1{print $NF}')
eval simple_dir=$simple_dir

# Run farm_simple.py
echo "Farming simple..."
farm_simple="$simple_dir/farm_simple.py"
#python $farm_simple

# wait...
#echo "Waiting for search to finish..."
#wait
#echo "Search finished."

# Run make_list.py
echo "Compiling candidate list..."
make_list="$simple_dir/make_list.py"
python $make_list

# Run farm_plots.py
echo "Farming plots..."
farm_plots="$simple_dir/farm_plots.py"
python $farm_plots

# TODO: Does not wait for all the cjobs to finish before continuing.
#       Need for all parallel jobs to complete before moving on.
# wait...
echo "Waiting for plots to finish..."
wait
echo "Plots finished."

# Tarball deliverable
echo "Tarring..."
save_dir=$(grep 'save_dir' config.yaml | awk 'NF>1{print $NF}')
eval save_dir=$save_dir
tar -czf simple_plots.tar.gz $save_dir
