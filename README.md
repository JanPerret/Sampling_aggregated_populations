
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6954499.svg)](https://doi.org/10.5281/zenodo.6954499)

# Sampling_aggregated_populations
Data and code for a paper on the relative precision of sampling methods for aggregated populations

This repository contains all the code to run the sampling survey simulations of spatially aggregated populations used in a paper that is actually under review. Please contact me if something goes wrong with the code, or for any question.

Mail: [jan.perret@cefe.cnrs.fr](mailto:jan.perret@cefe.cnrs.fr) or [jan.perret@gmail.com](mailto:jan.perret@gmail.com)

## Steps you need to follow to run the code

1. Download the current folder.

2. Download the ZIP file "output.zip" from the Zenodo link below, and unzip it in the root directory. <br>
Zenodo link : [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6954484.svg)](https://doi.org/10.5281/zenodo.6954484) <br>
I have stored this file separately from the current directory because it contains the virtual populations generated in advance which are over 1 Gb. <br>

3. If you want to run all the scripts, you have to run them in the order of the numbers at the beginning of the file names. You can also execute the master script called "make.R". The scripts starting with 1 can be run in any order, but they must all have been executed before moving on to the scripts starting with 2, etc. Be careful, running all the scripts will require a lot of computing power and/or time. Consider increasing the number of cores used for parallel computing in each script. Alternatively, you can run any script independently of the others. All the files resulting from the execution of the scripts are already in the "output" folder, so the input files needed for any of the scripts are already available.

## Description of what each script does

**"1_estimate_point_process_parameters[...].R"** : The scripts with names starting like this are used to estimate the values of the point process parameters needed to get the wanted dispersion index levels for the virtual populations. Indeed, we used the dispersion index to measure the aggregation of populations, and since we wanted a gradient of aggregation we had to estimate the point process parameter values allowing to get the dispersion index values we wanted. 

**"2_pattern_creation_[...].R"** : The scripts with names starting like this use the point process parameter values obtained in the previous step to create the virtual populations (i.e. point patterns) on which the sampling survey simulations will be carried out.


**"3_sampling_simulation_[...].R"** : The scripts with names starting like this perform the sampling survey simulations on the virtual populations created in the previous step or on the three natural populations.


**"4_compute_distance_between_sampling_units_SBS.R"** : This script is used to compute the average distance between the sampling units for SBS. Indeed, for SYS the average distance between sample units can easily be calculated by hand, but for SBS it requires running simulations.

**"5_merge_result_files.R"** : This script loads all the files resulting from the simulations and merges them into single files that can be used to make the plots.

**"6_Perret_et_al_FIGURES.R"** : This script makes the plots for the figures of the article.


