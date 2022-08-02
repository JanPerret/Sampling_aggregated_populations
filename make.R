#########################################################
# 
# Sampling simulation of spatially aggregated populations
# 
# make.R
# 
# jan.perret@gmail.com
#########################################################

### 
### Master script to launch all the other scripts.
### Be careful, some scripts require a lot of computation power and/or time.
### Consider increasing the number of cores used for parallel computing in each script.
### 

# clean workspace
rm(list = ls())

# estimate the values of the Matern cluster process parameters to get the desired levels of dispersion index for modality 1 (increasing density by increasing both kappa and mu)
source("./1_estimate_point_process_parameters_Matern_modality_1.R")

# same thing for modality 2 (increasing density by increasing only kappa)
source("./1_estimate_point_process_parameters_Matern_modality_2.R")

# same thing for modality 3 (increasing density by increasing only mu)
source("./1_estimate_point_process_parameters_Matern_modality_3.R")

# estimate the values of the SSI process parameters to get the desired levels of dispersion index
source("./1_estimate_point_process_parameters_SSI.R")

# create the point patterns for the aggregated virtual populations (Matern cluster process with modality 1)
source("./2_pattern_creation_Matern_modality_1.R")

# create the point patterns for the aggregated virtual populations (Matern cluster process with modalities 2 and 3)
source("./2_pattern_creation_Matern_modalities_2_3.R")

# create the point patterns for the additional densities (5, 15, 25 and 30 individuals/cell) for the aggregated virtual populations (Matern cluster process with modality 1)
source("./2_pattern_creation_Matern_additional_intensities.R")

# create the point patterns for the virtual populations with repulsion between individuals (SSI process)
source("./2_pattern_creation_SSI.R")

# simulate the sampling surveys for the virtual populations (Matern clustered populations with modality 1, poisson populations and SSI populations)
source("./3_sampling_simulation_modality_1.R")

# simulate the sampling surveys for the aggregated virtual populations (Matern clustered populations with modality 2)
source("./3_sampling_simulation_modality_2.R")

# simulate the sampling surveys for the aggregated virtual populations (Matern clustered populations with modality 3)
source("./3_sampling_simulation_modality_3.R")

# simulate the sampling surveys for the additional densities for the aggregated virtual populations (Matern clustered populations with modality 1)
source("./3_sampling_simulation_additional_intensities.R")

# simulate the sampling surveys for the three natural populations
source("./3_sampling_simulation_natural_populations.R")

# simulate the sampling surveys to test the performance of variance estimators for the virtual populations
source("./3_sampling_simulation_variance_estimators_virtual_populations.R")

# simulate the sampling surveys to test the performance of variance estimators for the three natural populations
source("./3_sampling_simulation_variance_estimators_natural_populations.R")

# compute the average distance between the sampling units for Spatially Balanced Sampling (SBS)
source("./4_compute_distance_between_sampling_units_SBS.R")

# merge all the files resulting from the simulations
source("./5_merge_result_files.R")

# make the figures for the article
source("./6_Perret_et_al_FIGURES.R")


