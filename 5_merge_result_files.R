#######################################################
# 
# Sampling simulation of spatially aggregated populations
# 
# 5_merge_result_files.R
# 
# jan.perret@gmail.com
#######################################################

###
### script objective : merge all csv files resulting from the sampling survey simulations and save a version averaged over the 60 patterns
###

### clean workspace
rm(list = ls())

### load packages
library(data.table)
library(dplyr)
library(readr)



#### load simulation results for MODALITY 1 ####
folder_path <- "./output/text/RUN_2022-06-08_simul_aggregation_intensity_results_modality_1/"

# load and merge the simulation output csv files
filenames <- paste0(folder_path, list.files(path = folder_path, pattern = "*.csv"))
files <- lapply(filenames, read.csv2)
merged_data_mod1 <- data.table::rbindlist(files, use.names = FALSE)

# sort the table
merged_data_mod1 <- merged_data_mod1[with(merged_data_mod1, order(theoretical_intensity, expect_disp_ind, sample_size, num_pattern)), ]

# join the cluster diameter
merged_data_mod1 <- cbind(merged_data_mod1, cluster_diam = merged_data_mod1$my_scale * 2)

# # save the dataframe
# write.csv2(merged_data_mod1, file = "./data/tidy/RUN_2022-06-08_simul_aggregation_intensity_results_merged_patterns_1-60.csv", row.names = FALSE)

# make a version of the dataframe averaged over the 60 patterns
merged_data_mod1 <- merged_data_mod1 %>% select(-num_pattern)

merged_data_mod1_averaged <- merged_data_mod1 %>% 
  group_by(pattern, theoretical_intensity, sample_size) %>% 
  summarise(across(my_kappa:cluster_diam, ~ mean(.x, na.rm = TRUE)))

# save the dataframe
write.csv2(merged_data_mod1_averaged, file = "./data/tidy/RUN_2022-06-08_simul_aggregation_intensity_results_merged_patterns_1-60_averaged.csv", row.names = FALSE)

### 25 minutes on my personal computer



#### load simulation results for the additional intensities of MODALITY 1 ####
folder_path <- "./output/text/RUN_2022-06-08_simul_aggregation_intensity_results_additional_intensities/"

# load and merge the simulation output csv files
filenames <- paste0(folder_path, list.files(path = folder_path, pattern = "*.csv"))
files <- lapply(filenames, read.csv2)
merged_data_add_intens <- data.table::rbindlist(files, use.names = FALSE)

# sort the table
merged_data_add_intens <- merged_data_add_intens[with(merged_data_add_intens, order(theoretical_intensity, expect_disp_ind, sample_size, num_pattern)), ]

# join the cluster diameter
merged_data_add_intens <- cbind(merged_data_add_intens, cluster_diam = merged_data_add_intens$my_scale * 2)

# # save the dataframe
# write.csv2(merged_data_add_intens, file = "./data/tidy/RUN_2022-06-08_simul_aggregation_intensity_results_additional_intensities_merged_patterns_1-60.csv", row.names = FALSE)

# save a version of the dataframe averaged over all the patterns
merged_data_add_intens <- merged_data_add_intens %>% select(-num_pattern)

merged_data_add_intens_averaged <- merged_data_add_intens %>% 
  group_by(pattern, theoretical_intensity, sample_size) %>% 
  summarise(across(my_kappa:cluster_diam, ~ mean(.x, na.rm = TRUE)))

# save the dataframe
write.csv2(merged_data_add_intens_averaged,
           file = "./data/tidy/RUN_2022-06-08_simul_aggregation_intensity_results_additional_intensities_merged_patterns_1-60_averaged.csv", row.names = FALSE)

### 47 minutes on my personal computer



#### load simulation results for MODALITY 2 ####
folder_path <- "./output/text/RUN_2022-06-08_simul_aggregation_intensity_results_modality_2/"

# load and merge the simulation output csv files
filenames <- paste0(folder_path, list.files(path = folder_path, pattern = "*.csv"))
files <- lapply(filenames, read.csv2)
merged_data_mod2 <- data.table::rbindlist(files, use.names = FALSE)

# sort the table
merged_data_mod2 <- merged_data_mod2[with(merged_data_mod2, order(theoretical_intensity, expect_disp_ind, sample_size, num_pattern)), ]

# join the cluster diameter
merged_data_mod2 <- cbind(merged_data_mod2, cluster_diam = merged_data_mod2$my_scale * 2)

# # save the dataframe
# write.csv2(merged_data_mod2, file = "./data/tidy/RUN_2022-06-08_simul_aggregation_intensity_results_modality_2_merged.csv", row.names = FALSE)

# save a version of the dataframe averaged over all the patterns
merged_data_mod2 <- merged_data_mod2 %>% select(-num_pattern)

merged_data_mod2_averaged <- merged_data_mod2 %>% 
  group_by(pattern, theoretical_intensity, sample_size) %>% 
  summarise(across(my_kappa:cluster_diam, ~ mean(.x, na.rm = TRUE)))

# join dispersion index levels 0.4, 0.6, 0.8 and 1 (for intensities 1, 5, 10, 20)
merged_data_mod1_averaged <- rbind(merged_data_mod1_averaged, merged_data_averaged_additional_intensities_total)
vect_patterns_to_join <- c("int1_disp0.4", "int1_disp0.6", "int1_disp0.8", "int1_disp1", "int5_disp0.4", "int5_disp0.6", "int5_disp0.8", "int5_disp1",
                           "int10_disp0.4", "int10_disp0.6", "int10_disp0.8", "int10_disp1", "int20_disp0.4", "int20_disp0.6", "int20_disp0.8", "int20_disp1")
merged_data_mod1_averaged_for_join <- subset(merged_data_mod1_averaged, merged_data_mod1_averaged$pattern %in% vect_patterns_to_join)

merged_data_mod2_averaged <- rbind(merged_data_mod2_averaged, merged_data_mod1_averaged_for_join)

# save the dataframe
write.csv2(merged_data_mod2_averaged,
           file = "./data/tidy/RUN_2022-06-08_simul_aggregation_intensity_results_modality_2_merged_averaged.csv", row.names = FALSE)



#### load simulation results for MODALITY 3 ####
folder_path <- "./output/text/RUN_2022-06-08_simul_aggregation_intensity_results_modality_3/"

# load and merge the simulation output csv files
filenames <- paste0(folder_path, list.files(path = folder_path, pattern = "*.csv"))
files <- lapply(filenames, read.csv2)
merged_data_mod3 <- data.table::rbindlist(files, use.names = FALSE)

# sort the table
merged_data_mod3 <- merged_data_mod3[with(merged_data_mod3, order(theoretical_intensity, expect_disp_ind, sample_size, num_pattern)), ]

# join the cluster diameter
merged_data_mod3 <- cbind(merged_data_mod3, cluster_diam = merged_data_mod3$my_scale * 2)

# # save the dataframe
# write.csv2(merged_data_mod3, file = "./data/tidy/RUN_2022-06-08_simul_aggregation_intensity_results_modality_3_merged.csv", row.names = FALSE)

# save a version of the dataframe averaged over all the patterns
merged_data_mod3 <- merged_data_mod3 %>% select(-num_pattern)

merged_data_mod3_averaged <- merged_data_mod3 %>% 
  group_by(pattern, theoretical_intensity, sample_size) %>% 
  summarise(across(my_kappa:cluster_diam, ~ mean(.x, na.rm = TRUE)))

# join dispersion index levels 0.4, 0.6, 0.8 and 1 (for intensities 1, 5, 10, 20)
merged_data_mod3_averaged <- rbind(merged_data_mod3_averaged, merged_data_mod1_averaged_for_join)

# save the dataframe
write.csv2(merged_data_averaged_modality_3,
           file = "./data/tidy/RUN_2022-06-08_simul_aggregation_intensity_results_modality_3_merged_averaged.csv", row.names = FALSE)


