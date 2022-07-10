#########################################################
# 
# Sampling simulation of spatially aggregated populations
# 
# 6_Perret_et_al_FIGURES.R
# 
# jan.perret@gmail.com
#########################################################

###
### script objective : make the figures for the article
###

# clean workspace
rm(list = ls())

# load packages
library(ggplot2)
library(ggpubr)
library(spatstat)
library(dplyr)
library(tidyr)

# load functions
files.sources <- list.files("./R", full.names = TRUE)
invisible(sapply(files.sources, source))

# load data
simul_result <- read.csv2(file = "./data/tidy/RUN_2022-06-08_simul_aggregation_intensity_results_merged_patterns_1-60_averaged.csv")
simul_result_additional_intensities <- read.csv2(file = "./data/tidy/RUN_2022-06-08_simul_aggregation_intensity_results_additional_intensities_merged_patterns_1-60_averaged.csv")
simul_result <- rbind(simul_result, simul_result_additional_intensities)

df_SBS_dist_neighbours <- read.csv2(file = "./data/tidy/BAS_mean_distance_nearest_units_14_neighbours.csv")

# convert sample_size to factor
simul_result$sample_size <- factor(simul_result$sample_size, levels = c("9", "15", "25", "49", "100", "150", "196", "300", "400"))
df_SBS_dist_neighbours$sample_size <- factor(df_SBS_dist_neighbours$sample_size, levels = c("9", "15", "25", "49", "100", "150", "196", "300", "400"))



##### FIGURE 1 : workflow diagram #####

# 1. plot 4 patterns with increasing aggregation for one intensity level  
# simulate the patterns to be plotted
my_window <- owin(c(0, 20), c(0, 20))

set.seed(20210921)
pattern1 <- rMatClust(kappa = 0.02236068, scale = 4, mu = 223.6068, win = my_window)
pattern2 <- rMatClust(kappa = 0.02236068, scale = 2.6, mu = 223.6068, win = my_window)
pattern3 <- rMatClust(kappa = 0.02236068, scale = 1.5, mu = 223.6068, win = my_window)
pattern_poisson <- rpoispp(lambda = 5, win = my_window)

### function to plot the point patterns
plot_mypattern2 <- function(mypattern, title) {
  
  # convert ppp object to dataframe for plotting with ggplot2
  pattern_df <- data.frame(x = mypattern$x, y = mypattern$y)
  
  # plot the points with ggplot
  myplot <- ggplot(pattern_df, aes(x = x, y = y)) + 
    geom_point(size = 0.4, alpha = 0.8, color = "black") + 
    xlim(min = my_window$xrange[1], max = my_window$xrange[2]) +
    ylim(min = my_window$yrange[1], max = my_window$yrange[2]) +
    coord_fixed() +
    theme_bw() +
    theme(panel.background = element_rect(colour = "black", fill = "NA",
                                          linetype = "solid", size = 0.3),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.line = element_blank(), plot.background = element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.y = element_blank(), axis.text.y = element_blank())
  
  return(myplot)
}

plot_pattern1 <- plot_mypattern2(mypattern = pattern1)
plot_pattern2 <- plot_mypattern2(mypattern = pattern2)
plot_pattern3 <- plot_mypattern2(mypattern = pattern3)
plot_pattern_poisson <- plot_mypattern2(mypattern = pattern_poisson)

# arrange the plots on a single line
workflow_patterns <- ggarrange(plot_pattern_poisson, plot_pattern1, plot_pattern2, plot_pattern3, 
                               ncol = 4, nrow = 1)

# save the plot
pdf(file = "./output/plot/Workflow_patterns.pdf", width = 7, height = 2.5)
print(workflow_patterns)
dev.off()
png(file = "./output/plot/Workflow_patterns.png", width = 7, height = 2.5, units = "in", res = 600)
print(workflow_patterns)
dev.off()


# 2. for one pattern with medium aggregation, 3 sample examples with same sample size and the 3 sampling methods
# draw the samples
set.seed(20211004)
mysample_random <- sample_random(pop = pattern1, nx = 20, ny = 20, sample_size = 16)
mysample_syst <- sample_systematic_fixed_size(pop = pattern1, nx = 20, ny = 20, sample_size = 16)
mysample_SBS <- sample_BAS(pop = pattern1, nx = 20, ny = 20, sample_size = 16)

# make the plots
plot_random_1 <- plot_sample(pop = pattern1, nx = 20, ny = 20, sample_result = mysample_random, show.grid = FALSE, quadrat.color = "#377eb8")
plot_SYS_1 <- plot_sample(pop = pattern1, nx = 20, ny = 20, sample_result = mysample_syst, show.grid = FALSE, quadrat.color = "#e41a1c")
plot_SBS_1 <- plot_sample(pop = pattern1, nx = 20, ny = 20, sample_result = mysample_SBS, show.grid = FALSE, quadrat.color = "#4daf4a")


# arrange the plots on a single line
workflow_samples <- ggarrange(plot_random_1, plot_SBS_1, plot_SYS_1,
                              ncol = 3, nrow = 1)

# save the plot
pdf(file = "./output/plot/Workflow_samples.pdf", width = 6, height = 2.5)
print(workflow_samples)
dev.off()
png(file = "./output/plot/Workflow_samples.png", width = 6, height = 2.5, units = "in", res = 600)
print(workflow_samples)
dev.off()


# 3. violin plot for this pattern with the 3 sampling methods
# draw the samples for the violin plot on a bigger pattern to have smoother violin plots

set.seed(20210921)
pattern1_BIG <- rMatClust(kappa = 0.02236068, scale = 2.6, mu = 223.6068, win = owin(c(0, 100), c(0, 100)))
pattern1_BIG_true_density <- intensity(pattern1_BIG)

nsim = 500
sample_size = 196
sample_names <- paste0("sample", c(1:nsim))
df_SRS_1 <- data.frame(matrix(NA, nrow = sample_size, ncol = nsim))
df_SYS_1 <- data.frame(matrix(NA, nrow = sample_size, ncol = nsim))
df_SBS_1 <- data.frame(matrix(NA, nrow = sample_size, ncol = nsim))
colnames(df_SRS_1) <- sample_names
colnames(df_SYS_1) <- sample_names
colnames(df_SBS_1) <- sample_names

for (i in 1:nsim) {
  df_SRS_1[, i] <- sample_random(pop = pattern1_BIG, nx = 100, ny = 100, sample_size = sample_size)
  df_SYS_1[, i] <- sample_systematic_fixed_size(pop = pattern1_BIG, nx = 100, ny = 100, sample_size = sample_size)
  df_SBS_1[, i] <- sample_BAS(pop = pattern1_BIG, nx = 100, ny = 100, sample_size = sample_size)
}

# make the recap tables for each pattern in long format
pattern1_BIG_df <- data.frame(Random = colMeans(df_SRS_1),
                              Systematic = colMeans(df_SYS_1),
                              Spatially_Balanced = colMeans(df_SBS_1))
pattern1_BIG_df <- pattern1_BIG_df %>% pivot_longer(cols = c(1:3), names_to = "sampling_method", values_to = "sample_mean")


# make the violin plot
violin_palette <- c("#377eb8", "#4daf4a", "#e41a1c")
#                   random     SBS       Systematic

# force the order of the factor levels so they will be ploted in this order
pattern1_BIG_df$sampling_method <- factor(pattern1_BIG_df$sampling_method,
                                          levels = c("Random", "Spatially_Balanced", "Systematic"), ordered = TRUE)

# rename the level "Spatially_Balanced" to replace the underscore by a space
levels(pattern1_BIG_df$sampling_method)[levels(pattern1_BIG_df$sampling_method) == "Spatially_Balanced"] <- "Spatially Balanced"

# make the plot
violin_pattern1_BIG <- ggplot(pattern1_BIG_df,
                              aes(x = sampling_method, y = sample_mean, fill = sampling_method)) + 
  geom_violin() +
  geom_hline(yintercept = pattern1_BIG_true_density, linetype = "dashed", size = 0.5) +
  xlab("") +
  ylab("Sample mean") +
  labs(fill = "Sampling method") +
  theme_light() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 12), plot.subtitle = element_text(size = 11),
        legend.title = element_text(size = 11),  legend.text = element_text(size = 11),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.y.left = element_text(size = 15),
        panel.spacing = unit(0.8, "lines"),
        strip.text.y = element_text(size = 8, colour = "black"),
        strip.background = element_rect(size = 2, colour = "grey88", fill = "grey88")) +
  scale_fill_manual(values = violin_palette)


# save the plot
pdf(file = "./output/plot/Workflow_violin.pdf", width = 4, height = 3)
print(violin_pattern1_BIG)
dev.off()
png(file = "./output/plot/Workflow_violin.png", width = 4, height = 3, units = "in", res = 600)
print(violin_pattern1_BIG)
dev.off()


# 4. plot one variance ratio curve
simul_result_example <- read.csv2(file = "./data/tidy/var_ratio_example_curve_for_Figure_1.csv")

# convert sample_size to factor
simul_result_example$sample_size <- factor(simul_result_example$sample_size, levels = c("9", "15", "25", "49", "100", "150", "196", "300", "400"))


var_ratio_curve <- ggplot(subset(simul_result_example, simul_result_example$theoretical_intensity == 1 & simul_result_example$sample_size == 400),
                          aes(x = dispersion_index, y = var_ratio_syst)) +
  geom_line(size = 0.8, color = "black") + 
  scale_x_continuous(limits = c(-0.1, 76), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.6, 1.08), expand = expansion(mult = c(0, 0))) +
  ggtitle("Var(Systematic) / Var(Random)") +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  theme_light() + 
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        plot.title = element_blank(), plot.subtitle = element_text(size = 10),
        legend.title = element_text(size = 11),  legend.text = element_text(size = 11),
        axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12),
        panel.spacing = unit(0.8, "lines"),
        strip.text.y = element_text(size = 8, colour = "black"),
        strip.background = element_rect(size = 2, colour = "grey88", fill = "grey88")) +
  guides(color = guide_legend(nrow = 1)) +
  labs(color = "Sample size") +
  xlab("Dispersion index") +
  ylab("Variance ratio")

# save the plot
pdf(file = "./output/plot/Workflow_var_ratio_curve.pdf", width = 6.5, height = 1.8)
print(var_ratio_curve)
dev.off()
png(file = "./output/plot/Workflow_var_ratio_curve.png", width = 6.5, height = 1.8, units = "in", res = 600)
print(var_ratio_curve)
dev.off()



##### FIGURE 2 : SYS and SBS variance ratio as a function of dispersion index #####

palette <- rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')) # palette for 9 sample size

figure2_SYS <- ggplot(subset(simul_result, simul_result$theoretical_intensity %in% c(1, 10, 20)),
                      aes(x = dispersion_index, y = var_ratio_syst, color = sample_size)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  geom_line(size = 0.6) + 
  scale_x_continuous(limits = c(-0.1, 76), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.15, 1.2), expand = expansion(mult = c(0, 0))) +
  ggtitle("Var(Systematic) / Var(Random)") +
  facet_grid(theoretical_intensity ~ .) + 
  theme_light() + 
  theme(legend.position = "bottom", axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.title = element_text(size = 11), plot.subtitle = element_text(size = 10),
        legend.title = element_text(size = 11),  legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        strip.text.y = element_text(size = 10, colour = "black"),
        strip.background = element_rect(size = 2, colour = "grey88", fill = "grey88")) +
  guides(color = guide_legend(nrow = 1)) +
  labs(color = "Sample size") +
  scale_color_manual(values = palette)

figure2_SBS <- ggplot(subset(simul_result, simul_result$theoretical_intensity %in% c(1, 10, 20)),
                      aes(x = dispersion_index, y = var_ratio_BAS, color = sample_size)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  geom_line(size = 0.6) + 
  scale_x_continuous(limits = c(-0.1, 76), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.15, 1.2), expand = expansion(mult = c(0, 0))) +
  ggtitle("Var(Spatially Balanced) / Var(Random)") +
  facet_grid(theoretical_intensity ~ .) + 
  theme_light() + 
  theme(legend.position = "bottom", axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.title = element_text(size = 11), plot.subtitle = element_text(size = 10),
        legend.title = element_text(size = 11),  legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        strip.text.y = element_text(size = 10, colour = "black"),
        strip.background = element_rect(size = 2, colour = "grey88", fill = "grey88")) +
  guides(color = guide_legend(nrow = 1)) +
  labs(color = "Sample size") +
  scale_color_manual(values = palette)

# arrange the plots on a single page
figure2 <- ggarrange(figure2_SYS, figure2_SBS,
                     ncol = 1, nrow = 2, common.legend = TRUE, legend = "top")

figure2 <- annotate_figure(figure2,
                           bottom = text_grob("Dispersion index", size = 11),
                           left = text_grob("Variance ratio", size = 11, rot = 90))

# save the plot
pdf(file = "./output/plot/Figure2_SYS_SBS_variance_ratio_disp_index.pdf", width = 6.4, height = 7.6)
print(figure2)
dev.off()
png(file = "./output/plot/Figure2_SYS_SBS_variance_ratio_disp_index.png", width = 6.4, height = 7.6, units = "in", res = 600)
print(figure2)
dev.off()



##### FIGURE 3 : SYS variance ratio as a function of cluster diameter #####

dist_df_syst <- data.frame(
  sample_size = c(     9,     15,     25,     49,    100,   150,   196,   300,   400),
  mean_dist =   c(40.232, 30.733, 24.142, 17.243, 12.071, 9.561, 8.621, 6.666, 6.035)) # accounting for direct neighbours and diagonal neighbours

# convert sample_size to factor with right order
dist_df_syst$sample_size <- factor(dist_df_syst$sample_size,
                                   levels = c("9", "15", "25", "49", "100", "150", "196", "300", "400"))

# restrict plots to levels given below
wanted_intensities = c(1, 5, 10, 15, 20, 25, 30)
wanted_sample_sizes = c(49, 100, 196, 300, 400)

palette2 <- c('#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026') # palette for 7 densities

# make the plot
figure3 <- ggplot(subset(simul_result, simul_result$theoretical_intensity %in% wanted_intensities &
                           simul_result$sample_size %in% wanted_sample_sizes),
                  aes(x = cluster_diam, y = var_ratio_syst, color = as.factor(theoretical_intensity))) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  geom_line(size = 0.6) + 
  scale_x_continuous(limits = c(-0.1, 26), expand = c(0, 0)) +
  ggtitle("Var(Systematic) / Var(Random)") +
  facet_grid(sample_size ~ .) +
  geom_vline(data = subset(dist_df_syst, dist_df_syst$sample_size %in% wanted_sample_sizes),
             aes(xintercept = mean_dist), size = 0.5) +
  theme_light() +
  theme(legend.position = "bottom", axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 11), plot.subtitle = element_text(size = 10),
        legend.title = element_text(size = 11),  legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        strip.text.y = element_text(size = 10, colour = "black"),
        strip.background = element_rect(size = 2, colour = "grey88", fill = "grey88"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, -6, -6, -6)) +
  guides(color = guide_legend(nrow = 1)) +
  labs(color = "Population density (individuals/cell)") +
  xlab("Cluster diameter") +
  ylab("Variance ratio") +
  scale_color_manual(values = palette2)


pdf(file = "./output/plot/Figure3_SYS_variance_ratio_cluster_diameter.pdf", width = 6.4, height = 7)
print(figure3)
dev.off()
png(file = "./output/plot/Figure3_SYS_variance_ratio_cluster_diameter.png", width = 6.4, height = 7, units = "in", res = 600)
print(figure3)
dev.off()



##### FIGURE 4 : results of the 3 plant populations #####
### make the heatmaps
# load data
matrix_SANMIN_long <- read.csv2(file = "./data/tidy/matrix_SANMIN_long.csv")
matrix_LIMGIR_long <- read.csv2(file = "./data/tidy/matrix_LIMGIR_long.csv")
matrix_BELSYL_long <- read.csv2(file = "./data/tidy/matrix_BELSYL_long.csv")

# viridis color palette
viridis_pal <- c("#440154FF", "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")

# make the plots
heatmap_SANMIN <- ggplot(matrix_SANMIN_long, aes(x = x, y = -y, fill = count)) + 
  geom_tile() +
  theme_bw() +
  scale_fill_gradientn(colors = viridis_pal) +
  coord_fixed() +
  theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(),
        legend.position = "bottom", panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.background = element_blank(), panel.background = element_blank(),
        legend.title = element_text(size = 17, vjust = 1, margin = margin(t = 7.5)),
        legend.text = element_text(size = 16),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"), legend.box.margin = margin(0, 0, 0, 0, unit = "cm"),
        plot.margin = margin(0, 0, 0, 0, unit = "cm"),
        legend.box.spacing = unit(-0.3, "cm"),
        legend.key.size = unit(1.1, "cm")) +
  labs(fill = "Individuals")

heatmap_LIMGIR <- ggplot(matrix_LIMGIR_long, aes(x = x, y = -y, fill = count)) + 
  geom_tile() +
  theme_bw() +
  scale_fill_gradientn(colors = viridis_pal) +
  coord_fixed() +
  theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(),
        legend.position = "bottom", panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.background = element_blank(), panel.background = element_blank(),
        legend.title = element_text(size = 17, vjust = 1, margin = margin(t = 7.5)),
        legend.text = element_text(size = 16),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"), legend.box.margin = margin(0, 0, 0, 0, unit = "cm"),
        plot.margin = margin(0, 0, 0, 0, unit = "cm"),
        legend.box.spacing = unit(-0.3, "cm"),
        legend.key.size = unit(1.1, "cm")) +
  labs(fill = "Individuals")

heatmap_BELSYL <- ggplot(matrix_BELSYL_long, aes(x = x, y = -y, fill = count)) + 
  geom_tile() +
  scale_fill_gradientn(colors = viridis_pal) +
  coord_fixed() +
  theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(),
        legend.position = "bottom", panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.background = element_blank(), panel.background = element_blank(),
        legend.title = element_text(size = 17, vjust = 1, margin = margin(t = 7.5)),
        legend.text = element_text(size = 16),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"), legend.box.margin = margin(0, 0, 0, 0, unit = "cm"),
        plot.margin = margin(0, 0, 0, 0, unit = "cm"),
        legend.box.spacing = unit(-0.3, "cm"),
        legend.key.size = unit(1.1, "cm")) +
  labs(fill = "Individuals")


# save the heatmaps
png(file = "./output/plot/Heatmap_SANMIN.png", width = 5, height = 5.1, units = "in", res = 600)
print(heatmap_SANMIN)
dev.off()
png(file = "./output/plot/Heatmap_LIMGIR.png", width = 5, height = 5.1, units = "in", res = 600)
print(heatmap_LIMGIR)
dev.off()
png(file = "./output/plot/Heatmap_BELSYL.png", width = 5, height = 5.1, units = "in", res = 600)
print(heatmap_BELSYL)
dev.off()


### make the variance ratio curves
# load data
tibble_natural_pop_var_ratio_long <- read.csv2(file = "./data/tidy/variance_ratio_table_natural_populations_60_sample_sizes.csv")
tibble_natural_pop_var_ratio_long <- tibble(tibble_natural_pop_var_ratio_long)

# get mean cluster diameter
cluster_diameter_SANMIN <- read.csv2(file = "./data/tidy/cluster_diameters_sanguisorba_minor.csv", header = TRUE)
cluster_diameter_LIMGIR <- read.csv2(file = "./data/tidy/cluster_diameters_limonium_girardianum.csv", header = TRUE)
cluster_diameter_BELSYL <- read.csv2(file = "./data/tidy/cluster_diameters_bellis_sylvestris.csv", header = TRUE)
mean_diameter_SANMIN <- mean(as.matrix(cluster_diameter_SANMIN[, 2:6]), na.rm = TRUE) # 158.12
mean_diameter_LIMGIR <- mean(as.matrix(cluster_diameter_LIMGIR[, 2:6]), na.rm = TRUE) # 104.04
mean_diameter_BELSYL <- mean(as.matrix(cluster_diameter_BELSYL[, 2:6]), na.rm = TRUE) #  50.45


# get mean inter-unit distance for SYS
sample_sizes_vect <- c(9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256, 289, 324, 361, 400, 441, 484, 529, 576, 625)

dist_direct_neighbours <- 100/sqrt(sample_sizes_vect)
dist_diagonal_neighbours <- sqrt((dist_direct_neighbours^2)*2)
mean_dist <- (dist_direct_neighbours + dist_diagonal_neighbours)/2

dist_df_SYS <- data.frame(
  sample_size = sample_sizes_vect,
  mean_dist = mean_dist)


# get mean inter-unit distance for SBS
df_SBS_dist_neighbours <- read.csv2(file = "./data/tidy/BAS_mean_distance_nearest_units_14_neighbours.csv")
df_SBS_dist_neighbours <- cbind(df_SBS_dist_neighbours, mean_dist = rowMeans(df_SBS_dist_neighbours[, 2:13]))


# # mean cluster diameter expressed in number of cells
# mean_diameter_SANMIN/20 # 7.90625
# mean_diameter_LIMGIR/20 # 5.202128
# mean_diameter_BELSYL/20 # 2.522727 --> clusters are too small for the sample sizes we simulated

# make a dataframe containing the sample size at which the mean cluster diameter is equal to the mean inter-unit distance (taken by linear approximation between the two sample sizes surrounding the mean cluster diameter)
df_SYS_vline_position <- data.frame(species_name = c("Bellis sylvestris", "Limonium girardianum", "Sanguisorba minor"),
                                    sample_size =  c(                 NA,                    539,               235))
df_SBS_vline_position <- data.frame(species_name = c("Bellis sylvestris", "Limonium girardianum", "Sanguisorba minor"),
                                    sample_size =  c(                 NA,                     NA,               367))

# join species complete name to the data frame
df_species_name <- data.frame(species = c("BELSYL", "LIMGIR", "SANMIN"),
                              species_name = c("Bellis sylvestris", "Limonium girardianum", "Sanguisorba minor"))
tibble_natural_pop_var_ratio_long <- left_join(tibble_natural_pop_var_ratio_long, df_species_name, by = "species")


## make one plot for every species
Fieldwork_var_ratio_BELSYL <- ggplot(subset(tibble_natural_pop_var_ratio_long, tibble_natural_pop_var_ratio_long$species == "BELSYL"),
                                     aes(x = sample_size, y = var_ratio, color = sampling_method)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) + 
  geom_line(size = 1) + 
  scale_color_manual(values = c("#4daf4a", "#e41a1c")) +
  scale_x_continuous(limits = c(-0.1, 625), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.21, 1.5), expand = c(0, 0)) +
  theme_light() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10),
        legend.title = element_text(size = 11),  legend.text = element_text(size = 11),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 11),
        panel.spacing = unit(2, "cm"),
        strip.text.y = element_text(size = 11, colour = "black"),
        strip.background = element_rect(size = 2, colour = "grey88", fill = "grey88"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, -6, -6, -6)) +
  guides(color = guide_legend(nrow = 1)) +
  labs(color = "Sampling method") +
  xlab("Sample size") +
  ylab("Variance ratio")


Fieldwork_var_ratio_LIMGIR <- ggplot(subset(tibble_natural_pop_var_ratio_long, tibble_natural_pop_var_ratio_long$species == "LIMGIR"),
                                     aes(x = sample_size, y = var_ratio, color = sampling_method)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) + 
  geom_line(size = 1) + 
  scale_color_manual(values = c("#4daf4a", "#e41a1c")) +
  scale_x_continuous(limits = c(-0.1, 625), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.21, 1.5), expand = c(0, 0)) +
  theme_light() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10),
        legend.title = element_text(size = 11),  legend.text = element_text(size = 11),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 11),
        panel.spacing = unit(2, "cm"),
        strip.text.y = element_text(size = 11, colour = "black"),
        strip.background = element_rect(size = 2, colour = "grey88", fill = "grey88"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, -6, -6, -6)) +
  geom_vline(data = subset(df_SYS_vline_position, df_SYS_vline_position$species_name == "Limonium girardianum"), aes(xintercept = sample_size), size = 0.5, color = "#4daf4a") +
  geom_vline(data = subset(df_SBS_vline_position, df_SBS_vline_position$species_name == "Limonium girardianum"), aes(xintercept = sample_size), size = 0.5, color = "#e41a1c") +
  guides(color = guide_legend(nrow = 1)) +
  labs(color = "Sampling method") +
  xlab("Sample size") +
  ylab("Variance ratio")


Fieldwork_var_ratio_SANMIN <- ggplot(subset(tibble_natural_pop_var_ratio_long, tibble_natural_pop_var_ratio_long$species == "SANMIN"),
                                     aes(x = sample_size, y = var_ratio, color = sampling_method)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) + 
  geom_line(size = 1) + 
  scale_color_manual(values = c("#4daf4a", "#e41a1c")) +
  scale_x_continuous(limits = c(-0.1, 625), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.21, 1.5), expand = c(0, 0)) +
  theme_light() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10),
        legend.title = element_text(size = 11),  legend.text = element_text(size = 11),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 11),
        panel.spacing = unit(2, "cm"),
        strip.text.y = element_text(size = 11, colour = "black"),
        strip.background = element_rect(size = 2, colour = "grey88", fill = "grey88"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, -6, -6, -6)) +
  geom_vline(data = subset(df_SYS_vline_position, df_SYS_vline_position$species_name == "Sanguisorba minor"), aes(xintercept = sample_size), size = 0.5, color = "#4daf4a") +
  geom_vline(data = subset(df_SBS_vline_position, df_SBS_vline_position$species_name == "Sanguisorba minor"), aes(xintercept = sample_size), size = 0.5, color = "#e41a1c") +
  guides(color = guide_legend(nrow = 1)) +
  labs(color = "Sampling method") +
  xlab("Sample size") +
  ylab("Variance ratio")


# save the plots
png(file = "./output/plot/Fieldwork_var_ratio_BELSYL.png", width = 5.1, height = 2.5, units = "in", res = 600)
print(Fieldwork_var_ratio_BELSYL)
dev.off()
png(file = "./output/plot/Fieldwork_var_ratio_LIMGIR.png", width = 5.1, height = 2.5, units = "in", res = 600)
print(Fieldwork_var_ratio_LIMGIR)
dev.off()
png(file = "./output/plot/Fieldwork_var_ratio_SANMIN.png", width = 5.1, height = 2.5, units = "in", res = 600)
print(Fieldwork_var_ratio_SANMIN)
dev.off()



##### FIGURE 5 : mechanism illustration #####
# simulate the pattern to be plotted
my_window <- owin(c(0, 50), c(0, 50))

set.seed(202100)
pattern1 <- rMatClust(kappa = 0.003, scale = 7, mu = 300, win = my_window)

# replace pattern1 by an empty matrix of dim 20*20
pattern1 <- as.data.frame.matrix(matrix(data = 0, nrow = 20, ncol = 20))

# set aspect for all the sample illustrations
my_pointsize = 0.0025
my_pointalpha = 1
my_quadratalpha = 1

# draw the samples and make the illustration plots
# sample size n = 9
plot_random_n9_1 <- plot_sample(pop = pattern1, nx = 20, ny = 20, 
                                sample_result = sample_random(pop = pattern1, nx = 20, ny = 20, sample_size = 9),
                                show.grid = FALSE, quadrat.color = "#377eb8", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_random_n9_2 <- plot_sample(pop = pattern1, nx = 20, ny = 20, 
                                sample_result = sample_random(pop = pattern1, nx = 20, ny = 20, sample_size = 9),
                                show.grid = FALSE, quadrat.color = "#377eb8", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_random_n9_3 <- plot_sample(pop = pattern1, nx = 20, ny = 20, 
                                sample_result = sample_random(pop = pattern1, nx = 20, ny = 20, sample_size = 9),
                                show.grid = FALSE, quadrat.color = "#377eb8", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_random_n9_4 <- plot_sample(pop = pattern1, nx = 20, ny = 20, 
                                sample_result = sample_random(pop = pattern1, nx = 20, ny = 20, sample_size = 9),
                                show.grid = FALSE, quadrat.color = "#377eb8", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)

plot_SYS_n9_1 <- plot_sample(pop = pattern1, nx = 20, ny = 20,
                             sample_result = sample_systematic_fixed_size(pop = pattern1, nx = 20, ny = 20, sample_size = 9),
                             show.grid = FALSE, quadrat.color = "#e41a1c", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_SYS_n9_2 <- plot_sample(pop = pattern1, nx = 20, ny = 20,
                             sample_result = sample_systematic_fixed_size(pop = pattern1, nx = 20, ny = 20, sample_size = 9),
                             show.grid = FALSE, quadrat.color = "#e41a1c", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_SYS_n9_3 <- plot_sample(pop = pattern1, nx = 20, ny = 20,
                             sample_result = sample_systematic_fixed_size(pop = pattern1, nx = 20, ny = 20, sample_size = 9),
                             show.grid = FALSE, quadrat.color = "#e41a1c", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_SYS_n9_4 <- plot_sample(pop = pattern1, nx = 20, ny = 20,
                             sample_result = sample_systematic_fixed_size(pop = pattern1, nx = 20, ny = 20, sample_size = 9),
                             show.grid = FALSE, quadrat.color = "#e41a1c", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)

# sample size n = 16
plot_random_n16_1 <- plot_sample(pop = pattern1, nx = 20, ny = 20, 
                                 sample_result = sample_random(pop = pattern1, nx = 20, ny = 20, sample_size = 16),
                                 show.grid = FALSE, quadrat.color = "#377eb8", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_random_n16_2 <- plot_sample(pop = pattern1, nx = 20, ny = 20, 
                                 sample_result = sample_random(pop = pattern1, nx = 20, ny = 20, sample_size = 16),
                                 show.grid = FALSE, quadrat.color = "#377eb8", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_random_n16_3 <- plot_sample(pop = pattern1, nx = 20, ny = 20, 
                                 sample_result = sample_random(pop = pattern1, nx = 20, ny = 20, sample_size = 16),
                                 show.grid = FALSE, quadrat.color = "#377eb8", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_random_n16_4 <- plot_sample(pop = pattern1, nx = 20, ny = 20, 
                                 sample_result = sample_random(pop = pattern1, nx = 20, ny = 20, sample_size = 16),
                                 show.grid = FALSE, quadrat.color = "#377eb8", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)

plot_SYS_n16_1 <- plot_sample(pop = pattern1, nx = 20, ny = 20,
                              sample_result = sample_systematic_fixed_size(pop = pattern1, nx = 20, ny = 20, sample_size = 16),
                              show.grid = FALSE, quadrat.color = "#e41a1c", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_SYS_n16_2 <- plot_sample(pop = pattern1, nx = 20, ny = 20,
                              sample_result = sample_systematic_fixed_size(pop = pattern1, nx = 20, ny = 20, sample_size = 16),
                              show.grid = FALSE, quadrat.color = "#e41a1c", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_SYS_n16_3 <- plot_sample(pop = pattern1, nx = 20, ny = 20,
                              sample_result = sample_systematic_fixed_size(pop = pattern1, nx = 20, ny = 20, sample_size = 16),
                              show.grid = FALSE, quadrat.color = "#e41a1c", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)
plot_SYS_n16_4 <- plot_sample(pop = pattern1, nx = 20, ny = 20,
                              sample_result = sample_systematic_fixed_size(pop = pattern1, nx = 20, ny = 20, sample_size = 16),
                              show.grid = FALSE, quadrat.color = "#e41a1c", pointsize = my_pointsize, pointalpha = my_pointalpha, quadratalpha = my_quadratalpha)

### function to make the circles
circleFun <- function(center = c(0, 0), diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0, 2*pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# data for the circles
dat1 <- circleFun(center = c(9, 10), diameter = 5, npoints = 1000)
dat2 <- circleFun(center = c(14, 15), diameter = 5, npoints = 1000)
dat3 <- circleFun(center = c(14.5, 6), diameter = 5, npoints = 1000)

vect_circles <- c(geom_polygon(data = dat1, aes(x,y), fill = "#4daf4a", colour = "black", alpha = 0.9, size = 0.1),
                  geom_polygon(data = dat2, aes(x,y), fill = "#4daf4a", colour = "black", alpha = 0.9, size = 0.1),
                  geom_polygon(data = dat3, aes(x,y), fill = "#4daf4a", colour = "black", alpha = 0.9, size = 0.1))

# add the circles to the plots with the samples
plot_random_n9_1$layers <- c(vect_circles, plot_random_n9_1$layers)
plot_SYS_n9_1$layers <- c(vect_circles, plot_SYS_n9_1$layers)
plot_random_n16_1$layers <- c(vect_circles, plot_random_n16_1$layers)
plot_SYS_n16_1$layers <- c(vect_circles, plot_SYS_n16_1$layers)
plot_random_n9_2$layers <- c(vect_circles, plot_random_n9_2$layers)
plot_SYS_n9_2$layers <- c(vect_circles, plot_SYS_n9_2$layers)
plot_random_n16_2$layers <- c(vect_circles, plot_random_n16_2$layers)
plot_SYS_n16_2$layers <- c(vect_circles, plot_SYS_n16_2$layers)
plot_random_n9_3$layers <- c(vect_circles, plot_random_n9_3$layers)
plot_SYS_n9_3$layers <- c(vect_circles, plot_SYS_n9_3$layers)
plot_random_n16_3$layers <- c(vect_circles, plot_random_n16_3$layers)
plot_SYS_n16_3$layers <- c(vect_circles, plot_SYS_n16_3$layers)
plot_random_n9_4$layers <- c(vect_circles, plot_random_n9_4$layers)
plot_SYS_n9_4$layers <- c(vect_circles, plot_SYS_n9_4$layers)
plot_random_n16_4$layers <- c(vect_circles, plot_random_n16_4$layers)
plot_SYS_n16_4$layers <- c(vect_circles, plot_SYS_n16_4$layers)


figure5 <- ggarrange(plot_random_n9_1, plot_SYS_n9_1, plot_random_n16_1, plot_SYS_n16_1, # first row
                     plot_random_n9_2, plot_SYS_n9_2, plot_random_n16_2, plot_SYS_n16_2, # 2nd row
                     plot_random_n9_3, plot_SYS_n9_3, plot_random_n16_3, plot_SYS_n16_3, # 3rd row
                     ncol = 4, nrow = 3, vjust = -0.1, hjust = -1.6) +
  theme(plot.margin = margin(1.8, 0.1, 0.1, 0.1, "cm"))

# save the plot
pdf(file = "./output/plot/Figure5_mecanism_illustration_circles.pdf", width = 6.4, height = 5.6)
print(figure5)
dev.off()
png(file = "./output/plot/Figure5_mecanism_illustration_circles.png", width = 6.4, height = 5.6, units = "in", res = 600)
print(figure5)
dev.off()


