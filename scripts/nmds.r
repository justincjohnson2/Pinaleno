# ------------------------------------------------------------------------------
# Description:  Script to run non-metric multidimensional scaling on plot-scale
#               data.
# Date:         1/29/21
# Author:       Justin Johnson
# ------------------------------------------------------------------------------
# Install and load packages
pkgs <- c('tidyverse', 'EnvStats')

check <- sapply(pkgs,require,warn.conflicts = TRUE,character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing,require,warn.conflicts = TRUE,character.only = TRUE)
}

#Inputs
filename <- "processed_data/Full_Pinalenos_2020.csv"

# ------------------------------------------------------------------------------
# Read in full dataset
full <- read_csv(filename, col_types = "cfnnfnnnnnnnnnnnnllnncccccc-nn-----n-nnnnn-")

# ------------------------------------------------------------------------------
# Get average value for each variable on each plot. Means are summary statistic
# unless otherwise mentioned.

# Calculate median WDPT for each depth at each point.
full1 <- full %>%
  rowwise() %>%
  mutate(WD_0cm_median.sec = median(WD_0cm_1.sec, 
                                    WD_0cm_2.sec, 
                                    WD_0cm_3.sec,
                                    na.rm = TRUE),
         WD_1cm_median.sec = median(WD_1cm_1.sec, 
                                    WD_1cm_2.sec, 
                                    WD_1cm_3.sec,
                                    na.rm = TRUE),
         WD_2cm_median.sec = median(WD_2cm_1.sec, 
                                    WD_2cm_2.sec, 
                                    WD_2cm_3.sec,
                                    na.rm = TRUE),
         WD_3cm_median.sec = median(WD_3cm_1.sec, 
                                    WD_3cm_2.sec, 
                                    WD_3cm_3.sec,
                                    na.rm = TRUE))
full$WD_4cm.hydrophobic
# Calculate means/proportions for most variables
full_means <- full1 %>%
  group_by(plot_ID) %>%
  summarize(burn_severity = sample(burn_severity, 1),
            WD_0cm_mean.sec = mean(WD_0cm_median.sec, na.rm = TRUE),
            WD_1cm_mean.sec = mean(WD_1cm_median.sec, na.rm = TRUE),
            WD_2cm_mean.sec = mean(WD_2cm_median.sec, na.rm = TRUE),
            WD_3cm_mean.sec = mean(WD_0cm_median.sec, na.rm = TRUE),
            WD_4cm.prop = mean(WD_4cm.hydrophobic, na.rm = TRUE),
            WD_5cm.prop = mean(WD_5cm.hydrophobic, na.rm = TRUE),
            litter_depth_mean.cm = mean(litter_depth.cm, na.rm = TRUE),
            theta_i_mean.percent = mean(theta_i.percent, na.rm = TRUE),
            S_geomean.mm_hr0.5 = geoMean(S.mm_hr0.5, na.rm = TRUE),
            kfs_geomean.mm_hr = geoMean(kfs.mm_hr, na.rm = TRUE),
            bulk_density_mean.gcm3 = mean(bulk_density.gcm3, na.rm = TRUE),
            sand_mean.percent = mean(sand.percent, na.rm = TRUE),
            silt_mean.percent = mean(silt.percent, na.rm = TRUE),
            clay_mean.percent = mean(clay.percent, na.rm = TRUE),
            total_organic_carbon_mean.percent = mean(total_organic_carbon.percent, na.rm = TRUE),
            total_nitrogen_mean.percent = mean(total_nitrogen.percent, na.rm = TRUE))

# Calculate proportion ground cover for each plot
ground_cover <- full %>%
  drop_na(ground_cover) %>%
  group_by(plot_ID) %>%
  count(ground_cover) %>%
  mutate(ground_prop = (n/sum(n)),
         ground_cover = recode(ground_cover, 
                               B = "grd_bare_soil.prop",
                               F = "grd_forb.prop",
                               L = "grd_litter.prop",
                               G = "grd_graminoid.prop",
                               R = "grd_rock.prop",
                               SD = "grd_standing_dead.prop",
                               WD = "grd_woody_dead.prop",
                               D = "grd_feces.prop",
                               M = "grd_cryptogram.prop",
                               T = "grd_tree.prop")) %>%
  select(-n) %>%
  spread(key = ground_cover, value = ground_prop) %>%
  replace(is.na(.), 0) %>%
  mutate(grd_basal.prop = grd_forb.prop + grd_graminoid.prop + grd_tree.prop)

# Calculate proportion understory canopy cover
canopy <- full %>%
  drop_na(ground_cover) %>% # drop points where cover was not measured
  group_by(plot_ID) %>%
  count(LF1) %>%
  rename(n_LF1 = n,
         cover = LF1)

canopy2 <- full %>%
  drop_na(ground_cover) %>% # drop points where cover was not measured
  group_by(plot_ID) %>%
  count(LF2) %>%
  rename(n_LF2 = n,
         cover = LF2)

total_canopy <- full %>%
  drop_na(ground_cover) %>% # drop points where cover was not measured
  mutate(foliar = case_when(
    LF1 == "G" | LF1 == "F" | LF1 == "T" | LF2 == "G" | LF2 == "F" | LF2 == "T" ~ TRUE,
    TRUE ~ FALSE)) %>%
  group_by(plot_ID) %>%
  summarize(can_total.prop = sum(!is.na(LF1))/n(),
            can_foliar.prop = mean(foliar))

canopy_join <- full_join(canopy, canopy2, by = c("plot_ID", "cover")) %>%
  rowwise() %>%
  mutate(total = sum(n_LF1, n_LF2, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(plot_ID) %>%
  mutate(cover_prop = (total/sum(n_LF1, na.rm = TRUE)),
         cover = recode(cover,
                        F = "can_forb.prop",
                        G = "can_graminoid.prop",
                        L = "can_litter.prop",
                        SD = "can_standing_dead.prop",
                        WD = "can_woody_dead.prop",
                        T = "can_tree.prop")) %>%
  drop_na(cover) %>%
  select(plot_ID, cover, cover_prop) %>%
  spread(key = cover, value = cover_prop) %>%
  replace(is.na(.), 0) %>%
  full_join(total_canopy, by = "plot_ID")

# Calculate proportion overstory cover
