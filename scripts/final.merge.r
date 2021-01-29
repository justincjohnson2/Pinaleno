# ------------------------------------------------------------------------------
# Description:  Code to merge all data from sampling in Pinalenos 2020 into a 
#               full dataset.
# Date:         1/28/21
# Author:       Justin Johnson
# ------------------------------------------------------------------------------
# Install and load packages
pkgs <- c('tidyverse')

check <- sapply(pkgs,require,warn.conflicts = TRUE,character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing,require,warn.conflicts = TRUE,character.only = TRUE)
}

#Inputs
bdfilename <- "raw_data/bulk_density_Pinalenos_2020.csv"
mdifilename <- "processed_data/MDI_Pinalenos_2020.csv"
soilfilename <- "raw_data/soiltexture_carbon_nitrogen_Pinalenos_2020.csv"
wdptfilename <- "raw_data/WDPT_cover_Pinalenos_2020.csv"
outfilename <- "processed_data/Full_Pinalenos_2020.csv"

# ------------------------------------------------------------------------------
# Read in bulk density data
bd <- read_csv(bdfilename, col_types = "c----n----------")

# ------------------------------------------------------------------------------
# Read in MDI data
mdi <- read_csv(mdifilename, col_types = "c----nnn--------nnnn------------")

# ------------------------------------------------------------------------------
# Read in soil texture, total organic carbon, and total nitrogen data
soil <- read_csv(soilfilename, col_types = "c----cnnnnn--c") %>%
  rename("notes_soil"= Notes)

# ------------------------------------------------------------------------------
# Read in WDPT, soil moisture, litter depth, and cover data
wdpt <- read_csv(wdptfilename, col_types = "c----nnnnnnnnnnnnllnnccccccc") %>%
  mutate(LF1 = replace(LF1, is.na(LF1), ""),
         LF2 = replace(LF2, is.na(LF2), ""),
         LF3 = replace(LF3, is.na(LF3), ""),
         overstory1 = replace(overstory1, is.na(overstory1), ""),
         overstory2 = replace(overstory2, is.na(overstory2), ""))

# ------------------------------------------------------------------------------
# Joins data sets
final <- full_join(wdpt, mdi, by="location_ID") %>%
  full_join(bd, by="location_ID") %>%
  full_join(soil, by="location_ID") %>%
  separate(location_ID, c("plot_ID", "x.m", "y.m"), sep="_", remove = FALSE) %>%
  mutate(plot_ID = as.numeric(plot_ID),
         burn_severity = case_when(plot_ID < 4 ~ "high",
                                   plot_ID > 3 & plot_ID < 7 ~ "moderate",
                                   plot_ID > 6 & plot_ID < 10 ~ "low",
                                   plot_ID > 9 ~"unburned")) %>%
  select(location_ID, plot_ID, x.m, y.m, burn_severity, everything())

# ------------------------------------------------------------------------------
# Write csv of cleaned data.
write_csv(final,outfilename)