# ------------------------------------------------------------------------------
# Description:  Runs all scripts in Pinaleno repository.
# Date:         1/29/21
# Author:       Justin Johnson
# ------------------------------------------------------------------------------
# Run MDI processing script
source("scripts/MDI_processing.r")

# Run script to merge datasets
source("scripts/final.merge.r")

# Run non-metric multidimensional scaling script
source("scripts/nmds.r")