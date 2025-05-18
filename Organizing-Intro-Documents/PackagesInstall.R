#' Run this script to make sure you have all the necessary packages for 
#' this week.  Then, render the Test.qmd file to demonstrate that everything
#' is working and upload the resulting .html file to the following google drive
#' folder: https://drive.google.com/drive/folders/19PJ11vSP2bQ8ObEj12fP0y_CcWcI139I

#  
# The code below will look for the various packages that we will use during 
# the workshop and make sure they are installed on your computer. 

# Install (or update) all other packages
pkg.required <- c("amt",
                  "sf",
                  "tidyverse",
                  "lubridate",
                  "TwoStepCLogit",
                  "glmmTMB",
                  "tictoc",
                  "broom",
                  "here",
                  "terra",
                  "moveHMM",
                  "boot",
                  "survival",
                  "ctmm",
                  "remotes",
                  "viridis",
                  "patchwork",
                  "readr")

pkg.available <- installed.packages()[, "Package"]

if (!all(pkg.required %in% pkg.available)) {
  pkg.missing <- pkg.required[!pkg.required %in% pkg.available]
  install.packages(pkg.missing, dependencies = TRUE)
} 

# install 'develop' version of momentuHMM
install.packages("momentuHMM", repos = "https://bmcclintock.r-universe.dev", dependencies = TRUE)
remotes::install_github("ecoisilva/movedesign")

pkg.required <- append(pkg.required, c("momentuHMM", "movedesign"))

sapply(pkg.required, require, character.only = TRUE)

