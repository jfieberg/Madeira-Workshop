#' This file can be used to create html reports from
#' each of the files in the repository and move them 
#' to the output directory.
#' 
#' Use the here package to locate each of the files 
#' using a relative path.
library(here)
library(rmarkdown)

#' Source directory for these files 
#' 
render(here("Resource-and-Step-Selection-Analyses/RCode", "amt_demo_iSSF.R"), output_dir = here("Resource-and-Step-Selection-Analyses/","Output"))
render(here("Resource-and-Step-Selection-Analyses/RCode", "Otters_SSF.R"), output_dir = here("Resource-and-Step-Selection-Analyses/","Output"))
render(here("Resource-and-Step-Selection-Analyses/RCode", "OttersTwoStep.R"), output_dir = here("Resource-and-Step-Selection-Analyses/","Output"))
render(here("Resource-and-Step-Selection-Analyses/RCode", "MultipleAnimalsRSFs.R"), output_dir = here("Resource-and-Step-Selection-Analyses/","Output"))
render(here("Resource-and-Step-Selection-Analyses/RCode", "MultipleAnimalsSSFs.R"), output_dir = here("Resource-and-Step-Selection-Analyses/","Output"))
