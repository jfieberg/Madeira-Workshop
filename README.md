# Madeira-Workshop

Presentations, code, and data associated with the workshop titled, Analysing animal movement in the marine environment, held in Madeira, Portugal from May 19-23, 2025.

## Preliminary Schedule

General format:

- We will cover a different topic on each day, starting at 9:00am and ending at 4:30pm.  
- At the end of the day, a select/volunteer group will meet to discuss how what we learned on that day may inform a broader perspective paper addressing: a) questions we want to tackle as marine scientists, b) currently available methods (and their relative accessibility to non-statisticians) for addressing these questions;  c) outstanding challenges that current methods are not fully capable of addressing.

Monday, May 19 (Day 1):

- 9:30 Introduction / registration
- 10:00 Welcome from the director of the institute
- 10:30 Discussion related to aims, objectives, and particularities of the marine environment
- 11:30 Data cleaning, dealing with measurement error
- 12:30 Lunch
- 14:00 Introduction to State-space models (Marie Auger-Méthé)
- 16:00 Wrap-up/discuss initial outline for discussion paper


Tuesday, May 20 (Day 2):

- Discrete time Hidden Markov Models (HMMs) (Brett McClintock)
- Continuous time HMM's
- State-switching habitat selection models in continuous time


Wednesday, May 21 (Day 3):

- Introduction to species distribution/resource selection models (John Fieberg)
- Introduction to step-selection analyses
- Addressing individual variability using 2-step approaches and mixed-effect models
- Model validation using simulations

Thursday, May 22 (Day 4):

- Introduction to ctmm (Inês Silva)
- Home ranges 
- Autocorrelated resource selection analyses
- Mean speed, accounting for measurement error, etc 

Friday, May 23 (Day 5):

- Advanced state-space models (Marie Auger-Méthé)
- Afternoon, outside speakers 
- Wrap up discussion (next steps for the perspective paper)


## Getting set up before the workshop

### R and RStudio

<!--- Latest version of R is 4.5 (released on April 11, 2025), prevoous was 4.4.3 (released on  Feb 28, 2025) --->

Please make sure you have [installed `R` version >4.2.0](https://cran.r-project.org/). We also recommend the latest version of [RStudio Desktop](https://www.rstudio.com/products/rstudio/download/). 

<!---
To build packages from source, you will need additional build tools; see details [here for Windows](https://cran.r-project.org/bin/windows/Rtools/) or [here for macOS](https://mac.r-project.org/tools/). 

*Note* that if you are upgrading to R 4.2 from a previous version on Windows, you will need RTools 4.2 as well. Your previous RTools installation will not be sufficient.
--->


### Download materials for this workshop.

There are two options here: 

a. You can download the repository from https://github.com/jfieberg/Madeira-Workshop. To do so, go to the URL, click on the green `< > Code` button, and select "Download ZIP". Then, unzip the files to a known location on your computer.

b. If you are familiar with working with git and github, you can clone the repository instead.
    

### Install packages and run test script

The script `PackagesInstall.R` in the `Organizing-Intro-Documents` subfolder can help you install and/or update the packages required for this workshop.  Please make sure to run it.

Once you have downloaded the course materials, open Rstudio and create a project associated with the directory holding all the workshop files. `File -> New Project`, then select `Existing directory`, find the `Madeira-Workshop` directory` on your computer and click OK.

Then, in Rstudio, open the `Test.qmd` file contained in the `Organizing-Intro-Documents` subfolder. Click `Render` (near the top of the screen in Rstudio).  This will:

- Run the R code in the `Test.qmd` file
- Ensure that your computer is ready to go for Monday when the workshop begins.

If you are successful, this will create a `Test.hmtl` file in the same subdirectory (Organizing-Intro-Documents). Upload this file to [this google drive](https://drive.google.com/drive/folders/19PJ11vSP2bQ8ObEj12fP0y_CcWcI139I) to demonstrate you are ready to go at the start of the workshop!

<!---

## Other Resources



 
**ESA Ecological Forecasting Initiative**: webinar on iSSA by Tal Avgar and Brian Smith. You can find a [recording of the webinar on YouTube](https://youtu.be/jiY9N-TNRjs). You can find the [lecture slides, R code, and Q&A on GitHub](https://github.com/eco4cast/Statistical-Methods-Seminar-Series/tree/main/avgar-smith_issa). You can find the Q&A markdown in the GitHub repo, or [just follow this link](https://github.com/eco4cast/Statistical-Methods-Seminar-Series/blob/main/avgar-smith_issa/Q_and_A.md).

The webinar includes some coded examples of iSSFs that include interactions with the movement parameters that we did not demonstrate in this workshop.
--->

