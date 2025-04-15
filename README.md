# Madeira-Workshop

Presentations, code, and data associated with the workshop titled, Analysing animal movement in the marine environment, held in Madeira, Portugal from May 19-23, 2025.

## Preliminary Schedule

General format:

- We will cover a different topic on each day, starting at 9:00am and ending at 4:30.  
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
- Limitations of discrete-time approach

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


### R Packages

<!--- This part still needs some work --->

The script `packages.R` in the `Organizing-Intro-Documents` subfolder can help you install and/or update the packages required for this workshop.  Please make sure to run it.

Once you install the packages needed for the workshop, please either: 

- Open `Test.R` (contained in the `Organizing-Intro-Documents` subfolder). Then, in Rstudio, go to the File menu and select File -> Compile Report -> html
- Or, knit/compile `Test.rmd` or `Test.qmd` to create the same html file

Then, upload the resulting Test.hmtl file to [this google drive - still need to provide link](provide link here!)

 
## Other Resources


<!--- 
**ESA Ecological Forecasting Initiative**: webinar on iSSA by Tal Avgar and Brian Smith. You can find a [recording of the webinar on YouTube](https://youtu.be/jiY9N-TNRjs). You can find the [lecture slides, R code, and Q&A on GitHub](https://github.com/eco4cast/Statistical-Methods-Seminar-Series/tree/main/avgar-smith_issa). You can find the Q&A markdown in the GitHub repo, or [just follow this link](https://github.com/eco4cast/Statistical-Methods-Seminar-Series/blob/main/avgar-smith_issa/Q_and_A.md).

The webinar includes some coded examples of iSSFs that include interactions with the movement parameters that we did not demonstrate in this workshop.
--->

