# Microbiome

## Description of Files (9/19)

* **XiaChenFungLiModel.R:** Complete implementaion of the MC-EM model from Xia et al. (2013)
* **XiaTestRun.R:** Imports other scripts and demonstrates how to run the MC-EM function. Also includes code for generating the fitting graphs, as in Figure 4 Xia et al.
* **getData.R:** Loads in raw data files provided by Amy. Includes functions to purturb, take logratio, and select most prevelent OTUs.
* **run0A.R:** Shows how to select a subset of the data based on experimental conditions. (Not currently used.)
* **HQQModel.R:** Implementation of a simple posterior predictive Dirichlet-Multinomial distribution.
* **fittingPlots.R:** Compares DayAmdmt 11 and 21, contains code for various plots.

## Notes

* Currently implemented to work specifically on my desktop with filepaths.
* Includes parallelization, make sure to check if run on different machines.
