# MultiCity_Squirrels

A repository that contains the data and code for:

Larson, R.N., H. A. Sander, M. Fidino, J. L. Angstmann, L. Barczak, A. Davidge, D. Drake, L. Hartley, K. Moore, P. R. Sanchez, A. Robey, C. Salsbury, E. Sadkin, T. Snyder, T. Stankowich, K. Tombs, D. Will, J. Williamson, A. J. Zellmer, and S. B. Magle. City-specific responses to local land cover alter spatial overlap among tree squirrels in the United States. *Global Change Biology*. (In prep)


This `README` file includes information on the various scripts and datasets used for this analysis. Not every data source is saved in this repository (e.g., GIS data). The manuscript contains citations for where to access the geospatial data.

---

<div align="center"> <h3>data</h3> </div>

---

There are 5 files in here. These files generally use a four-character code for each city, which are:

| City                                | Code    |
| ----------------------------------- | ------- |
| Austin, Texas                       | `autx`  |
| Chicago, Illinois                   | `chil`  |
| Denver, Colorado                    | `deco`  |
| Indianapolis, Indiana               | `inin`  |
| Iowa City, Iowa                     | `ioio`  |
| Madison, Wisconsin                  | `mawi`  |
| Los Angeles Metro Area, California  | `mela`  |
| Rochester, New York                 | `rony`  |
| Wilmington, Delaware                | `wide`  |

**./data/FS.csv:** The detection/non-detection data for fox squirrels
**./data/GS.csv:** The detection/non-detection data for eastern gray squirrels
**./data/RS.csv:** The detection/non-detection data for red squirrels

| Column | Type      | Description |
| ------ | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Site   | character | The code for the site name                                                                                                                                    |
| City   | character | The city code for a given city                                                                                                                                |
| ju18_1 | numeric   | The detection/non-detection (i.e., 1/0/NA) data for that species for that site for the first week of sampling in July 2018                                    |
| ju18_2 | numeric   | The detection/non-detection (i.e., 1/0/NA) data for that species for that site for the second week of sampling in July 2018                                   |
| ju18_3 | numeric   | The detection/non-detection (i.e., 1/0/NA) data for that species for that site for the third week of sampling in July 2018                                    |
| ju18_4 | numeric   | The detection/non-detection (i.e., 1/0/NA) data for that species for that site for the fourth week of sampling in July 2018                                   |
| oc18_1 | numeric   | The detection/non-detection (i.e., 1/0/NA) data for that species for that site for the first week of sampling in October 2018                                 |
| ...    | ...       | Continuation of the detection/non-detection data. Seasons are abbreviated JA = January (winter), AP = April (spring), JU = July (summer), OC = October (fall). The number following a season abbreviation is the year of sampling (e.g., AP19 = April 2019) |

**./data/cooc_covs.csv:** This contains site-level covariates. 

| Column         | Type      | Description |
| ------         | --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Site           | character | The code for the site name                                                                                                                               |
| City           | character | The city code for a given city                                                                                                                           |
| meanCan        | numeric   | The mean percentage of tree canopy cover for the metropolitan region in which a site is located                                 |
| Canopy         | numeric   | The percent tree canopy cover for each site                                    |
| gmcCanopy      | numeric   | The mean-centered tree canopy cover for each site (i.e., Canopy - meanCan)                                    |
| meanImpervious | numeric   | The mean percent impervious cover for the metropolitan region in which a site is located                                |
| Impervious     | numeric   | The percent impervious canopy cover for each site                                  |
| gmcImpervious  | numeric   | The mean-centered impervious cover for each site (i.e., Impervoius - meanImpervious |
| meanGrass      | numeric   | The mean percent grass cover for the metropolitan region in which a site is located |
| Grass          | numeric   | The percent grass cover for each site |
| gmcGrass       | numeric   | The mean-centered grass cover for each site (i.e., Grass - meanGrass) |

**./data/cooc_obsCovs.csv:** This contains observation-level covariates.

| Column | Type      | Description |
| ------ | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Site   | character | The code for the site name                                                                                                                                    |
| City   | character | The city code for a given city                                                                                                                                |
| ju18_1 | numeric   | The number of days the camera was active at this site during first week of sampling in July 2018                                    |
| ju18_2 | numeric   | The number of days the camera was active at this site during second week of sampling in July 2018                                   |
| ju18_3 | numeric   | The number of days the camera was active at this site during third week of sampling in July 2018                                    |
| ju18_4 | numeric   | The number of days the camera was active at this site during fourth week of sampling in July 2018                              |
| oc18_1 | numeric   | The number of days the camera was active at this site during first week of sampling in October 2018                                |
| ...    | ...       | Continuation of the camera data. Seasons are abbreviated JA = January (winter), AP = April (spring), JU = July (summer), OC = October (fall). The number following a season abbreviation is the year of sampling (e.g., AP19 = April 2019) |

---

<div align="center"><h3>JAGS</h3></div>

---

This contains 5 files. These files are code for each Bayesian multi-city multi-species occupancy model we fit to the data.

**./JAGS/multiSpec_PCA.R:** Model using PCA-derived variables as occupancy covariates

**./JAGS/multiSpec_Can.R:** Model using tree canopy cover as an occupancy covariate

**./JAGS/multiSpec_Grass.R:** Model using grass cover as an occupancy covariate

**./JAGS/multiSpec_Imperv.R:** Model using impervious cover as an occupancy covariate

---

<div align="center"><h3>R</h3></div>

---

This contains 1 subfolder, `functions`, and two files.

**./R/CoOccurrence_Model.R:** R script for running the models

**./R/graphing.R:** R script for producing Figures 2 - 4

<h2>functions</h2>

**./R/functions/wide_to_stacked.R:** Script for taking wide-format detection histories and stacking them.

**./R/functions/split_mcmc.R:** Script for splitting the output of mcmc chains into a matrix. Written by @MFidino

**./R/functions/my_inits.R:** Script that stores initial values for JAGS model implementation

**./R/functions/calc_cpo.R:** Script for calculating the conditional predictive coordinate (CPO) of each model. Written by @MFidino
