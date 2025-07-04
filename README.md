# GENERAL INFORMATION
This README.txt file was updated on 21 May 2025

## A. Paper associated with this archive
Citation: Omitted for double-blind peer review

Brief abstract: 

Coevolution requires reciprocal genotype-by-genotype (GXG) interactions for fitness, which occur when the fitness of a genotype in one species depends on the genotype it interacts with in another species, and vice versa. However, in mutualisms, when GXG interactions are mutually beneficial, simple models predict that the resulting positive feedbacks will erode genetic variation. Here, we explore how genotype-by-environment (GXE) interactions, which occur when the fitnesses of different genotypes respond differently to different environments, and GXGXE interactions, which occur when the environment changes the outcome of GxG interactions, maintain variation in mutualisms. We build a spatial population genetic model in which the fitnesses of two partners depend on mutually-beneficial GXG, GXE, and GXGXE interactions. Our analysis shows that variation will always be maintained via migration-selection balance with stronger GXE than GXG interactions. However, strong GXG interactions can erode variation by allowing genotypically-matched partners to fix and, more surprisingly, weak GXG interactions can erode variation by allowing genotypically-mismatched partners to fix at high dispersal rates, leading to apparent maladaptation between partners. We parameterize our model using data from three published reciprocal transplant experiments, infer the relative strengths of GXE and GXG, and discuss the implications for the maintenance of genetic variation.

This repository contains the Mathematica notebooks for full analysis of the model and numerical simulations, as well as the original data from a published reciprocal transplant experiment, and R and Python scripts for data analysis and visualization.

## B. Originators
Christopher I. Carlson, Matthew M. Osmond, Megan E. Frederickson

"Johnson_data.csv" was provided via email from Nancy Collins Johnson (Nancy.Johnson@nau.edu) and Gail Wilson (gail.wilson@okstate.edu) as original data from the paper listed under REFERENCES, below. 

## C. Contact information
Christopher I. Carlson
25 Willcocks St, Toronto, ON M5S 3B2, Canada
christopheriancarlson@gmail.com

## D. Dates of data collection
See original publications, cited in REFERENCES, below.

## E. Geographic Location(s) of data collection
See original publications, cited in REFERENCES, below.

## F. Funding Sources
This research was supported by Natural Sciences and Engineering Research Council of Canada (NSERC) Discovery Grants to M.E.F. and M.O., a NSERC Discovery Accelerator Supplement to M.E.F., and a University of Toronto Connaught International Scholarship to C.C. 

## ACCESS INFORMATION
1. Licenses/restrictions placed on the data or code
CC BY-SA 4.0 Creative Commons Attribution-ShareAlike 4.0 International

2. Data derived from other sources
Please contact Nancy Collins Johnson (Nancy.Johnson@nau.edu) and Gail Wilson (gail.wilson@okstate.edu) to reuse "Johnson_data.csv"

3. Recommended citation for this data/code archive
n/a 

## DATA & CODE FILE OVERVIEW
This data repository consists of 1 data file, 4 code scripts, and this README document, with the following data and code filenames and variables:

Data files and variables
1. Johnson_data.csv:
SOIL = soil source site, A = Fermi labs, K = Konza, and C = Cedar Creek. 
PLANT = plant source site, A = Fermi labs, K = Konza, and C = Cedar Creek. 
INOC = mycorrhizal inoculum source cite, A = Fermi labs, K = Konza, and C = Cedar Creek. Site A, K, or C preceded by an S indicate a sterile application of inocula. 
Hyphae =  hyphal density, in m g soil-1. 
Plant_Biomass = plant aboveground biomass, in g. 

Code scripts and workflow
1. GXE Analysis.nb: Comprehensive Mathematica notebook containing full analysis of our model along with numerical simulations. 

2. nuismer.nb: Mathematica notebook containing a brief supplemental analysis containing the same assumptions as the model from Chapter 8 of Nuismer (2017), referenced in REFERENCES, below. 

3. Henry_oviposition_analysis.Rmd: File containing analysis of oviposition data from Henry et al. (2025), which measured the number of eggs laid per fly for flies evolved in control or high sugar media with and without microbes. Their original data is available ("2536_aceto_ovidata.csv") at the Henry et al. (2025) Dryad repository cited in REFERENCES, below. For details on exact statistical analysis see Henry et al. (2025). We used the statistical analysis from this notebook when stating whether differences between fitness proxies was significant. We used means and standard errors from the "m1.emmB" object to generate reaction norms, using our python code "GXE_delta_norms.py" to generate plots. 

4. GXE_delta_norms.py: Python file containing code used to generate reaction norm plots using data from Johnson et al. 2010, Ricks et al. 2023, and Henry et al. 2025 (see REFERENCES, below). 

For the Johnson et al. (2010) study: We imported data directly from "Johnson_data.csv" to generate reaction norm plots for that study. 

For the Ricks et al. (2023) study: Ricks et al. (2023) built aster models to combine multiple fitness components to estimate  fitness. We used means and standard errors (squared because they had been square-root transformed) from their analysis (from the "aster_sum" object generated by their models.R script and "greenhouse_results.csv" from the Ricks et al. (2023) data repository referenced below) corresponding to our model's fitness advantages into our python code "GXE_delta_norms.py" to generate reaction norm plots for this interaction. We used statistical analysis of the aster models in Ricks et al. (2023) to denote significant fitness advantages. 

For the Henry et al. (2025) study: See above description of Henry_oviposition_analysis.Rmd.

## SOFTWARE VERSIONS
R version 4.3
Python 3.12
Mathematica versions 9.0 and 14.1

R packages: tidyverse 2.0.0, car 3.1-2, DHARMa 0.4.6, glmmTMB 1.1.9, lme4 1.1-35.3, emmeans 1.10.1, performance 0.12.4, glmm 1.4.5
Python libraries: numpy 1.24.2, scipy 1.10.1, csv 1.0, os 4.3, matplotlib 3.7.1

## REFERENCES
L. Henry, M. Fernandez, A. Webb, & J. Ayroles. Data from: Microbial solutions to dietary stress: Experimental evolution reveals host-microbiome interplay in Drosophila melanogaster [dataset]. Dryad (2025).

L. Henry, M. Fernandez, A. Webb, & J. Ayroles.  Microbial solutions to dietary stress: experimental evolution reveals host–microbiome interplay in Drosophila melanogaster. Proceedings of the Royal Society B: Biological Sciences, 292(2043) 20242558 (2025).

N.C. Johnson, G.W.T. Wilson, M.A. Bowker, J.A. Wilson, & R.M. Miller, Resource limitation is a driver of local adaptation in mycorrhizal symbioses, Proc. Natl. Acad. Sci. U.S.A. 107 (5) 2093-2098, https://doi.org/10.1073/pnas.0906710107 (2010).

S.L. Nuismer, Introduction to Ceovolutionary Theory. W.H. Freeman and Company, New York (2017). 

K. Ricks, Contribution of microbial interactions to the local adaptation of Bromus tectorum to variable salinity conditions [dataset]. Dryad, https://doi.org/10.5061/dryad.4f4qrfjdq (2023).

K.D. Ricks, N.J. Ricks, & A. Yannarell, Patterns of plant salinity adaptation depend on interactions with soil microbes. The American Naturalist, 202(3) 276–287 (2023). 




