Estimates of forest stand age distributions through time based upon land-use change and disturbance datasets

## Author
Thomas A. M. Pugh, University of Birmingham, UK, t.a.m.pugh@bham.ac.uk


## Description
The scripts herein underlie the forest stand age distribution calculations in Figures 2 and 3 of the following paper:
McDowell, N. et al. (2020) Pervasive shifts in forest dynamics in a changing world. 


## Methods

1) Estimates of the influence of human activities only on forest stand age distribution in all forests.

Calculations were based upon a historical dataset of land-use change (LUHv2). Stand age was tracked from 1 to 140 years, with stands older than 140 years assigned to an old-growth category (Pugh et al., 2019a). Forest cover fractions were initialised in 1750 according to the “primary” and “secondary” categories in LUHv2. All forest cover was initialised as old-growth. Each year, all forest stand ages between 1 and 149 years were incremented by one year. Then annual forest cover losses, as defined in LUHv2h for both land-use change and wood harvest, were subtracted from the existing forest fractions. In the baseline estimate, an equal probability of loss across all stand age classes was assumed. Sensitivity simulations were also made in which the probability of disturbance changed linearly with age, such that,
(a) the probability of an old-growth (nominally 140 y old) stand being transitioned was five times that of a 1-year old stand,
(b) the probability of a 1-year old stand being transitioned was five times that of an old-growth stand.
Forest cover gain given by LUHv2 was assigned to the youngest “secondary” age class. Total forest area was based on the “forested” fractions given in the LUHv2 dataset. LUH2 “primary” and “secondary” classes were maintained throughout the calculations to avoid losing transition information, but grouped together in the outputs. These calculations were run annually between 1750 and 2015. Given that old-growth is defined as greater than 140 years, at least 140 years of simulation are required to generate a fully-developed age distribution. Therefore, output for the period 1750-1900 was discarded as spin-up. The final age distributions are grouped to decadal age classes for display purposes. All calculations were made at 1° x 1° spatial resolution.
To make these calculations set, inc_luh=true, inc_luh2futscen=false, inc_woodharv=true, inc_dist=false, inc_agesens=true, inc_distsens=false, ccanopymask=false, in age_class_reconstruction.m.

2) Estimates of forest stand age distribution from human activities and natural disturbances in closed-canopy forests.

The same basic method as in (1) was employed, with the following differences:
- At each timestep, after the calculation of the effect of land-use change on forest age based on LUHv2, an additional perturbation representing natural and wood harvest disturbances was applied. These used a generic gridded disturbance rate (Pugh et al., 2019b). This disturbance rate is representative of the period 2000-2014 and was applied throughout in the absence of information on how disturbance rates other than due to LUC have changed over time. The same assumptions and sensitivity tests for how disturbance likelihood scales with age were used as described for the LUHv2 transitions.
- Forest transitions in the LUHv2 dataset due to wood harvest were not included to avoid double-counting with the disturbance dataset. 
- All results are presented for closed-canopy forest area only, defined as at least 50% canopy cover at 0.01° x 0.01° resolution (Pugh et al., 2019b). Land-use changes based on LUHv2 were assumed to be proportionally distributed across open- and closed-canopy forest within a grid cell.
- Simulations were extended to 2100 using the GCAM RCP 3.4 scenario from LUHv2f.
- In a sensitivity simulation, non-land-use-change disturbance rates were incremented linearly to 200% of the 2001-2014 values over the period 2015 to 2050 and held constant at that level thereafter.
To make these calculations set, inc_luh=true, inc_luh2futscen=true, inc_woodharv=false, inc_dist=true, inc_agesens=true, inc_distsens=true, ccanopymask=true, in age_class_reconstruction.m.


## Files

age_class_reconstruction.m - master script for the reconstruction. Options allow simulations in the above paper to be recreated.
age_class_reconstruction_func.m - function making the reconstruction calculation
luh2_forstates_read.m - function reading forest state from LUHv2
luh2_forloss_read.m - function reading forest loss transitions from LUHv2
luh2_forgain_read.m - function reading forest gain transitions from LUHv2
make_figure_2.m - script to make Figure 2 in the above paper
make_figure_3.m - script to make Figure 3 in the above paper
olson_biom_to_raster.m - script to aggregate the Olson biomes to the rasters used in age_class_reconstruction.m
global_grid_area_1deg.m - function to calculation area of grid cells.

Locations of data files required for the calculations are given in the individual scripts.


## References

Pugh, T.A.M., Lindeskog, M., Smith, B., Pouleter, B., Arneth, A., Haverd, V., Calle, L. (2019a) Role of forest regrowth in global carbon sink dynamics. Proc. Nat. Acad. Sci. USA 116(10), 4382-4387.

Pugh, T.A.M., Arneth, A., Kautz, M., Poulter, B., Smith, B. (2019b) Important role of forest disturbances in the global biomass turnover and carbon sinks. Nature Geoscience, doi: 10.1038/s41561-019-0427-2


05.01.20
