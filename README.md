# KokDustDRE2017
This repository holds the code used to generate the results in Kok et al., Nature Geoscience, 2017.

The main functions in this repository of MATLAB code are main_bootstrap_global_dust_cycle.m and main_calc_DRE.m. The first function obtains constraints on the dust extinction efficiency (Fig. 1b), the emitted dust size distribution (Fig. 1c), and the atmospheric dust lifetime (Fig. 1d). It combines those results with constraints on the global dust aerosol depth from Ridley et al. (ACP, 2016; see Fig. 1a) to constrain the size-resolved global dust loading (Figs. 2 and 3), as well as the size distribution of the global dust optical depth (Fig. 2c). The second function, main_calc_DRE.m, then combines those constraints with global model simulations of the dust radiative effect per unit dust optical depth to obtain constraints on the global direct radiative effect of dust aerosols. All results reported in Kok et al. (Nature Geoscience, 2017) can be reproduced by calling main_bootstrap_global_dust_cycle.m and main_calc_DRE.m. 

Below follows a brief description of every script, function, and data set included in this code.

Scripts and functions:
•	Bin.m: function that sorts the data in var into bins of size bin_size, and writes out the result to filetext
•	Calc_AeroCom_to_DRE_model_AOD.m: script that calculates the AOD of each AeroCom model bin when mapped onto the bins of the four global models that calculated the radiative effect efficiency. This script also calculates the DRE from the size-resolved dust loading of each AeroCom model.
•	Calc_PSD_corr_fact.m: this function calculates the factor by which each data sets needs to be corrected to maximize agreement against the fit, which is necessary because the non-continuity of dust PSD measurements in D make it such that there is no single normalization factor
•	dustPSD_power_law_fit.m: function that fits the dust PSD data to a power law following Kok (PNAS, 2011)
•	Fit_PSD_and_lifetime_functions.m: Function that fits the lifetime and emitted dust PSD data sets using the statistical model described in the supplement to Kok et al. (2017), obtaining the MLE
•	Gaussian.m: draws a realization from the normal distribution with specified mean and SD
•	linearfit.m: this function performs a linear fit to the x, y, and y_err data
•	load_lifetime_data.m: function that loads the model data on the dust lifetime
•	load_PSD_data.m: function that loads the experimental data on the emitted dust PSD
•	main_bootstrap_global_dust_cycle_size.m: Function that uses a bootstrap procedure to constrain the size-resolved global dust loading, using Eqs. (1) – (4) in Kok et al. (2017).
•	Main_calc_DRE.m: function that obtains constraints on the global dust DRE, using the size-resolved dust loading obtained through main_bootstrap_global_dust_cycle_size.m.
•	Mie.m, Mie_ab.m, Mie_pt.m, Mie_S12.m, Miecoated.m, Miecoated_ab1.m, Miecoated_S12.m: Mie code from http://www.hiwater.org/Mie_calcs_files/Miecalc-Maetzler-Subset.zip (last accessed on 1/27/2017), which is based on code in Bohren and Huffman (1983). Citation for reuse: Mätzler, C., MATLAB Functions for Mie Scattering and Absorption, Institut für Angewandte Physik, Research Report No. 2002-08, Bern, Switzerland, 2002.
•	Write_data.m: Script that writes the important results from main_bootstrap_global_dust_cycle_size.m to files

Data sets:
•	Aerocom_data.mat: data file containing the size-resolved loading and DAOD of AeroCom models
•	DAOD.mat: data file containing the pdf of the DAOD, after Ridley et al. (2016)
•	Lifetime_and_PSD_fits.mat: data file containing the fits to the lifetime and PSD data
•	REE.mat: data file containing the results of the radiative effect efficiency from the four global climate models
•	TAMU_data.mat: data file with extinction coefficient info from the TAMU single-scattering database (Meng et al., 2010)
•	TAMU_database_parameters.mat: data file containing the parameters for which the TAMU database was generated

