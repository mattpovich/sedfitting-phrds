# `phrd_age_analysis_procedure`

**Authors: M. S. Povich, J. T. Maldonado, & E. H. Nuñez**

## Description

Step-by-step guide to producing "probabilistic" H-R diagrams (pHRDs) and associated mass-age histograms. Follow this recipe only after completing the SED fitting recipe described in `sedfitting_procedure_xpms`.

## Version History
* CURRENT (v1.0) – 2019-04-12; original public version


## Making pHRDs, Mass and Age Histograms

*Run IDL with $TARGET/sedfitter as your working directory.*
	
The SED fitting recipe `sedfitting_procedure_xpms` created a subdirectory `results_xpms` containing one IDL save file for each source well-fit by the `models_pms` SED model set. The file `results_xpms/sourcelist.txt` is a simple list of all the IDL save file names. You may create new `sourcelist_new.txt` files containing various subsets of interest by directly editing `sourcelist.txt` (but be warned that lists containing <100 sources may not produce reliable results!). All plots created using a given `sourcelist*.txt` file will be given
`sourcelist*` as their filename prefixes.

**IDL>**

  	target_pms = 'results_xpms'  ;Naming convention; these are X-ray selected PMS stars
  	sourcelist = 'sourcelist.txt' ;Default name, can be replaced by a subset of the original list.
  	region_name = '<short target name>' ; e.g. 'M17'

Plot the composite pHRD for ALL stars in sourcelist.txt. 
If pHRDs for individual stars are also desired, set `psplots=3` instead.
If you are NOT using an X-ray selected source sample, set `restrict=2` instead.
	
  	pdisthrd_pms,target_pms,psplots=1,restrict=3,sourcelist=sourcelist,region_name=region_name

The above command produces a plot of the pHRD in `sourcelist*_composite_HRD_agecons.eps`.

Plot mass and age histograms, finding the best isochrone for the population in the process. This may require interactive selection of the optimal histogram bin to characterize the age distribution. Experiment with adjusting the `mass_offset[12]=` and `age_offset=` keyword values and reviewing the effect on the printed and plotted output of the command below:

  	massage_plots_pms,target_pms,sourcelist=sourcelist,region_name=region_name,mass_offset1=0,mass_offset2=0,age_offset=0

*Explanatory notes:*
* The above command produces two plots: 
	1. `sourcelist*_MassAgeHist_m.m-o_A.A_ISO.eps` — `m.m` encodes the value of *MC1*, the minimum mass (Msun) of stars used for the age determination; `o` the value of `mass_offset1`; and `A.A` the characteristic age in Myr (Povich et al. 2019).
	1. `sourcelist*_MassAgeHist_M.M-O_A.A_ISO.eps` — `M.M` encodes the value of *MCS*, the cutoff mass (Msun) below which the final model mass distribution turns over from a Salpeter power law; `O` the value of `mass_offset2`; and `A.A` the adopted isochronal age in Myr (Povich et al. 2019).
* `mass_offset1` specifies the location of the histogram bin *MC1* (Povich et al. 2019), counting from the peak of the distribution to the right (increasing mass) in plot (1).
* `mass_offset2` specifies the location of the histogram bin *MC2* (Povich et al. 2019), counting from the peak of the distribution to the right (increasing mass) in plot (2).
* `age_offset` specifies the location of the histogram bin giving the optimal isochronal age of the population (tau_SF). This should *only* be set differently from 0 if the model age distribution exhibits an exceptionally broad peak.


