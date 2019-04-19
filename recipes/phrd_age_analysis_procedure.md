# `phrd_age_analysis_procedure`

**Authors: M. S. Povich (@mattpovich), J. T. Maldonado, & E. H. Nuñez**

Step-by-step guide to producing pHRDs and mass-age histograms
This recipe must be run AFTER the SED fitting procedure, see
sedfitting_procedure_xpms.txt 

Current version:
v1.0 2019 April 12 - ORIGINAL public version

Custom IDL library codes called (make sure these are on your IDL
!PATH):
PDISTHRD, PROBDIST2D, OPLOT_SIESS, MAD_HIST2, FIND_MPCM, FIND_CLUSTSRC, SIESS_ZAMS.TXT

------------------------------------------------------------------------------------------------------------

Making pHRDs, Mass and Age Histograms

   Launch IDL from your $TARGET_DIR/sedfitter directory. The lines below can be copied and pasted line-by-line or all at once into an IDL terminal.
	
;****************
;pHRDs. These will by default be placed within the target_pms directory
;with sourcelist_pms as the filename prefix

  target_pms = 'results_xpms'  ;Naming convention, these are X-ray selected PMS stars
  sourcelist = 'sourcelist.txt' ;Default name, various source subsets can be created and saved with alternate names in the target_pms directory.
  region_name = '<enter short target name>' ; e.g. 'M17'

;Plot the composite pHRD for ALL stars in sourcelist.txt. If pHRDs for individual stars are desired, set psplots=3 instead.

  pdisthrd_pms,target_pms,psplots=1,restrict=3,sourcelist=sourcelist,region_name=region_name

;Plot mass and age histograms, find best isochrones. This may require interactive selection of the best bin to characterize the age distribution. This is accomplished by adjusting the MASS_OFFSET[12] and AGE_OFFSET values in the command below:

  massage_plots_pms,target_pms,sourcelist=sourcelist,region_name=region_name,mass_offset1=0,mass_offset2=0,age_offset=0
;****************


END OF RECIPE

