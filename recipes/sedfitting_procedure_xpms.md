sedfitting_procedure_xpms.txt
Original 2011-11-22 by M. S. Povich
UPDATED to use PYTHON sedfitter software 2018-08-28 by M. S. Povich

PRODUCTION (V1.0)  2019-04-17 M. S. Povich

This recipe guides you through the steps of performing a basic SED fitting classification for the combined Spitzer, 2MASS, and UKIDSS COUNTERPART SEDs to X-ray sources in a specified target field.

The final science product is a directory called results_xpms
containing IDL save files of the SED model parameters for each source
in the population.

=============================================================================
INITIAL SETUP 

Get Tom Robitaille's sedfitter software and install it following the instructions
at http://sedfitter.readthedocs.io/en/stable/installation.html

Download the custom pms SED model set (constructed using Siess et al. 2000 and Bernasconi & Maeder 1996 evolutionary tracks) from <ZENODO URL>.

#I am using this version of python, which I need to
 run in a bash shell:

 Python 3.6.5 :: Anaconda, Inc.
 
 IDL commands will need to be executed in tcsh,
 so we require two different terminal windows open in different shells.
 
  ### Python commands/blocks are preceded by ">>>"
  ### tcsh commands/blocks are preceded by "%"

>>>
TARGET='m17'  #python EXAMPLE. Choose your own target name.

%
setenv TARGET m17  # tcsh EXAMPLE. Choose your own target name (must match above)

#It is a good idea to name your working directory after the target and then cd into it in both shells! E.g. ~/sedfitting/m17

%
mkdir sedfitter
cd sedfitter 


# PREPARE SOURCE PHOTOMETRY FOR SED FITTING

You are responsible for assembling a sample of young stellar sources with a distance and extinction distribution that can be reasonably well constrained. My own preferred method is X-ray selection plus parallax-based cleaning of remaining field-star contaminants using Gaia (see Povich et al. 2019 for details).

Prepare the input fitter data file following the guidelines at https://sedfitter.readthedocs.io/en/stable/data.html. This recipe assumes you have a total of n=10 wavelengths/filters, in this order:

UKIDSS/VVV: JHK (filters 1-3)
2MASS:  JHKs (filters 4-6)
Spitzer/IRAC: [3.6], [4.5], [5.8], [8.0] (filters 7-10)

Other surveys could be substituted (e.g. WISE for Spitzer/IRAC), but the ORDER that photometry points appear in the data file should not be changed.

Data file name: data_xir

Estimate the maximum reddening to stars in the sample using the J-H vs. H-K color-color diagram.

%
IDL>
   	data = 'data_xir'
	twomass = 0   ; Set to 1 if your JHKs colors are on the 2MASS system
	if not twomass then mk = 1
	magfromdata,data,0,j,nwav=10,mk=mk
	magfromdata,data,1,h,nwav=10,mk=mk
	magfromdata,data,2,k,nwav=10,mk=mk
	plot_nircc_rv,j,h,k,twomass=twomass
	
;Estimate the maximum reddening in Av by comparing the locus
  of stars to the reddening vectors (marked at intervals of Av=5
  mag). Note the value that you used.
  
=============================================================================
FIT WITH DISKLESS PMS MODELS

  >>>
  python
  from astropy import units as u
  from sedfitter import fit
  from sedfitter.extinction import Extinction

	# Define path to models
  model_dir_pms = <your path here> + '/models_pms'

	# Read in extinction law. The two extinction laws that I used are packaged in the ex_law subdirectory of this repository. COPY ex_law_gl.par into your sedfitter directory!

  extinction = Extinction.from_file('ex_law_gl.par', columns=[0, 3], wav_unit=u.micron, chi_unit=u.cm**2 / u.g)

  # Define filters and apertures. Apertures are tricky but the PMS models are aperture-independent, so it doesn't really matter!
  # Note these filter names are SPECIFIC to the pre-convolved model SEDs. These names must match the filenames in the models_pms/convolved directory.
  # The example below is for UKIDSS, 2MASS and IRAC. VISTA[ZYJHK] are also available.

  filters = ['UDSSJ', 'UDSSH', 'UDSSK', '2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4']
  apertures = [1.5, 1.5, 1.5, 3., 3., 3., 3., 3., 3., 3.] * u.arcsec

# Fit!
# Make sure the distance_range parameter is set to the minimum and maximum distance, in kpc, to your target source population, and the av_range reflects the maximum extinction estimated from the JHK color-color diagram above.

  fit('data_xir', filters, apertures, model_dir_pms, 'xpms.fitinfo', n_data_min=4, extinction_law=extinction, distance_range=[1.5, 2.] * u.kpc, av_range=[0., 25.], output_format=('F',1))

###Split the output into well-fit versus poorly-fit. I strongly recommend using the chi^2/Ndata = 1 cutoff for well-fit models, see Povich et al. (2019).

  >>>
  from sedfitter import filter_output
  filter_output('xpms.fitinfo', cpd=1.) 

  ### OPTIONAL Plot up some SEDs for sanity checks.
  ##CANNOT DO THIS AT PRESENT, SED .fits files for PMS models are in the wrong format for the python fitter!##
  >>>
  from sedfitter import plot
  plot('xpms.fitinfo_good', 'plots_xpms_good') 
  plot('xpms.fitinfo_bad', 'plots_xpms_bad')

  ######

  from sedfitter import write_parameters
  write_parameters('xpms.fitinfo_good','pars_xpms_good.txt',select_format=('A',-1))

  # Make fitter data file and SAOImage ds9 regionfile of the well-fit sources, and create results_xpms directory contining model fitting results as IDL save files.
  %
  IDL>
  parfile_good = 'pars_xpms_good.txt' 
  data_parent='data_xir'
  fitinfo2data,parfile_good,data_parent,'data_xpms_good'
  ds9regfromdata,'data_xpms_good','data_xpms_good.reg',color='dodgerblue'

  target_pms ='results_xpms'
  pms_pars2idl,parfile_good,target_pms,data_parent=data_parent,filter='F1'	
  
### NEXT STEPS:
  (I) results_pms is starting point for advanced analysis of the PMS model results. See phrd_age_analysis_procedure.txt.
  (II) Sources that were poorly fit by "naked" PMS models may be YSOs. Follow our recipe sedfitting_procedure_xyso.txt to fit YSO models from Robitaille (2017) to the sources in xpms.fitinfo_bad. 


END OF RECIPE