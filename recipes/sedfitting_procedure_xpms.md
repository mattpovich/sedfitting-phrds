# `sedfitting_procedure_xpms`

## Description
This recipe guides you through the steps of performing a basic SED fitting classification for the combined Spitzer, 2MASS, and UKIDSS counterpart SEDs to X-ray sources in a specified target field. The final science product is a directory called `results_xpms`
containing IDL save files of the SED model parameters for each source in the target field.

## Version History

* Original `sedfitting_procedure_xpms.txt` – 2011-11-22 by M. S. Povich
* Updated to use  [`python sedfitter`](https://github.com/astrofrog/sedfitter)  – 2018-08-28 by M. S. Povich
* CURRENT (v1.0)  2019-04-17 M. S. Povich

## Important Notes 
*This recipe requires both `python` and `IDL`. I have verified the software using `Python 3.6.5` run in the `bash` shell and `IDL 8.7.1` in the `tcsh` shell via the MacOS X 10.13.6 Terminal app. I have NOT yet tested it on other platforms or other software versions, but I suspect `python` may be pickier about version compatibility than `IDL` or MacOS X.*

 I recommend using two different terminal windows simultaneously, one in `bash` and the other in `tcsh`. The correct shell in which to execute each block of command is indicated as follows:
 
  * `python` command blocks are preceded by **>>>**
  * `tcsh/IDL` command blocks are preceded by **%/IDL>**


## INITIAL SETUP 

Get Tom Robitaille's sedfitter software and install it following the instructions
at http://sedfitter.readthedocs.io/en/stable/installation.html

Download the custom `models_pms' SED model set from **GET ZENODO URL**

**>>>**

    TARGET='m17'  # python EXAMPLE. Choose your own target name.


**%**

    setenv TARGET m17  # tcsh EXAMPLE. Choose your own target name (must match above)


It is a good idea to name your working directory `$TARGET/sedfitter`, e.g. (in `tcsh`):

**%**

    mkdir -p $TARGET/sedfitter
    cd $TARGET/sedfitter


## PREPARE SOURCE PHOTOMETRY FOR SED FITTING

You are responsible for assembling a sample of young stellar sources with distributions of distance and extinction that can be reasonably well constrained by independent metrics. My own preferred method is X-ray selection plus parallax-based cleaning of remaining field-star contaminants using Gaia DR2 (see Povich et al. 2019 for details).

Prepare the input fitter data file following the guidelines at https://sedfitter.readthedocs.io/en/stable/data.html. This recipe assumes you have a total of `n=10` wavelengths/filters, in this order:

* UKIDSS: *JHK* (filters 1-3)
* 2MASS:  *JHKs* (filters 4-6)
* *Spitzer*/IRAC: [3.6], [4.5], [5.8], [8.0] (filters 7-10)

Other surveys could be substituted (e.g. *WISE* for *Spitzer*/IRAC or VVV for UKIDSS), but the ORDER that photometry points appear in the data file should not be changed. If available, you may also use VVV/UKIDSS *YZ*

**Name the fitter data file `data_xir`.**

## Estimate the maximum reddening to stars in the sample 

Plot and visually inspect the *J – H* vs. *H - K* color-color diagram of the source population:

**IDL>**

   	data = 'data_xir'
	twomass = 0   ; Set to 1 if your JHKs colors are on the 2MASS system
	if not twomass then mk = 1
	magfromdata,data,0,j,nwav=10,mk=mk
	magfromdata,data,1,h,nwav=10,mk=mk
	magfromdata,data,2,k,nwav=10,mk=mk
	plot_nircc_rv,j,h,k,twomass=twomass
	
Estimate the maximum reddening in *Av* magnitudes  by comparing the locus
  of stars to the reddening vectors (marked at intervals of *Av* = 5
  mag). Record the maximum value you estimated. *OPTIONAL: estimate and record a minimum Av if this is not zero.*
  
## FIT WITH DISKLESS PMS MODELS

**Preliminary step:** The two extinction laws that I use are packaged in the `ex_law` subdirectory of this repository. Place a copy of `ex_law_gl.par` into your `$TARGET/sedfitter` directory.

**>>>**

	python
 	from astropy import units as u
  	from sedfitter import fit
  	from sedfitter.extinction import Extinction

	# Define path to models
  	model_dir_pms = <your path here> + '/models_pms'

	# Read in extinction law. 
	extinction = Extinction.from_file('ex_law_gl.par', columns=[0, 3], wav_unit=u.micron, chi_unit=u.cm**2 / u.g)

  	# Define filters and apertures. 
	filters = ['UDSSJ', 'UDSSH', 'UDSSK', '2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4']
  	apertures = [1.5, 1.5, 1.5, 3., 3., 3., 3., 3., 3., 3.] * u.arcsec
	
Apertures are tricky, but the SED `models_pms` set is aperture-independent, so it doesn't really matter here! Note these filter names are *specific* to the pre-convolved model SEDs and must match the filenames in the `models_pms/convolved` directory. The example above is for the combination of UKIDSS, 2MASS, and IRAC filters. Vista filter convolutions are also available, named `VVV[ZYJHK]`.

Make sure the `distance_range` parameter is set to the minimum and maximum distance (in kpc) to your target source population, and the `av_range` reflects the maximum (and minimum if nonzero) extinction estimated from the *JHK* color-color diagram above. Then run the `sedfitter`:

**>>>**
	fit('data_xir', filters, apertures, model_dir_pms, 'xpms.fitinfo', n_data_min=4, extinction_law=extinction, distance_range=[1.5, 2.] * u.kpc, av_range=[0., 25.], output_format=('F',1))

Split the output into well-fit versus poorly-fit. I *strongly recommend* using chi^2/Ndata = 1 as the cutoff for well-fit models from the `models_pms` set, see Povich et al. (2019).

**>>>**
  
  	from sedfitter import filter_output
 	filter_output('xpms.fitinfo', cpd=1.) 

*OPTIONAL: Plot up some SEDs for sanity checks. I personally find plotting the `*bad` SEDs to be more informative than the `*good`.*

**>>>**

  	from sedfitter import plot
  	plot('xpms.fitinfo_good', 'plots_xpms_good') 
  	plot('xpms.fitinfo_bad', 'plots_xpms_bad')

Write SED fit parameters to a text file:

**>>>**

	from sedfitter import write_parameters
  	write_parameters('xpms.fitinfo_good','pars_xpms_good.txt',select_format=('A',-1))

## Make data products from the SED fitting results

Use IDL to make a new fitter data file and SAOImage ds9 regionfile containing only the well-fit sources. Create the `results_xpms` directory contining model fitting results as IDL save files.
  
  **IDL>**
  
  	parfile_good = 'pars_xpms_good.txt' 
  	data_parent='data_xir'
  	fitinfo2data,parfile_good,data_parent,'data_xpms_good'
  	ds9regfromdata,'data_xpms_good','data_xpms_good.reg',color='dodgerblue'

  	target_pms ='results_xpms'
  	pms_pars2idl,parfile_good,target_pms,data_parent=data_parent,filter='F1'	
  
## NEXT STEPS:

  1. The files within the `results_xpms` directory are the starting point for advanced analysis of the PMS model results. Use `phrd_age_analysis_procedure` to create probabilistic HR diagrams plus mass and age distributions for these sources.
  1. Sources that were poorly fit by "naked" PMS models may be YSOs with circumstellar disks/envelopes. *COMING SOON:* a recipe `sedfitting_procedure_xyso` to fit YSO models from Robitaille (2017) to the sources in `xpms.fitinfo_bad `. 


# END OF RECIPE
