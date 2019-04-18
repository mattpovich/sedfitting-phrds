;Create a 2-D probability distribution for a specified combination of 2 PMS star model parameters
;from the SED fitter output. NEW FEATURE IN BETA: For these models one can
;restrict the range of fits used in the distribution based on external
;knowledge of the AGE range for the sources. This is done by
;setting the RESTRICT_AGES keyword, and the constraints are defined
;explicitly in the code (for now).

pro probdist2d_pms, filename, probable, xpbins, ypbins, parx_mean, pary_mean, parx_best, pary_best, parameterx=parameterx,parametery=parametery,parxmin=parxmin,parxmax=parxmax,parymin=parymin,parymax=parymax,nxbin=nxbin,nybin=nybin, restrict_ages=restrict_ages,spec_teff=spec_teff, age_iso=age_iso

;INPUT
;       FILENAME    'string'     -- name of IDL file containing
;                                   parameters (created by FITS2IDL_PMS).
;
;OUTPUT
;       PROBABLE    float[nxbin,nybin] -- 2-D probability distribution of the
;                                    parameter pair. 
;       XBINS,YBINS float[nxbin],float[nybin] -- abcissa values for the probability
;                                   distribution histogram. 
;       PAR[XY]_MEAN      float        -- weighted (by chi2 AND age
;                                        restrictions, if used) parameter values.
;       PAR[XY]_BEST     float        -- parameter values for the best
;                                        fit model (NOTE: These should
;                                        be regarded as LESS RELIABLE
;                                        than PAR[XY]_MEAN!).
;   
;KEYWORDS
;       PARAMETER[XY]=   'string'    -- The YSO fitter parameters to be
;                                   analyzed. DEFAULTS: PARAMETERX='TEFF',PARAMETERY='LBOL'. 
;Accepted parameters are:
;                                     AV: Reddening to star (external parameter)
;                                     DKPC: Distance to star (external parameter)
;                                     AGE:  Stellar evolutionary age (grid parameter)
;                                     MASS: Stellar mass (grid parameter)
;                                     LBOL: Stellar bolometric luminosity (derived parameter)
;                                     REFF: Stellar radius (derived parameter)
;                                     TEFF: Stellar temperature (derived parameter)
;                                     LOGG: Stellar surface gravity (derived parameter)
;
;       PAR[XY]MIN,PAR[XY]MAX  float -- Specifies the minimum and maximum
;                                    values for [XY]BINS. If NOT set,
;                                    then the min and max parameter values from
;                                    the fits will be used.
;       N[XY]BIN=        integer      -- The number of bins to use in
;                                    constructing the probability
;                                    distribution. DEFAULT=20.
;       RESTRICT_AGES= integer    -- Reduce the weight of models from the
;                                    distribution that fall outside
;                                    the mass-dependent age range
;                                    specified in the code.
;                                    IF YSO is NOT set, use integer:
;                                    1 = Use only the older age
;                                    constraints
;                                    2 = Use only the lower age
;                                    constraints 
;                                    3 = Use both constraints
;                                    
;        SPEC_TEFF=  float(2) -- ADVANCED, not prompted at call! [Teff,Teff_err] 2-element vector containing stellar effective temperature and 1-sigma ucertainty. NOTE: SPEC_TEFF OVERRIDES RESTRICT_AGES and AGE_ISO
;        AGE_ISO=    float(2) -- EXPERIMENTAL! Input an assumed isochronal age and uncertainty [age, age_err] in Myr and use this to constrain other parameters.  NOTE: OVERRIDES RESTRICT_AGES
  
;
;WARNINGS
;  The full possible range for "external" model parameters AV and DKPC
;  is not sampled in the SED fitting procedure, hence probability
;  distributions of these parameters should be viewed with suspicion.
;  The "derived" parameters LBOL and TEFF are the "observables" that are
;  most constrained by the data, however given our ASSUMED PMS
;  evolutionary models these are EQUIVALENT to MASS and AGE, the
;  fundamental parameters sampled by the grid. The remaining "derived"
;  parameters REFF and LOGG are similarly not "free" parameters.
;
; ORIGINAL 10 November 2006  -- M.S. Povich
;      30-3-2007 Changed handling of CHI2 probability weighting to avoid
;      floating underflows for high CHI2 values.
;      Adapted PROBDIST to make PROBDIST_PMS 21 June 2011
;      Adapted PROBDIST_PMS to make PROBDIST2D_PMS 23 June 2011
;      Completely renovated functionality of RESTRICT_AGES keyword 18
;      Aug 2011
;      Added PAR[XY]_MEAN output parameters and YSO keyword on 28
;      December 2011
; JTM EDITS: 26 July 2016
;      Adapted from PROBDIST2D_PMS (found in mspovich idl library) to
;      make PROBDIST2D_JTM
;      Added restrict ages if statements to specify constraints for
;      pdisthrd_jtm and mad_hist2
; EHN EDITS: 24 June 2017
;     Adapted from probdist2d_jtm (in jtmaldanado, idl directory) to 
;     make probist2d_ehn
;     Added constraining expression for agecritx and px
;     for pdisthrd_ehn and affiliated scripts
;MSP edits June 29 2017
;  Updated RESTRICT_AGES with YSO keyword set to use age distributions
;  of diskless PMS as posterior constraint on the model fits. EXPERIMENTAL!!
;  Fixed issues with IMPS keyword funcionality.
;  -- July 12 2017 changed logic yet again, to enable X-ray based age
;     constraints for YSOs with X-ray counterparts. Added AGE_DIST keyword
;MSP edits November 2 2018
;  Added SPEC_TEFF keyword to utlize constraints from spectroscopy.

;PRODUCTION VERSION (v1.0)  MSP  -- 15 April 2019
; Removed YSO and AGE_DIST keywords
; Cleaned up code
  
if n_params() lt 4 then begin 
    print,'syntax: probdist2d_pms, filename, probable, xbins, ybins [, parxmean, parymean, parxbest, parybest,' 
    print,'            parameterx=, parametery=, parxmin=, parxmax=, parymin, parymax, nxbin=, nybin=, restrict_ages=, spec_teff=]'
    return
endif

  restore,filename

  if not keyword_set(nxbin) then nxbin = 20
  if not keyword_set(nybin) then nybin = 20
  if not keyword_set(parameterx) then parameterx = 'TEFF'
  if not keyword_set(parametery) then parametery = 'LBOL'
  if parameterx eq parametery then begin
     print,'Hey, it is really boring to plot one parameter against itself! RETURNING.' ;I now see the humor in this since I see what it plots
     return
  endif 

  ;;;;; RESTRICT AGES OF ACCEPTABLE MODELS ;;;;;;;

;High-mass cutoff above which no age-weighting is applied
  himass = 15.
  
;We don't exclude models outright, instead we
;WEIGHT them by probability distribution functions.
  if keyword_set(restrict_ages) and not keyword_set(spec_teff) then begin      
     age = p.age
     mass = p.mass

     if restrict_ages gt 0 then begin
        
        ; Low-age weighting define fixed parameters from the literature.
        mcritdisk = 1.
        agecritdisk = 3.d6
        mcritdiskhi = 10.
        agecritdiskhi = 0.1d6
        alpha = (alog(agecritdisk) - alog(agecritdiskhi)) / (alog(mcritdiskhi) - alog(mcritdisk))
                                ;Set up pdl (relative probability that
                                ;the source has no IR excess from
                                ;inner dust disk
        pdl = replicate(1.,n_elements(age))

       ;Defining agecritx (age where a star's x-ray emission begins to decay due to
           ;the development of a radiative core) as a piecewise array
        agecritx = replicate(1., n_elements(age))
        ind_constant = where (mass le 1)                       ;Low-mass stars
        ind_function = where(mass gt 1 and mass le himass) ;Intermediate-mass stars
        agecritx[ind_constant] = ((1.494)^(2.364))*(10.^6) ;Converting Myr to yr, ~2.6Myr, Gregory et. al 2016
        agecritx[ind_function] = ((1.494/mass[ind_function])^(2.364))*(10.^6) ;Converting Myr to yr, Gregory et. al (2016);

        ;Set up px (relative probability that the source is detected in xrays)
        px = replicate(1., n_elements(age))
        
        ind_hiage_onefive = where(age gt agecritx and mass le 1.5, n_ind_hiage_onefive)
        if n_ind_hiage_onefive ne 0 then $
           px[ind_hiage_onefive] = (age[ind_hiage_onefive]/agecritx[ind_hiage_onefive])^(-0.75) ;Gregory et. al (2016)

           ;px for 1.5 < mass <= 2 Msun
        ind_hiage_onefive_two = where(age gt agecritx and mass gt 1.5 and mass le 2., n_ind_hiage_onefive_two)
        if n_ind_hiage_onefive_two ne 0 then $
           px[ind_hiage_onefive_two] = (age[ind_hiage_onefive_two]/agecritx[ind_hiage_onefive_two])^(-0.86)  

           ;px for 2 < mass
        ind_hiage_two_eight = where(age gt agecritx and mass gt 2. and mass le himass, n_ind_hiage_two_eight)
        if n_ind_hiage_two_eight ne 0 then $
           px[ind_hiage_two_eight] = (age[ind_hiage_two_eight]/agecritx[ind_hiage_two_eight])^(-1.19)

        ind_lomass_loage = where(mass lt mcritdisk and age lt agecritdisk,nlomass)
        if nlomass ne 0 then begin
           agecritsig = agecritdisk/5. ;Empirical! 
           pdl[ind_lomass_loage] = 0.5*(1 + erf((age[ind_lomass_loage] - agecritdisk/2.)/(sqrt(2)*agecritsig))) 
        endif
;         Intermediate-mass weighting (NO high-mass weighting)
        ind_im = where(mass ge mcritdisk and mass le 10. and age lt agecritdisk*(1./mcritdisk)*mass^(-1*alpha),nim)
        if nim ne 0 then begin
           agecritm = agecritdisk*(mass[ind_im]/mcritdisk)^(-1.*alpha)
           agesigm = agecritm/5.
           pdl[ind_im] = 0.5*(1 + erf((age[ind_im] - agecritm/2.)/(sqrt(2)*agesigm))) 
        endif
                  ;End lower age constraint

     endif else begin ;Use previously-computed age distribution from DISKLESS source population instead
        if not keyword_set(age_dist) then begin
           age_dist = ''
           read,age_dist,prompt='Please specify a path to the IDL save file containing the age distribution variables (AGE_BINS, AGE_HIST): '
        endif 
        if not file_test(strtrim(age_dist,2)) then begin
           print,'WARNING: The file requested by RESTRICT_AGES, '+strtrim(restrict_ages,2)+' does not exist! RETURNING.'
           return
        endif
        restore,age_dist
        pyso = replicate(1.,n_elements(age))
        pyso_all = interpol(age_hist,age_bins*1.d6,age)
        peak = max(pyso_all,peak_loc)
        agecrityso = age[peak_loc]
        ind_hiageyso = where(age ge agecrityso,n_hiageyso)
        if n_hiageyso ne 0 then pyso[ind_hiageyso] = pyso_all[ind_hiageyso]/peak
        
     endelse 
  endif  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; Select model weights using Teff measured from spectroscopy
  if keyword_set(spec_teff) then begin
     ;idiot-proofing..
     if n_elements(spec_teff) ne 2 then begin
        print,'SPEC_TEFF MUST be a 2-element vector = [TEFF,TEFF_ERR]! RETURNING.'
        return
     endif  
     par_t = p.teff

;Initialize pt (relative probability that the source has a given teff)
     pt = gaussian( par_t, [1, spec_teff[0,0], spec_teff[0,1] ] )
     
  endif ;SPEC_TEFF restrictions

;; Select model weights using isochronal age and uncertainty
  if keyword_set(age_iso) then begin
     ;idiot-proofing..
     if n_elements(age_iso) ne 2 then begin
        print,'AGE_ISO MUST be a 2-element vector = [AGE,AGE_ERR]! RETURNING.'
        return
     endif  
     par_a = p.age
     
;Initialize pa (relative probability that the source has a given age)
     pa = gaussian( par_a, [1, 1.e6*age_iso[0], 1.e6*age_iso[1] ] )
     
  endif ;End AGE_ISO restrictions

  
;Find BEST-FIT model for comparison
  minchi2 = min(p.chi2,ind_best)

  case 1 of  ;X AXIS selection
      parameterx eq 'MASS': begin
         parx = p.mass
         logxbins = 1
         parx_best = p[ind_best].mass
      end
      parameterx eq 'AGE': begin
         parx = p.age
         logxbins = 1
         parx_best = p[ind_best].age
      end
      parameterx eq 'LBOL': begin
         parx = p.lbol 
         logxbins = 1
         parx_best = p[ind_best].lbol
      end 
      parameterx eq 'REFF': begin
         parx = p.reff
         logxbins = 0
         parx_best = p[ind_best].reff
      end 
      parameterx eq 'TEFF': begin
         parx = p.teff
         logxbins = 1
         parx_best = p[ind_best].teff
      end 
      parameterx eq 'LOGG': begin
         parx = p.reff
         logxbins = 0
         parx_best = p[ind_best].logg
      end 
      parameterx eq 'AV': begin
         parx = p.av
         logxbins = 0
         parx_best = p[ind_best].av
      end 
      parameterx eq 'DKPC': begin
         parx = 10^(p.logd)
         logxbins = 0
         parx_best = 10^(p[ind_best].logd)
      end 
      else: begin
         print,'Please specify a valid PMS parameterX for examination.'
         print,' Options are MASS, AGE, LBOL, TEFF, REFF, LOGG. RETURNING.'
         return
      end 
   endcase
  
  case 1 of  ;Y AXIS selection
      parametery eq 'MASS': begin
         pary = p.mass
         logybins = 1
         pary_best = p[ind_best].mass
      end
      parametery eq 'AGE': begin
         pary = p.age
         logybins = 1
         pary_best = p[ind_best].age
      end
      parametery eq 'LBOL': begin
         pary = p.lbol 
         logybins = 1
         pary_best = p[ind_best].lbol
      end 
      parametery eq 'REFF': begin
         pary = p.reff
         logybins = 0
         pary_best = p[ind_best].reff
      end 
      parametery eq 'TEFF': begin
         pary = p.teff
         logybins = 1
         pary_best = p[ind_best].teff
      end 
      parametery eq 'LOGG': begin
         pary = p.reff
         logybins = 0
         pary_best = p[ind_best].logg
      end 
      parametery eq 'AV': begin
         pary = p.av
         logybins = 0
         pary_best = p[ind_best].av
      end 
      parametery eq 'DKPC': begin
         pary = 10^(p.logd)
         logybins = 0
         pary_best = 10^(p[ind_best].logd)
      end 
      else: begin
         print,'Please specify a valid PMS parameterY for examination.'
         print,' Options are MASS, AGE, LBOL, TEFF, REFF, LOGG. RETURNING.'
         return
      end 
   endcase

  if not keyword_set(parxmin) then parxmin = min(parx)
  if not keyword_set(parxmax) then parxmax = max(parx)
  if not keyword_set(parymin) then parymin = min(pary)
  if not keyword_set(parymax) then parymax = max(pary)

  x0 = double(parxmin)
  xn = double(parxmax)
  if logxbins eq 1 then begin
     dlogpx = (alog10(xn)-alog10(x0))/nxbin
     xbins = x0*10^(dlogpx*findgen(nxbin+1))
  endif else begin
     dpx = (xn - x0)/nxbin
     xbins = x0 + findgen(nxbin+1)*dpx
  endelse
  y0 = double(parymin)
  yn = double(parymax)
  if logybins eq 1 then begin
     dlogpy = (alog10(yn)-alog10(y0))/nybin
     ybins = y0*10^(dlogpy*findgen(nybin+1))
  endif else begin
     dpy = (yn - y0)/nybin
     ybins = y0 + findgen(nybin+1)*dpy
  endelse

  
  ;make the 2-D probability histogram
 ;Compute probability from Chi2
  chi2 = double(p.chi2)  ;It's scary how little this seems to matter...
  probs = double(exp(-chi2/2.))

  ;APPLY AGE-BASED WEIGHTING FUNCTIONS ;;;;;;;;;
  if keyword_set(restrict_ages) $
     and not keyword_set(spec_teff) then begin
     case 1 of
        restrict_ages eq 1: probs = probs*px ;only old-age weighting from X-rays
        restrict_ages eq 2: probs = probs*pdl ;only young-age weighting from diskless
        restrict_ages eq 3: probs = probs*px*pdl ;Combine all weighting functions - probability of age from x-ray and dl from diskless stars
        restrict_ages eq -1: probs = probs*pyso
        else: begin
           print,'RESTRICT_AGES must have a value of -1, 1, 2, or 3. RETURNING.'
           return
        end
     endcase
  endif
  ;;;;;;;;;;;;;;

  
;APPLY SPEC_TEFF GAUSSIAN WEIGHTING FUNCTION
  if keyword_set(spec_teff) then begin
     if max(probs*pt) eq 0 then begin
        print,'WARNING: SPEC_TEFF weighting found zero probability! IGNORING for source '+ s.desig
     endif else probs = probs*pt
  endif

    ;APPLY AGE_ISO GAUSSIAN WEIGHTING FUNCTION;;
  if keyword_set(age_iso) then begin
     if max(probs*pa) eq 0 then begin
        print,'WARNING: AGE_ISO weighting found zero probability! IGNORING for source '+ s.desig
     endif else probs = probs*pa
  endif 
  
  
  ;Compute weighted mean parameter values
  if logxbins eq 1 then $
     parx_mean = 10.^(total(probs*alog10(parx))/total(probs)) else $
        parx_mean = total(probs*parx)/total(probs)
  if logybins eq 1 then $
     pary_mean = 10.^(total(probs*alog10(pary))/total(probs)) else $
        pary_mean = total(probs*pary)/total(probs)


  probable = dblarr(nxbin,nybin) ;2 dimensions, stupid!

  xpbins = fltarr(nxbin)
  ypbins = fltarr(nybin)
  for i=0,nxbin-1 do begin      ;Loop over X bins
     xpbins[i] = 0.5*(xbins[i]+xbins[i+1])
     whrx = where(parx ge xbins[i] and parx lt xbins[i+1],nwhrx)
     for j=0,nybin-1 do begin                 ;Loop over Y bins
        ypbins[j] = 0.5*(ybins[j]+ybins[j+1]) ;This is NOT efficient or clever.           
        if nwhrx ne 0 then begin
           whry = where(pary[whrx] ge ybins[j] and pary[whrx] lt ybins[j+1],nwhr)
           if nwhr ne 0 then probable[i,j] = total(probs[whrx[whry]]) else probable[i,j] = 0.
        endif else probable[i,j] = 0.
     endfor                     ;j
  endfor                        ;i
  
  probable = probable/total(probable) ;Normalize!     

end
