pro mad_hist2_pms, target_dir, agebins, massbins, age_hist, mass_hist, age, age_err, nabin=nabin, nmbin=nmbin, sourcelist=sourcelist, percentile=percentile, offset=offset, age_offset=age_offset, restrict_ages=restrict_ages,region_name=region_name,yso=yso, spec_teff_file=spec_teff_file, age_iso=age_iso, input_mass=input_mass, input_age=input_age, noplot=noplot

; INPUT
;       TARGET_DIR     'string'       -- Directory containg source files
;
; OUTPUT/INPUT
;       AGEBINS        float[nabin]   -- abcissa values for age histogram
;       MASSBINS       float[nmbin]   -- abcissa values for age
;                                        histogram
;       AGEHIST        float[nabin]   -- age histogram
;       MASSHIST       float[nmbin]   -- mass histogram
;
;OUTPUT
;       AGE            float          -- age at cutoff 
;       AGE_ERR        float          -- uncerainty on AGE (based on
;                                        histogram binsize)
;                                       
;
; KEYWORDS
;       NABIN,NMBIN=    integer        -- binsize for each axis
;       
;       SOURCELIST=     'string'       -- name of the text file under the TARGET_DIR directory containing a list of sources
;
;       RESTRICT_AGES= (integer OR 'string')    -- Reduce the weight of models from the
;                                    distribution that fall outside
;                                    the mass-dependent age range
;                                    specified in the code.
;                                    IF YSO is NOT set, use integer:
;                                    1 = Use only the older age
;                                    constraints
;                                    2 = Use only the lower age
;                                    constraints 
;                                    3 = Use both constraints
;                                    IF YSO is set, this keyword
;                                    should be a string giving the
;                                    path to an IDL save file
;                                    containing the desired weighting
;                                    distribution.
;                
;       OFFSET=      integer   -- If set, selects the histogram bin
;                                 OFFSET places to the right of the
;                                 histogram peak for the mass cutoff.
;       AGE_OFFSET=      integer   -- If set, selects the AGE histogram bin
;                                 OFFSET places to the LEFT of the
;                                 histogram peak for the AGE determination.
;      AGE_ISO=      float(2)  -- [Age,uncertainty] in Myr used to
;                                 constraint SED model parameters
;                                 based on a single, isochronal age.
;      /INPUT_MASS,/INPUT_AGE (switch) -- If set, use the input mass and/or age
;                             histograms for plotting, rather than
;                             calculating these internally.
;
;;
; CALLS: PDISTHRD_PMS
;
;
;VERSION HISTORY
;ORIGINAL  -- JTM Summer 2016
;EHN NOTE:28 June 2017
;Added switch for inset plots and updated comments  
;on called scripts
;MSP added /YSO keyword functionality June 30, 2017
;   MSP improved plotting functionality October 2018
;   MSP added ability to iterate using pre-determined mass/age
;   distributions (AGEBINS, MASSBINS, AGEHIST, MASSHIST input/output parameters), and removed the deprecated /INSET_PLOT functionality

;PRODUCTION VERSION (v1.0) MSP -- 15 April 2019
  ;Removed YSO, PERCENTILE keywords
  ;Cleaned up code and comments.
  
  if n_params() lt 1 then begin
     print,'syntax -- mad_hist2_pms, target_dir[, agebins, massbins, agehist, masshist, age,'
     print,'  nabin=, nmbin=, psplots=, sourcelist=, offset=, age_offset=,'
     print,'  restrict_ages=,region_name=]'
     return
  endif 
  
;Define some keywords up front
  if not keyword_set(sourcelist) then sourcelist = 'sourcelist.txt'
  if not keyword_set(offset) then offset = 0
  if not keyword_set(age_offset) then age_offset = 0
  if not keyword_set(region_name) then region_name = sourcelist.Remove(-4)
  if keyword_set(input_mass) then nmbin = n_elements(massbins)
  if keyword_set(input_age) then nabin = n_elements(agebins)


;Optimize the histogram binning

  nsrc = file_lines(target_dir+'/'+sourcelist)
  if not keyword_set(nxbin) or not keyword_set(nybin) then begin
     case 1 of
        nsrc gt 250: begin
           nabin = 40/(1+keyword_set(offset))
           nmbin = 50
        end
        nsrc gt 100 and nsrc le 250: begin
              nabin = 30/(1+keyword_set(offset))
              nmbin = 40
        end
        nsrc le 100: begin 
           nabin = 20/(1+keyword_set(offset))
           nmbin = 30
        end
     endcase
  endif 

  
  pdisthrd_pms, target_dir, cumhrd, agebins, massbins, nxbin=nabin, nybin=nmbin, sourcelist=sourcelist, /mad, restrict_ages=restrict_ages, spec_teff_file=spec_teff_file, age_iso=age_iso

; Make histogram of stellar mass   
  if not keyword_set(input_mass) then mass_hist = total(cumhrd,1)

; Find maximum value of the number of stars per mass bin & over plot Salpeter relation for a loglog plot
   maxn_mass = max(mass_hist,ind_max)    
   ind_use = ind_max + offset            
   max_mass = massbins[ind_use]
   salpeter = mass_hist[ind_use]/(massbins[ind_use]^(-1.35))*massbins^(-1.35)    

   ; Define histogram of stellar age, omitting high-mass stars
  high_mass = where(massbins GT 8.) 
  if not keyword_set(input_age) then age_hist = total(cumhrd[*,ind_use:high_mass[0]],2)
 
  numstars = total(age_hist)
  numstarsmass = total(mass_hist[ind_use:high_mass[0]])
  
; Compute stellar deficit in X-ray desert
   frac_desert = float(numstarsmass)/float(total(salpeter[ind_use:high_mass[0]]))

;Print some key metrics   
   print, 'Number of Stars in Region',fix(total(mass_hist))
   print, 'Mass Cutoff', max_mass
   print, 'Intermediate-Mass Stars in AGE distribution',numstars
   print, 'Intermediate-Mass Stars in MASS distribution',numstarsmass
   print, 'Intermediate-Mass Stars in Scaled IMF',total(salpeter[ind_use:high_mass[0]])
   print, 'Fraction of X-ray selected IMPS and AB stars',frac_desert
   maxn_age = max(age_hist,index,/nan)
   index = index - age_offset
   age = sqrt(agebins[index]*agebins[index+1])
   age_err = (sqrt(agebins[index]*agebins[index+1]) - agebins[index])*(age_offset + 1)
   print, 'Age at break point',age,'+/-',age_err
     
   if age lt 10. then age_format = '(F3.1)' else age_format = '(F4.1)'

   if not keyword_set(noplot) then begin
      srclen = strlen(sourcelist)
      sourceliststub = strmid(sourcelist,0,srclen-4)
      if not keyword_set(age_iso) then $
         cgps_open,sourceliststub+'_MassAgeHist_'+string(max_mass,format='(F3.1)')+'-'+strtrim(offset,2)+'_'+string(age,format=age_format)+'.eps',/encap,/nomatch else $
            cgps_open,sourceliststub+'_MassAgeHist_'+string(max_mass,format='(F3.1)')+'-'+strtrim(offset,2)+'_'+string(age,format=age_format)+'_ISO.eps',/encap,/nomatch 
         
      device,xsize=5,ysize=8,/inches
      !p.multi = [0,1,2]
      chsize = 1
      th = 7
      xmarg = [7,2]

      cgplot, agebins, age_hist, ps=10,xtitle='Age (Myr)',ytitle='Number of Stars',/xlog,xr=[0.1,30], thick=th, charthick=th, xcharsize=chsize, ycharsize=chsize, xthick=th, ythick=th, /xsty, xmargin=xmarg,ymargin=[4,1]
      cgoplot,[replicate(age,2)],[0.01, 1000], linestyle=3, thick=th
      cgoplot,[replicate(age-age_err,2)],[0.01, 1000], linestyle=1, thick=th
      cgoplot,[replicate(age+age_err,2)],[0.01, 1000], linestyle=1, thick=th
      if age gt 1 then $
         cgtext, age - 2*age_err, 0.9*!y.crange[1],string(age,format='(F4.1)')+cgsymbol('+-')+string(age_err,format='(F3.1)')+' Myr ', charsize=1.5*chsize, charthick=th, /align $
      else $
         cgtext, age + 2*age_err, 0.9*!y.crange[1],string(age,format='(F4.1)')+cgsymbol('+-')+string(age_err,format='(F3.1)')+' Myr ', charsize=1.5*chsize, charthick=th

      cgplot, massbins, mass_hist, ps=10,xtitle='Mass (M'+sunsymbol()+')',ytitle='Number of Stars',/ylog,/xlog, yrange=[1.0,salpeter[ind_max]], xrange=[1.0,20.0], thick=th, charthick=th, xcharsize=chsize, ycharsize=chsize, xthick=th, ythick=th, /xsty, xmargin = xmarg, ymargin=[4,1]
      cgoplot, massbins,salpeter, th=th, /NoErase
      cgoplot,[replicate(massbins[ind_use],2)],[0.01, 1000], linestyle=3, thick=th, /NoErase
      cgtext, massbins[ind_use+1], 1.2, string(max_mass,format='(F0.1)')+' M'+sunsymbol(), charsize=2, charthick=th
      if keyword_set(region_name) then cgtext, 0.9, 0.4, region_name, charthick=th, charsize=1.2*chsize,/norm, /align
      
      cgps_close  
   endif
   
end
