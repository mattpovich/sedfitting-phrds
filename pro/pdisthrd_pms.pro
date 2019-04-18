;Make an HR diagram (logLBOL versus logTEFF) from probability distributions of SED fitting parameter results for an ensemble of diskless PMS stars. NOTE: This routine works with the 
;models_pms grid parameters ONLY! 

pro pdisthrd_pms, target_dir, cumhrd, xbins, ybins, nxbin=nxbin, nybin=nybin, psplots=psplots, mad=mad, sourcelist=sourcelist, restrict_ages=restrict_ages, spec_teff_file=spec_teff_file, age_iso=age_iso, oplot_mean=oplot_mean,oplot_best=oplot_best,IMPS=IMPS,region_name=region_name

;INPUT
;      TARGET_DIR   'string'     -- Directory containing IDL save
;                                   files of fit parameters and
;                                   sourcelist produced by
;                                   FITS2IDL_PMS. Example: 'results_newstars_good'
;
;OUTPUT
;      CUMHRD       float[nxbin,nybin] -- The HR diagram as a 2-D histogram of model parameters.
;      XBINS,YBINS float[nxbin],    -- The abcissa and ordinate values for HRD.
;
;KEYWORDS
;      NXBIN=,NYBIN=     integer    -- The number of bins to use in constructing each axis of the 2D histogram.
;
;      PSPLOTS=      integer      -- =1: Make PostScript plot of the probability distribution function.
;                                    =2: Make PostScript plots of the probability distributions for EACH source in SOURCELIST.
;                                    >2: Make all of the PostScript plots above, saves 2 separate PS files.
;      /MAD          (switch) -- If set, plots a MASS-AGE Diagram instead of the traditional (DEFAULT) HR Diagram.
;      SOURCELIST=   'string' -- The name of a text file containing a list of
;                                IDL save files, 1 for each source in
;                                TARGET_DIR to be
;                                analyzed. DEFAULT =
;                                target_dir+'sourcelist.txt'
;      MODELS_DIR=   'string' -- <NOT USED> Allows manual specification of the
;                                 path to the models_pms
;                                 directory. DEFAULT is to let
;                                 MODELS_DIR = TARGET_DIR.
;       RESTRICT_AGES= integer    -- Reduce the weight of models from the
;                                    distribution that fall outside
;                                    the mass-dependent age range
;                                    specified in the code.
;                                    1 = Use only the older age
;                                    constraints
;                                    2 = Use only the lower age
;                                    constraints 
;                                    3 = Use both constraints
;                                    
;      SPEC_TEFF_FILE= 'string'   -- Path to an IDL aave file
;                                    containing a floating-point
;                                    arrage SPEC_TEFF(nsrc,2) giving
;                                    the temperature and uncertainty
;                                    for each source in SOURCELIST
;      AGE_ISO=      float(2)    -- Estimated isochronal age and
;                                   uncertainty for the ensemble
;                                   distribution in SOURCELIST. Passed
;                                   along to PROBDIST2D_PMS. 

;      /OPLOT_MEAN    (switch) -- If set, the weighted average
;                                parameter values for each source are
;                                 overplotted as WHITE dots on the output HRDs
;                                 (note, if PSPLOTS is NOT set,
;                                 /OPLOT_BEST does nothing).
;      /OPLOT_BEST    (switch) -- If set, the parameters for the
;                                 BEST-FIT model to each source are
;                                 overplotted as RED dots on the output HRDs
;                                 (note, if PSPLOTS is NOT set,
;                                 /OPLOT_BEST does nothing).
;       /IMPS         (switch)    -- Classify objects as IMPS, AB, OB,
;                                    or T Tauri, and make new sourcelists.
;
;
;CALLS: PROBDIST2D_PMS, OPLOT_SIESS, OPLOT_SIESS_ZAMS, OPLOT_BM_ZAMS
;       CGPS_OPEN, CGPS_CLOSE, CGCONTOUR, CGOPLOT
;
;MSP wrote 23 June 2011
;    BASED ON PDISTFUNC_PMS.PRO
;    Finished 11 July, 2011 after paternity leave. :-)
;          Added OPLOT_MEAN, OPLOT_BEST, YSO keywords 28 December 2011
;          Added IMPS keyword 9 May 2012
;JTM edits June-July 2016
;  Improved contour scaling and overplots for more visually-appealing output.
;MSP edits June 29 2017
;  Removed MODELS_DIR keyword as it inconveniently required one to always need the
;  model parameter.fits files handy, for the sole purpose of computing
;  minmax parameter values! THESE ARE NOW HARD-CODED because for the
;  time being there are only two grids that we use.
;  Updated RESTRICT_AGES with YSO keyword set to input age
;  distributions passed along to PROBDIST2D_PMS. EXPERIMENTAL!!
;  Fixed issues with IMPS keyword functionality...
;MSP edits November 2 2018
;  Added SPEC_TEFF_FILE keyword functionality
;
;PRODUCTION version 1.0    MSP -- 15 April 2019
;   * Adopted Coyote Graphics routine for (hopefully) better handling of
;        color tables across platforms
;   * Removed YSO, AGE_DIST, and FORTRAN keywords  
;   * Edited default N[XY]BINS to be chosen pased on number of sources
;    in SOURCELIST
;   * Changed some variable names to make more general
;   * Removed unnecessary comments
  
if n_params() lt 1 then begin 
    print,'syntax: pdisthrd_pms, target_dir, hrd, xbins, ybins, '
    print, '         nxbin=, nybin=, psplots=, /mad, sourcelist=, restrict=,/oplot_mean,/oplot_best,/imps, spec_teff_file='
    return
endif

  ;get source names
  if not keyword_set(sourcelist) then sourcelist = 'sourcelist.txt'
  sourcelistuse=target_dir + '/' + sourcelist
  readcol,sourcelistuse,format='A',names
  files = target_dir+'/'+names

  nsrc = n_elements(files)

  if not keyword_set(nxbin) or not keyword_set(nybin) then begin
     case 1 of
        nsrc gt 250: begin
           if not keyword_set(mad) then begin
              nxbin = 30
              nybin = 30
           endif else begin
              nxbin = 40
              nybin = 50
           endelse 
        end
        nsrc gt 100 and nsrc le 250: begin
           if not keyword_set(mad) then begin
              nxbin = 20
              nybin = 20
           endif else begin
              nxbin = 30
              nybin = 40
           endelse 
        end
        nsrc le 100: begin 
           if not keyword_set(mad) then begin
              nxbin = 15
              nybin = 15
           endif else begin
              nxbin = 20
              nybin = 30
           endelse 
           print,'WARNING! <100 sources detected in '+strtrim(sourcelist,2)+'.' 
           print,'   Results may be unreliable.'
        end
     endcase
  endif 

  hrd = dblarr(nxbin,nybin,nsrc)

;If called for, load and verify arrays for SPEC_TEFF weighting
  if keyword_set(spec_teff_file) then begin
     restore,spec_teff_file
     if not keyword_set(spec_teff_arr) then begin
        print,spec_teff_file + ' MUST contain a 2-D floating-point array named SPEC_TEFF_AR. RETURNING.'
        return
     endif
     sztarr = size(spec_teff_arr)
     if sztarr[0] ne 2 or sztarr[1] ne nsrc or sztarr[2] ne 2 then begin
        print, spec_teff_file + ' MUST contain a 2-D floating-point array named SPEC_TEFF_ARR.'
        print,'  Columns are TEFF, TEFF_ERR and rows must match '+sourcelistuse
        print,'       RETURNING.'
        return
     endif  
  endif  
  
;Set up plotting of probability distribution for EACH source, if desired
  if keyword_set(psplots) then begin
     loadct,39                  ;Rainbow + white is nice
     psplots = fix(psplots) 
  endif else psplots = 0

  srclen = strlen(sourcelist)
  sourceliststub = strmid(sourcelist,0,srclen-4) 

  if psplots ge 2 then begin
     !p.multi=[0,3,4]
     if not keyword_set(mad) then plotname = 'individual_HRDs' else plotname = 'individual_MADs'
     if keyword_set(restrict_ages) then plotname = plotname + '_agecons'
     if keyword_set(age_iso) then plotname = plotname + string(age_iso[0],format='(F3.1)')+'Myr'
;     ps_open, sourceliststub + '_' + plotname,/port,/co 
     cgps_open, sourceliststub + '_' + plotname+'.ps',/nomatch
     device,xsize=6.5,ysize=9,/inches
  endif

;Calculate the probability distributions for EACH source.
  if not keyword_set(mad) then begin
     parxname = 'TEFF'   
     paryname = 'LBOL'   
     parymin = 0.0157300
     parymax = 743950.
     parxmin = 2535.60
     parxmax = 45970.0
     parxr = [4.6,3.55] ;NOTE: EXPLICITLY setting temperature range for plotting!
     paryr = [-1,4.5]  ;NOTE: EXPLICITLY setting luminosity range for plotting!
  endif else begin
     parxname = 'AGE'
     paryname = 'MASS'
     parymin = 0.1
     parymax = 49.9960
     parxmin = 1001.90
     parxmax = 1.5e+07
     parxr = [0.01,10.]       ;NOTE: EXPLICITLY setting age range for plotting!
     paryr = [0.15,20.]  ;NOTE: EXPLICITLY setting mass range for plotting!
  endelse
  
  ;Arrays to hold weighted mean parameter values
  if keyword_set(oplot_mean) then begin  
     xmeanarr = fltarr(nsrc)
     ymeanarr = fltarr(nsrc)
  endif
  ;Arrays to hold best-fit parameter values
  if keyword_set(oplot_best) then begin  
     xbestarr = fltarr(nsrc)
     ybestarr = fltarr(nsrc)
  endif 

  ;Find IMPS -- setup
  if keyword_set(IMPS) then begin
     class = [ 'TTauri', 'OB', 'IMPS', 'AB', 'unc.' ]
     nclass = n_elements(class)
     uclass = indgen(nclass)+1
     count = intarr(nclass) 
     for jcl=0,nclass-1 do openw,uclass[jcl],sourceliststub + '_' + class[jcl] + '.txt'
  endif 

  for isrc=0L, nsrc-1 do begin
     if keyword_set(spec_teff_arr) then spec_teff = spec_teff_arr[isrc,*]
     probdist2d_pms,files[isrc],probhrd,xt,yl,tmean,lmean,tbest,lbest,restrict_ages=restrict_ages,parxmin=parxmin,parxmax=parxmax,parymin=parymin,parymax=parymax,nxbin=nxbin,nybin=nybin,parameterx=parxname,parametery=paryname,spec_teff=spec_teff,age_iso=age_iso
     hrd[*,*,isrc] = probhrd

     if keyword_set(imps) then begin
        probdist2d_pms,files[isrc],tmdist,tt,mm,tmean,mmean,restrict_ages=restrict_ages,parameterx='TEFF',parametery='MASS',parxmax=47965.319,parymax=59.9680,yso=keyword_set(yso),nxbin=2*nxbin,nybin=2*nybin, spec_teff=spec_teff, age_iso=age_iso 
        
       ;Classify source using IMPS scheme
        teffcool = 7300        ;K, subject to change, but this is for F0 star..
        massim = 2.
        ind_cool = where(tt lt teffcool,n_cool,complement=ind_hot,ncomplement=n_hot)
        ind_lomass = where(mm lt massim,n_lomass)
        ind_immass = where(mm ge massim and mm lt 8.,n_immass)
        ind_himass = where(mm ge 8.,n_himass)
         
       ;NEW: properly carve up probability distributions
        Prel_TTauri = total(tmdist[ind_cool[0]:ind_cool[n_cool-1],ind_lomass[0]:ind_lomass[n_lomass-1]])
        Prel_OB = total(tmdist[*,ind_himass[0]:ind_himass[n_himass-1]])
        Prel_IMPS = total(tmdist[ind_cool[0]:ind_cool[n_cool-1],ind_immass[0]:ind_immass[n_immass-1]])
         Prel_AB = total(tmdist[ind_hot[0]:ind_hot[n_hot-1],ind_immass[0]:ind_immass[n_immass-1]])

                                ;Emulating the decision rule from
                                ;Broos et al. 2011 for CCCP;
                                ;classification if Prel(max) >
                                ;Prel(next)
         Prel = [Prel_OB,Prel_AB,Prel_IMPS,Prel_TTauri]
         psrt = sort(Prel)
         if Prel[psrt[3]] gt 2*Prel[psrt[2]] then begin
            case 1 of
               psrt[3] eq 3: jsel = 0 ;T Tauri
               psrt[3] eq 0: jsel = 1 ;OB
               psrt[3] eq 2: jsel = 2 ;IMPS
               psrt[3] eq 1: jsel = 3 ;AB
            endcase
         endif else jsel = 4    ;unclassified
         
         printf,uclass[jsel],names[isrc]  
         count[jsel]++   
      endif ;Done finding IMPS

      if keyword_set(oplot_mean) then begin
         if not keyword_set(mad) then begin
            xmeanarr[isrc] = alog10(tmean)
            ymeanarr[isrc] = alog10(lmean)
         endif else begin
            xmeanarr[isrc] = tmean/1.e6
            ymeanarr[isrc] = lmean            
         endelse 
      endif 
      if keyword_set(oplot_best) then begin
         if not keyword_set(mad) then begin
            xbestarr[isrc] = alog10(tbest)
            ybestarr[isrc] = alog10(lbest)
         endif else begin
            xbestarr[isrc] = tbest/1.e6
            ybestarr[isrc] = lbest            
         endelse 
      endif 

      if psplots ge 2 then begin
         ;Get optimal contour levels
         iHRDsrt = probhrd[sort(probhrd)]
         iHRDsrt = iHRDsrt[where(iHRDsrt gt 0)]
         ipdfHRD = total(iHRDsrt,/cumulative)/total(iHRDsrt)
         itoplev = 1.0           ; top contour
         ibotlev = 0.05          ; bottom contour
         iwhrtop = where(ipdfHRD ge itoplev,intop)
         if intop ne 0 then ilevmax = iHRDsrt[iwhrtop[0]] else ilevmax = itoplev
         iwhrbot = where(ipdfHRD le ibotlev,inbot)
         if inbot ne 0 then ilevmin = iHRDsrt[iwhrbot[inbot-1]] else ilevmin = ibotlev
; Linear + 0 Contour Levels...     
         numlevs = 30
         idlogc = ((ilevmax)-(ilevmin))/numlevs
         ilevs = [0, ilevmin+(idlogc*findgen(numlevs+1))]
         if not keyword_set(mad) then begin
            cgcontour,probhrd,alog10(xt),alog10(yl),xtitle='log(T!Deff!N) [K]',ytitle='log(L!Dbol!N) [L'+sunsymbol()+']',title=strmid(names[isrc],0,strlen(names[isrc])-4),levels=ilevs,/fill,xr=[parxr[0],3.45],/xsty,/ysty,xthick=6,ythick=6,charth=2,charsize=0.8
         ; Overplot SDF00 tracks and isochrones
            oplot_siess,/tracks,color='black', thick=4, /annotate, charsize=0.5, charthick=2
            oplot_siess,/tracks,color='white', thick=3, /annotate, charsize=0.5, charthick=2

            cgaxis,!x.crange[0],!y.crange[0],xr=xr,/xsty,color='white',xtickformat='(A1)',xth=6 
            cgaxis,!x.crange[0],!y.crange[1],xr=xr,/xsty,color='white',xtickformat='(A1)',xth=6,xaxis=1
            cgaxis,!x.crange[0],!y.crange[0],yr=yr,/ysty,color='white',ytickformat='(A1)',yth=6,yaxis=0
            cgaxis,!x.crange[1],!y.crange[0],yr=yr,/ysty,color='white',ytickformat='(A1)',yth=6,yaxis=1
            if keyword_set(oplot_best) then begin
               filledcirc
               cgoplot,[tbest],[lbest],ps=8,symsize=0.6,co=254
               opencirc
               cgoplot,[tbest],[lbest],ps=8,symsize=0.6
            endif 
            if keyword_set(oplot_mean) then begin
               filledcirc
               cgoplot,[alog10(tmean)],[alog10(lmean)],ps=8,symsize=0.6,co=255
               opencirc
               cgoplot,[alog10(tmean)],[alog10(lmean)],ps=8,symsize=0.6
            endif 
         endif else begin
            cgcontour,probhrd,xt/1.e6,yl,xtitle='Age [Myr]',ytitle='Mass [M'+sunsymbol()+']',title=strmid(names[isrc],0,strlen(names[isrc])-4),levels=ilevs,/fill,/xsty,/ysty,xthick=6,ythick=6,charth=2,/xlog,/ylog,xr=[parxr],yr=paryr,charsize=0.8
            cgaxis,parxr[0],paryr[0],/xsty,color='white',xtickformat='(A1)',xth=6,/xlog 
            cgaxis,parxr[0],paryr[1],/xsty,color='white',xtickformat='(A1)',xth=6,xaxis=1,/xlog
            cgaxis,parxr[0],paryr[0],/ysty,color='white',ytickformat='(A1)',yth=6,yaxis=0,/ylog
            cgaxis,parxr[1],paryr[0],/ysty,color='white',ytickformat='(A1)',yth=6,yaxis=1,/ylog
;Overplot ZAMS
            oplot_siess_zams,color=0,thick=6
            oplot_bm_zams,color=0,thick=5
            oplot_siess_zams,color=255,thick=4
            oplot_bm_zams,color=255,thick=3
;Overplot disk boundary lines (THESE ARE SUBJECT TO CHANGE)
            oplot,replicate(3.,2),[0.1,1],co=0,th=6,li=1
            oplot,[3.,0.1],[1.,10.],co=0,th=6,li=1
            oplot,replicate(3.,2),[0.1,1],co=255,th=4,li=1
            oplot,[3.,0.1],[1.,10.],co=255,th=4,li=1
            
            if keyword_set(oplot_best) then begin
               filledcirc
               cgoplot,[tbest/1.e6],[lbest],ps=8,symsize=0.6,co=254
               opencirc
               cgoplot,[tbest/1.e6],[lbest],ps=8,symsize=0.6
            endif 
            if keyword_set(oplot_mean) then begin
               filledcirc
               cgoplot,[tmean/1.e6],[lmean],ps=8,symsize=0.6,co=255
               opencirc
               cgoplot,[tmean/1.e6],[lmean],ps=8,symsize=0.6
            endif 
         endelse
         if keyword_set(IMPS) then $  ;Overplot IMPS classification label
            if not keyword_set(mad) then xyouts,parxr[0]+0.1*(parxr[1]-parxr[0]),paryr[0]+0.9*(paryr[1]-paryr[0]),class[jsel],co=254,charsize=1.3,charth=2 else $
               cgtext,0.85*(parxr[1]),10.,class[jsel],co=254,charsize=1.3,charth=2,/align
      endif ;psplots ge 2
  endfor ;isrc

  if keyword_set(IMPS) then begin
     for jcl=0,nclass-1 do close,uclass[jcl]
     print,class+':  ',count
  endif 

  if psplots ge 2 then begin
     cgps_close
     !p.multi=0
  endif 


;Cumulative HRD -- trivial!
  cumHRD = total(hrd,3)
  if not keyword_set(mad) then begin
     xbins = alog10(xt)
     ybins = alog10(yl)
  endif else begin
     xbins = xt/1.e6
     ybins = yl
  endelse 


;Plot ensemble distribution function, if desired.
  if psplots eq 1 or psplots ge 3 then begin
     srclen = strlen(sourcelist)
     sourceliststub = strmid(sourcelist,0,srclen-4)
     if not keyword_set(mad) then plotname1 = 'composite_HRD' else plotname1 = 'composite_MAD'
     if keyword_set(restrict_ages) then plotname1 = plotname1 + '_agecons'
     if keyword_set(age_iso) then begin
        if age_iso[0] ge 10. then fmt = '(F4.1)' else fmt = '(F3.1)'
        plotname1 = plotname1 + string(age_iso[0],format=fmt)+'Myr'
     endif 
;     ps_open, sourceliststub + '_' + plotname1,/encap,/co
     cgps_open, sourceliststub + '_' + plotname1+'.eps',/encap
     device,xsize=7.5,ysize=7,/inches

;Get optimal contour levels
     HRDsrt = cumHRD[sort(cumHRD)]
     HRDsrt = HRDsrt[where(HRDsrt gt 0)]
     pdfHRD = total(HRDsrt,/cumulative)/total(HRDsrt)
     toplev = 1.0   ; top contour
     botlev = 0.05   ; bottom contour
     whrtop = where(pdfHRD ge toplev,ntop)
     if ntop ne 0 then levmax = HRDsrt[whrtop[0]] else levmax = max(HRDsrt)
     whrbot = where(pdfHRD le botlev,nbot)
     if nbot ne 0 then levmin = HRDsrt[whrbot[nbot-1]] else levmin = min(HRDsrt)
; Linear + 0 Contour Levels...     
     nlevs = 58
     dlogc = ((levmax)-(levmin))/nlevs
     clevs = [0, levmin+(dlogc*findgen(nlevs+1))]
     if not keyword_set(mad) then begin ;Make pHRD
        xr = parxr
        yr = paryr  
        cgcontour,cumHRD,xbins,ybins,xtitle='log(T!Deff!N) [K]',ytitle='log(L!Dbol!N) [L'+sunsymbol()+']',xthick=9,ythick=9,charthick=7,levels=clevs,/fill,xr=xr,/xsty,/ysty,yr=yr,position=[0.11,0.11,0.85,0.95];,charsize=2
        numstars = string(nsrc)
        numstars = strtrim(numstars,1)
        if keyword_set(region_name) then begin
           cgtext, 3.6, 3.95,region_name, color=193, charsize=2, charthick=6,/align
           cgtext, 3.6, 3.6,numstars+' stars', color=193, charsize=2, charthick=6,/align
        endif else $
           cgtext, 3.875, 3.95,numstars+' stars', color=193, charsize=2, charthick=6
        oplot_siess,/tracks,color=0, thick=7, /annotate, charsize=1.8, charthick=8
        oplot_siess,/tracks,color=255, thick=5, /annotate, charsize=1.8, charthick=5

        cgaxis,!x.crange[0],!y.crange[0],xr=xr,/xsty,color=255,xtickformat='(A1)',xth=8 
        cgaxis,!x.crange[0],!y.crange[1],xr=xr,/xsty,color=255,xtickformat='(A1)',xth=8,xaxis=1
        cgaxis,!x.crange[0],!y.crange[0],yr=yr,/ysty,color=255,ytickformat='(A1)',yth=8,yaxis=0
        cgaxis,!x.crange[1],!y.crange[0],yr=yr,/ysty,color=255,ytickformat='(A1)',yth=8,yaxis=1
 
     endif else begin  ;Make Mass-Age diagram
  ;NOTE explicitly set plotting ranges!
        xr = parxr
        yr = paryr
        cgcontour,cumHRD,xbins,ybins,xtitle='Age [Myr]',ytitle='Mass [M'+sunsymbol()+']',xthick=8,ythick=8,charthick=7,levels=clevs,/fill,/xsty,/ysty,/xlog,/ylog,xr=xr,yr=yr,position=[0.11,0.11,0.85,0.95];charsize=2
; Shade Main Sequence
        numstars = string(nsrc)
        numstars = strtrim(numstars,1)
        xstarcountpos = 1.0
        cgtext, xstarcountpos, 0.2,numstars+' stars', color=193, charsize=2, charthick=6
        cgtext, 0.15,15, 'ZAMS', color=255, charsize=1.75, charthick=5
;Overplot ZAMS
        oplot_siess_zams,color=0,thick=10
        oplot_bm_zams,color=0,thick=10
        oplot_siess_zams,color=255,thick=8
        oplot_bm_zams,color=255,thick=8
;Overplot disk boundary lines
        cgoplot,replicate(3.,2),[0.1,1],co=0,th=10,li=1
        cgoplot,[3.,0.1],[1.,10.],co=0,th=10,li=1
        cgoplot,replicate(3.,2),[0.1,1],co=255,th=8,li=1
        cgoplot,[3.,0.1],[1.,10.],co=255,th=8,li=1

;        axis,xr[0],yr[0],/xsty,color=255,xtickformat='(A1)',xth=8,/xlog 
;        axis,xr[0],yr[1],/xsty,color=255,xtickformat='(A1)',xth=8,xaxis=1,/xlog
;        axis,xr[0],yr[0],/ysty,color=255,ytickformat='(A1)',yth=8,yaxis=0,/ylog
;        axis,xr[1],yr[0],/ysty,color=255,ytickformat='(A1)',yth=8,yaxis=1,/ylog
     endelse 

;Overplot best-fit models, if desired
     if keyword_set(oplot_best) then begin
        filledcirc
        cgoplot,xbestarr,ybestarr,ps=8,symsize=0.4,co=254
        opencirc
        cgoplot,xbestarr,ybestarr,ps=8,symsize=0.4
     endif 

;Overplot weighted mean parameter values, if desired
     if keyword_set(oplot_mean) then begin
        filledcirc
        cgoplot,xmeanarr,ymeanarr,ps=8,symsize=0.4,co=255
        opencirc
        cgoplot,xmeanarr,ymeanarr,ps=8,symsize=0.4
     endif 

 ;Overplot histogram binsize 
     xnudge = nxbin/20
     ynudge = nybin/5
     ind_xplot = where(xbins ge min(xr) and xbins le max(xr),nxplot)
     ind_yplot = where(ybins ge min(yr) and ybins le max(yr),nyplot)
        
     if nxplot gt xnudge + 1 and nyplot gt ynudge + 1 then begin
        xbinbox = [ xbins[ind_xplot[xnudge]], xbins[ind_xplot[xnudge+1]] ]
        ybinbox = [ ybins[ind_yplot[nyplot-ynudge]], ybins[ind_yplot[nyplot-ynudge+1]] ]
        cgoplot,xbinbox,replicate(mean(ybinbox),2),co=0,th=8
        cgoplot,replicate(mean(xbinbox),2),ybinbox,co=0,th=8
        cgoplot,xbinbox,replicate(mean(ybinbox),2),co=255,th=6
        cgoplot,replicate(mean(xbinbox),2),ybinbox,co=255,th=6
     endif else print,'WARNING! No valid points for plotting binsize!'

     ;Add a colorbar 

      ncol = 2
      colorbar = dblarr(ncol,n_elements(clevs))
      for i=0,ncol-1 do colorbar[i,*] = clevs
      cgcontour,colorbar,[0,1],clevs,/fill,levels=clevs,/noerase,position=[0.88,0.15,0.91,0.9],ysty=1,xsty=4,ytickform='(A1)',ythick=5,yticklen=0.001
      cgoplot,!x.crange,[0,0],th=2
      cgaxis,1,0,yaxis=1,ytitle='N per bin',/ysty,yticklen=0.25,ythick=5,charth=6, charsize=1

     if keyword_set(psplots) then cgps_close
  endif ;End plotting distribution function for the ensemble 

end
