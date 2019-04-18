pro pms_pars2idl,infile,target_dir,data_parent=data_parent,filter=filter

;Read the (expicitly-formatted) ASCII output (INFILE) of the python sedfitter,
;run using the custom PMS stellar models, into IDL structures. Creates one structure for each source and saves in an IDL save
;file in the specified directory (TARGET_DIR).
;Also creates a list of the output filenames for guiding further automated
;analysis routines.

;INPUTS
;         INFILE     'string' -- Path to ASCII parameter file containing SED
;                     fitting results. EXAMPLE: 'pars_xpms_g4.txt'  
;         TARGET_DIR 'string': Path to the output directory for the IDL save
;                       files. EXAMPLE: 'results_xpms'
;
;KEYWORD PARAMETERS
;        DATA_PARENT    'string' -- Path to a common fitter data format
;                                   file used for all fits. If this
;                                   keyword is set, the output save
;                                   files will contain a source S
;                                   structure in addition to the
;                                   parameter P structure.
;         FILTER=         'string' -- Filter the model fits saved
;                                     using the chi2 criteria defined
;                                     in the SED fitter
;                                     manual. Examples: 'F1','N10','E2'
;                                  
;
;
;WARNING: Unless you are working in a directory (normally the
;sedfitter directory) with the models_pms symlink active, you will
;need to specify the path to the models_pms directory using the MODELS keyword.
;
;MSP adapted FITS2IDL_PMS using R17_PARS2IDL 29 August 2018


if n_params() lt 2 then begin
    print,'syntax: pms_fits2idl, infile, target_dir [, data_parent=, filter=]'
    return
endif


;create output directory
  if file_test(target_dir) then exists = 1
  if keyword_set(exists) then begin
     check = ''
     read,check,prompt=target_dir+': Directory exists. Overwrite? y/[n]  '
     if check ne 'y' then begin
        print,'Aborting to avoid overwriting files.'
        return
     endif else file_delete,target_dir,/recursive
  endif 
  spawn,'mkdir '+target_dir
  target = target_dir + '/'

;initialize ASCII I/O
 
  if keyword_set(data_parent) then begin
     n_datalines = file_lines(data_parent)
     datalines = strarr(n_datalines)
     openr,w,data_parent,/get_lun
     line = '' 
     for j=0L,n_datalines-1 do begin
        readf,w,line
        datalines[j] = line
        cells = line.Split(' ')
        cells = cells[where(cells ne '',n_cols)]
        if j eq 0 then begin  ;What's wrong with initializing arrays this way?
           n_bands = n_cols/3 - 1
           datanames = strarr(n_datalines)
           l = fltarr(n_datalines)
           b = fltarr(n_datalines)
           valid = intarr(n_datalines,n_bands)
           flux = fltarr(n_datalines,n_bands)
           flux_error = fltarr(n_datalines,n_bands)
        endif
        datanames[j] = cells[0]
        l[j] = float(cells[1])
        b[j] = float(cells[2])
        valid[j,*] = fix(cells[3:3+n_bands-1])
        flux[j,*] = float(cells[3+n_bands:*:2])
        flux_error[j,*] = float(cells[4+n_bands:*:2])
     endfor 
     close,w
     free_lun,w

  endif 

  
  n_source = 0L                 ;RUNNING counter of sources processed
  openw,v,target_dir+'/sourcelist.txt',/get_lun ;print sourcelist
  
  n_srcinset = 0L
  n_head = 3                ;number of header lines in parameter file for each source
  h0 = file_lines(infile) - n_head
  h = h0
  head = strarr(n_head)
  headline = ''
  infoline = ''
  modelline = ''
  openr,u,infile,/get_lun

  print,'Processing parameter file: '+infile
     
  for i=0, n_head-1 do begin
     readf,u,headline
     head[i] = headline
  endfor 
  colnames = head[1].Split(' ')
  colnames = colnames[where(colnames ne '',n_colnames)]
     
  while keyword_set(h) do begin

     readf,u,infoline
     fit_info = infoline.Split(' ')
     fit_info = fit_info[where(fit_info ne '',n_info)]
     if n_info ne 3 then message,'BUG Alert! Should be 3 items in the last header line for each source in parfile.'
      
     source_name = fit_info[0]
     ndata = fix(fit_info[1])
     n_fits = long(fit_info[2])
        
                                ;Construct source structure (OPTIONAL)
     if keyword_set(data_parent) then begin
        ind_data = where(datanames eq source_name,n_datamatch)
        if n_datamatch eq 1 then begin
           s = {DESIG:source_name,L:l[ind_data],B:b[ind_data],VALID:valid[ind_data,*],$
                F:flux[ind_data,*],DF:flux_error[ind_data,*]}
        endif else begin
           print,'WARNING! No unique match to source '+source_name+' found in '+data_parent +'!'
           print,'No S structure was created.'
           stop
        endelse 
     endif

                                ;Create structure table for PMS-specific model
                                ;parameters
     f_nan    = !VALUES.F_NAN
     template_row = {$
                    MODEL_NAME       :''   ,$
                    CHI2             :f_nan,$
                    AV               :f_nan,$
                    SCALE            :f_nan,$
                    AGE              :f_nan,$
                    MASS             :f_nan,$
                    LBOL             :f_nan,$
                    REFF             :f_nan,$
                    TEFF             :f_nan,$
                    LOGG             :f_nan $
                    }
     
     pars = replicate(template_row, n_fits)
     
     for j=0L,n_fits-1 do begin
        readf,u,modelline
        modelcells = modelline.Split(' ')
        modelcells = modelcells[where(modelcells ne '',n_modcells)]
        pars[j].model_name = modelcells[1]
        pars[j].chi2 = modelcells[2]
        pars[j].av = modelcells[3]
        pars[j].scale = modelcells[4]
        pars[j].age = modelcells[5]
        pars[j].mass = modelcells[6]
        pars[j].lbol = modelcells[7]
        pars[j].reff = modelcells[8]
        pars[j].teff = modelcells[9]
        pars[j].logg = modelcells[10]
     endfor

     p = temporary(pars)
     
;FILTER fits by chisq, if desired. NOTE: This works ACROSS ALL MODEL
;SETS as recommended by R17.
     if keyword_set(filter) then begin
        ftype = strmid(filter,0,1)
        fnum = float(strmid(filter,1))
        chi2 = p.chi2
        chi2_best = min(chi2)
        case 1 of 
           ftype eq 'N': begin
              fnum = fix(fnum)
              srtchi2 = sort(chi2)
              ind_good = srtchi2[0:fnum-1]
              n_fits2 = fnum
           end 
           ftype eq 'C': begin
              ind_good = where(chi2 lt fnum,n_fits2)
           end 
           ftype eq 'D': begin
              ind_good = where(chi2 - chi2_best lt fnum,n_fits2)
           end 
           ftype eq 'E': begin
              ind_good = where(chi2/ndata lt fnum,n_fits2)
           end 
           ftype eq 'F': begin
              ind_good = where((chi2 - chi2_best)/ndata lt fnum,n_fits2)
           end 
           else: print,'Invalid FILTER type specified (N C D E or F required)! SKIPPING filtering on Chisq.' 
        endcase
        p = p[ind_good]
     endif 
        
     if isa(s) then save,s,p,file=target_dir + '/' + source_name + '.dat' $
     else save,p,file=target_dir + '/' + source_name + '.dat'
     n_source++
     h = h - n_fits - 1

     printf,v,source_name + '.dat'
     if n_source ne 0 and n_source mod 100 eq 0 then print,n_source,' sources completed.'
  endwhile                      

  close,v
  free_lun,v
  print,'Successfully generated parameter files for',n_source,' sources.'

end
