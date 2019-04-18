;Compute an isochrone for 0.1-7 Msun at a given age using the (solar-metallicity)
;models of Siess et al. 2000

pro isochrone_siess, iso_age, logteff, logl

iso_age = 1.e6*iso_age   ;Convert Myr to yr

if n_params() lt 3 then begin
   print,'syntax - isochrone, age, logteff, logl'
   return
endif

  filelist = file_search(STRSPLIT(!PATH, PATH_SEP(/SEARCH_PATH),/EXTRACT),'OV02/hrdlist.txt')

  filelist = filelist[0]
  filedir = strmid(filelist,0,strlen(filelist)-strlen('hrdlist.txt'))

  nmass = file_lines(filelist)
  line = ''
  filenames = strarr(nmass)
  openr,u,filelist,/get_lun
  for j=0,nmass-1 do begin
     readf,u,line
     filenames[j] = line
  endfor 
  close,u
  free_lun,u
;  readcol,format='A',filelist,filenames
  filenames = filedir + '/' + filenames

;Arrays to hold isochrone data (may need to interpolate these later)
  logteff = fltarr(nmass)
  logl = fltarr(nmass)
  
  ;Read each file in turn
  for i=0,nmass-1 do begin
;     readcol,filenames[i],model,l,mbol,reff,rstar,teff,rho,logg,m,age
     nlines = file_lines(filenames[i])
     lines = strarr(nlines)
     line=''
     openr,u,filenames[i],/get_lun
     for j=0L,nlines-1 do begin
        readf,u,line
        lines[j] = line
     endfor 
     close,u
     free_lun,u
     lines = lines[3:*]
     teff = float(strmid(lines,46,6))
     l = float(strmid(lines,10,10))
     age = float(strmid(lines,86,16))
     jnk = min(abs(age-iso_age),agei)
;     if agei ne 1 then message,'Problem with AGE selection!'
     logteff[i] = alog10(teff[agei])
     logl[i] = alog10(l[agei])
  endfor 

end
