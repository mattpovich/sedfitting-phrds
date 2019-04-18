;Compute an PMS track for 0.1-7 Msun over specified age range using the (solar-metallicity)
;models of Siess et al. 2000

pro pmstrack_siess, mass, max_age, logteff, logl

max_age = 1.e6*max_age   ;Convert Myr to yr

if n_params() lt 3 then begin
   print,'syntax - pmstrack, mass, max_age, logteff, logl'
   return
endif

  filelist = file_search(STRSPLIT(!PATH, PATH_SEP(/SEARCH_PATH),/EXTRACT),'OV02/hrdlist.txt')
  filelist = filelist[0]
  
  filedir = strmid(filelist,0,strlen(filelist)-strlen('hrdlist.txt'))

  nmass = file_lines(filelist)
  nmass = nmass[0]
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

  file_mass = float(strmid(filenames,1,3))
  jnk = min(abs(file_mass - mass),massi)
  
  loadfilename = filedir + '/'+ filenames[massi]
  
  ;Read  file for specified MASS
  nlines = file_lines(loadfilename)
  lines = strarr(nlines)
  line=''
  openr,u,loadfilename,/get_lun
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

  whrtrack = where(age le max_age,nage)
     if nage eq 0 then message,'Problem with AGE selection!'
     logteff = alog10(teff[whrtrack])
     logl = alog10(l[whrtrack])

end
