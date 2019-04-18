pro oplot_siess,tracks=tracks,color=color,thick=thick, annotate=annotate, charsize=charsize, charthick=charthick

                                ;Production version, updated to use
                                ;Coyote Graphics commands  15 April 2019
  
  ages = [0.1,1.,3.,10.,30.]

  niso = n_elements(ages)
  for j=0,niso-1 do begin
     age = ages[j]
     isochrone_siess,age,logtiso,logliso                 
     cgoplot,logtiso,logliso,color=color,thick=thick
  endfor 

  if keyword_set(tracks) then begin
     trackmass = [7.,5.5,4,3.,2.,1.]
     nmass = n_elements(trackmass)
     for jt=0,nmass-1 do begin
        mass = trackmass[jt]
        max_age = max(ages)
        pmstrack_siess,mass,max_age,logttrack,logltrack
        cgoplot,logttrack,logltrack,li=2,color=color,thick=thick
        if keyword_set(annotate) then begin
           annot = ['7 M','5.5 M','4 M','3 M','2 M','1 M']+sunsymbol()
           maxlogt = max(logttrack,ind)
           cgtext, 4.55, -0.7,'Isochrones: 30, 10, 3, 1, 0.1 Myr', charsize=charsize, charthick=charthick, color=color
           cgtext, logttrack[ind]+0.04, logltrack[ind], annot[jt], charsize=charsize, charthick=charthick, color=color, /align
;           xyouts, 3.875, 3.95,numstars+' stars', color=193, charsize=2, charthick=6
        endif   
     endfor 
  endif 

;stop
end
