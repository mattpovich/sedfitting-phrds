;Wrapper program for MAD_HIST2_PMS. Make 2 histograms from SEd fitting
;results, one each for
;stellar age and mass. Uses two iterations, first to identify the best
;isochrone, second to produce the mass function based on that isochrone.

pro massage_plots_pms,target_dir, sourcelist=sourcelist, mass_offset1=mass_offset1, mass_offset2=mass_offset2, age_offset=age_offset, region_name=region_name

  if n_params() lt 1 then begin
     print,'syntax -- massage_plots_pms, target_dir,sourcelist=[sourcelist.txt],mass_offset=[0],age_offset=[0],region_name=[]'
     return
  endif
  
  if not keyword_set(sourcelist) then sourcelist='sourcelist.txt'
  if not keyword_set(mass_offset) then mass_offset = 0
  if not keyword_set(age_offset) then age_offset = 0
  
     ;Get ages and age distributions
     mad_hist2_pms, target_dir, agebins, massbins, age_hist, mass_hist, age, age_err, sourcelist=sourcelist,restrict=3, region_name=region_name, offset=mass_offset1,age_offset=age_offset
     age_iso = [age,age_err]
     
     ;Get mass functions using isochronal ages and PLOT
     mad_hist2_pms, target_dir, agebins, massbins, age_hist, mass_hist, age, age_err, sourcelist=sourcelist,restrict=3, region_name=region_name, offset=mass_offset2,age_offset=age_offset,/input_age,age_iso=age_iso
  
end 
