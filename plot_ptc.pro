  ;+
  ; :Author: bob.zavala
  ;-
  
  ; PRO ptc 
  ; 
  ; Procedure to perform Photon Transder Curve analysis following:
  ;   "Photon Transfer DN -> Lambda" by James R. Janesick 

pro ptc

  ; The PTC procedure so the COMMON block can be made. Eventually 
  ; I want to create a PTC 'object' and do all this by inheritance
  ; in an object oriented manner. 
  
  ; The user needs to COMPILE this part, calling ptc or running 
  ; ptc at the command .
  
  COMMON SHARED,imagesDir,rawFilter,oneImstatLine 

  ; The path specificaiton under which the individual image date
  ; directories are stored
  imagesDir = '/mnt/nofs/projects/solarSystemEphem/images/'
  ; The filter to use in the DIALOG_PICKFILE to locate the 
  ; selections for the raw signal data in DN
  rawFilter = 'raw*.txt'

  ; Create the onme line structure to hold the imstat data for 
  ; PTC analysis. As this is needed it may be rpelicated the 
  ; required number of times to accomodate the size required.
  ; ULONG = unsigned LONG integer used as those values are positive.
  
  oneImstatLine = {filename:' ', npix:ULONG(0), mean:0.0, median:0.0, sigma:0.0, $ 
                  min:ULONG(0), max:ULONG(0)}
                  
  print,'**********'
  print,'Photon Transfer Curve precursor items set.'
  print,'**********'

end 
  
pro ptc_temp

  ; Make an ASCIII template to read the statisitcs files created by 
  ; imstatistics
  
  COMMON SHARED,imagesDir,rawFilter
  
  ptc
  
  astatFile  =  DIALOG_PICKFILE(FILTER = rawFilter)
  stat_templ = ASCII_TEMPLATE(astatFile)

  SAVE, stat_templ, FILENAME=imagesDir+'raw.imstat.template.sav',/VERBOSE

END

PRO plot_ptc

  ; Procedure to plot a PTC

  COMMON SHARED,imagesDir,rawFilter,oneImstatLine 
  
  astatFile =  DIALOG_PICKFILE(FILTER = rawFilter)
  RESTORE,FILENAME=imagesDir+'raw.imstat.template.sav'
  statData = READ_ASCII(astatFile, TEMPLATE=stat_templ)
  
  ; Extract the year and day from a filename
  
  yearDay = STRMID(statData.FILENAME[0],0,7)
  
  ; set up the PTC plot
  
  xs = DINDGEN(400000, START=1.0d)
  ys=xs
  
  ; Average the median of the first 100 biases SIGNAL and NOISE
  ; as an x-y pair
  biases = WHERE(statData.obstype EQ 'BIAS')
  avg_bias = [MEAN(statData.median[biases]),MEAN(MEAN(statData.sigma[biases]))]
    
  ; Check if this is the short exposure time day
  
  ; IF (yearDay EQ 'g23d313') || (yearDay EQ 'g23d314')THEN BEGIN
    
    ; The odd numbered indices in uniqDomeFlat need to be averaged to plot the
    ; noise to signal for the PTC for the low light level linearity
    ; Also make an array to hold the short linearity data
    domeFlats = WHERE(statData.OBSTYPE EQ 'Dome_Flat')
    ; Get the indices of the unique Dome_Flat exposure times
    ; These are indices to PUT INTO the domeFlats array
    uniqDomeFlat_expt = UNIQ(statData.EXPTIME[domeFlats])
    ; short_lin_index = [1,3,5,7,9]
    ; Do not hard code the above but select logically for
    ; g23d313
    short_lin_index = WHERE(statdata.exptime[domeFlats[uniqdomeflat_expt]] LT 1.5)
    short_signal = FLTARR(n_elements(short_lin_index))
    short_noise  = FLTARR(n_elements(short_lin_index))
    
    
    ; stop
    ; OK, I do know there are 5 indices, use that knowledge
    ; Find the 0.1 sec exposures: first index inshort_lin_index
    loc_point1 = WHERE(statData.EXPTIME EQ statData.exptime[domeFlats[uniqDomeFlat_expt[1]]])
    ; Find the 0.2 sec exposures: second index in short_lin_index
    loc_point2 = WHERE(statData.EXPTIME EQ statData.exptime[domeFlats[uniqDomeFlat_expt[3]]])
    ; Find the 0.4 sec exposures: thrid index in short_lin_index
    loc_point3 = WHERE(statData.EXPTIME EQ statData.exptime[domeFlats[uniqDomeFlat_expt[5]]])
    ; Find the 0.8 sec exposures
    loc_point4 = WHERE(statData.EXPTIME EQ statData.exptime[domeFlats[uniqDomeFlat_expt[7]]])
    ; Find the 1.0 sec exposures
    loc_point5 = WHERE(statData.EXPTIME EQ statData.exptime[domeFlats[uniqDomeFlat_expt[9]]])
    
    ; Insert the means into the array short_lin_dat
   
    short_signal[0] = [MEAN(statData.median[loc_point1])] & short_noise[0] = [MEAN(statData.sigma[loc_point1])]
    short_signal[1] = [MEAN(statData.median[loc_point2])] & short_noise[1] = [MEAN(statData.sigma[loc_point2])]
    short_signal[2] = [MEAN(statData.median[loc_point3])] & short_noise[2] = [MEAN(statData.sigma[loc_point3])]
    short_signal[3] = [MEAN(statData.median[loc_point4])] & short_noise[3] = [MEAN(statData.sigma[loc_point4])]
    short_signal[4] = [MEAN(statData.median[loc_point5])] & short_noise[4] = [MEAN(statData.sigma[loc_point5])]
    
    stop
   ; ENDIF
   
  ; Find the long linearity signal Dome_Flat data  
  long_lin_index = WHERE((statdata.exptime[domeflats] GT 1.0) AND (statdata.exptime[domeflats] NE 10.0)) 
  IF long_lin_index[0] NE -1 THEN BEGIN
    ; Do these ONLY IF WHERE does not return a -1
    ; Define the arrays to receive the single exposure linearity data
    long_signal = FLTARR(n_elements(long_lin_index))
    long_noise  = FLTARR(n_elements(long_lin_index))
    
    ; Now fille those arrays by selecting from the domeflats indices the 
    ; indices that have expoures GT 1.0 AND NE 10.0 seconds.
    long_signal = statdata.median[domeflats[long_lin_index]]
    long_noise  = statdata.sigma[domeflats[long_lin_index]]
    ;stop
  ENDIF
   
  p1 = PLOT(xs, ys, XLOG=1, YLOG=1, XRANGE=[1.0,1e6], YRANGE=[1,3*MAX(statdata.sigma)], $
    XTITLE='Signal (DN)', YTITLE='NOISE (DN)', FONT_SIZE=12, $
    FONT_STYLE='Bold', FONT_NAME='Times', $
    Title = yearDay+' PTC curve for PI Sophia CCD Camera', $
    XTHICK=2, YTHICK=2, /NODATA)
 
  ; Begin plotting the raw signal and noise 
  p2 = PLOT(short_signal,short_noise,'2D',NAME='Raw Signal (DN)',SYM_FILLED=1,/OVERPLOT)
  
  ; Now get the bias on there
  bias_x = [avg_bias[0],avg_bias[0]]
  bias_y = [avg_bias[1],avg_bias[1]]
  
  ; Plot the bias
  
  p3 = PLOT(bias_x, bias_y,"2o",COLOR='blue',NAME='Offset (DN)',/OVERPLOT)

  ; Check the 10.0 second lamp levels
  ; Find the lamp level check frames
  lamp_lev_index = WHERE(statData.exptime EQ 10.0)
  lamp_signal = statData.median[lamp_lev_index]
  lamp_noise  = statData.sigma[lamp_lev_index]
  lamp_count = INDGEN(n_elements(lamp_lev_index),START=1,/ULONG) 
  
  ; If long exposres that are not average 1.5 seconds or more are present 
  ; plot these
    
  IF long_lin_index[0] NE -1 THEN BEGIN
    
    ; Overplot the single exposure raw signal and noise PTC values
    p4 = PLOT(long_signal, long_noise,'2D',NAME='Raw Signal (DN)',SYM_FILLED=1,/OVERPLOT)     
    
  ENDIF
 
  ; Plot the lamp levels

  plamp = ERRORPLOT(lamp_count,lamp_signal,lamp_noise,'-2H',SYM_COLOR='red', $ 
                     SYM_FILLED=1,XRANGE=[0,n_elements(lamp_count)+1],$
                     YRANGE=[0.9*MIN(lamp_signal),1.1*MAX(lamp_signal)], $ 
                     FONT_STYLE='Bold', FONT_NAME='Times', $
                     XTITLE='Lamp Check #',YTITLE='Lamp Signal (DN)', $
                     TITLE=yearDay+' Linearity Lamp Level Tests')
                     
  ; Use the TEXT function to add some UT times
  the_times = STRARR(5)
  time_indices = [0,1,2,3,10]
  the_times = [STRMID(statdata.utcobs[lamp_lev_index[time_indices]],0,8)]
  
  ;t1 = TEXT(0.3,4710,the_times[0],/DATA,FONT_SIZE=10,FONT_STYLE='Italic')
  ;t2 = TEXT(1.3,4310,the_times[1],/DATA,FONT_SIZE=10,FONT_STYLE='Italic')
  ;t3 = TEXT(2.3,4600,the_times[2],/DATA,FONT_SIZE=10,FONT_STYLE='Italic')
  ;t4 = TEXT(3.3,4280,the_times[3],/DATA,FONT_SIZE=10,FONT_STYLE='Italic')
  ;t11 = TEXT(10.3,4270,the_times[4],/DATA,FONT_SIZE=10,FONT_STYLE='Italic')

  ; save the progress!
  
  SAVE,/VARIABLES,/VERBOSE,FILENAME=yearDay+'lin_dat.sav'

  stop

END