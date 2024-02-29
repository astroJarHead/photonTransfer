  ;+
  ; :Author: bob.zavala
  ;-
  
  ; PRO ptc 
  ; 
  ; Procedure to perform Photon Transder Curve analysis following:
  ;   "Photon Transfer DN -> Lambda" by James R. Janesick 

PRO ptc

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

END 

;******************************
  
PRO ptc_temp

  ; Make an ASCIII template to read the statisitcs files created by 
  ; imstatistics
  
  ptc
  
  COMMON SHARED,imagesDir,rawFilter
  
  ptc
  
  astatFile  =  DIALOG_PICKFILE(FILTER = rawFilter)
  stat_templ = ASCII_TEMPLATE(astatFile)

  SAVE, stat_templ, FILENAME=imagesDir+'raw.imstat.template.sav',/VERBOSE

END

;******************************

PRO plot_ptc

  ; Procedure to plot a PTC

  ptc
  
  COMMON SHARED,imagesDir,rawFilter,oneImstatLine 
  
  astatFile =  DIALOG_PICKFILE(FILTER = rawFilter)
  RESTORE,FILENAME=imagesDir+'raw.imstat.template.sav'
  statData = READ_ASCII(astatFile, TEMPLATE=stat_templ)
  
  ; Extract the year and day from a filename
  
  yearDay = STRMID(statData.FILENAME[0],0,7)
  
  ; set up the PTC plot dummy indices
  
  xs = DINDGEN(400000, START=1.0d)
  ys=xs
  
  ; If darks were taken a linearity test might be possible. Let's see.
  
  domeFlats = WHERE(statData.OBSTYPE EQ 'Dome_Flat', domeFlatCount)
   
  IF domeFlatCount GT 1 THEN BEGIN 
    
    ; If this is a short exposure time day the short exposures need 
    ; to be averaged, e.g. many 0.1 second domeFlats averaged. 
      
    ; The odd numbered indices in uniqDomeFlat need to be averaged to plot the
    ; noise to signal for the PTC for the low light level linearity
    ; Also make an array to hold the short linearity data
  
    ; Get the indices of the unique Dome_Flat exposure times
    ; These are indices to PUT INTO the domeFlats array
    uniqDomeFlat_expt = UNIQ(statData.EXPTIME[domeFlats])
    ; short_lin_index = [1,3,5,7,9]
    ; Do not hard code the above but select logically for
    ; dome flats with multiple exposures. I know these from the 
    ; datFile are exposure times < 1.5 seconds. 
    short_lin_index = WHERE(statdata.exptime[domeFlats[uniqdomeflat_expt]] LT 1.5)
    short_signal = FLTARR(n_elements(short_lin_index))
    short_noise  = FLTARR(n_elements(short_lin_index))
    
    ; I know there are five sets of multiple exposure times to average so I 
    ; will use this a priori knowledge.
    loc_point1 = WHERE(statData.EXPTIME EQ statData.exptime[domeFlats[uniqDomeFlat_expt[short_lin_index[0]]]])
    ; Find the 0.2 sec exposures: second index in short_lin_index
    loc_point2 = WHERE(statData.EXPTIME EQ statData.exptime[domeFlats[uniqDomeFlat_expt[short_lin_index[1]]]])
    ; Find the 0.4 sec exposures: thrid index in short_lin_index
    loc_point3 = WHERE(statData.EXPTIME EQ statData.exptime[domeFlats[uniqDomeFlat_expt[short_lin_index[2]]]])
    ; Find the 0.8 sec exposures
    loc_point4 = WHERE(statData.EXPTIME EQ statData.exptime[domeFlats[uniqDomeFlat_expt[short_lin_index[3]]]])
    ; Find the 1.0 sec exposures
    loc_point5 = WHERE(statData.EXPTIME EQ statData.exptime[domeFlats[uniqDomeFlat_expt[short_lin_index[4]]]])
    
    ; Insert the means into the array short_lin_dat
   
    short_signal[0] = [MEAN(statData.median[loc_point1])] & short_noise[0] = [MEAN(statData.sigma[loc_point1])]
    short_signal[1] = [MEAN(statData.median[loc_point2])] & short_noise[1] = [MEAN(statData.sigma[loc_point2])]
    short_signal[2] = [MEAN(statData.median[loc_point3])] & short_noise[2] = [MEAN(statData.sigma[loc_point3])]
    short_signal[3] = [MEAN(statData.median[loc_point4])] & short_noise[3] = [MEAN(statData.sigma[loc_point4])]
    short_signal[4] = [MEAN(statData.median[loc_point5])] & short_noise[4] = [MEAN(statData.sigma[loc_point5])]
    
    ; stop
     
    ; Find the long linearity signal Dome_Flat data  
    long_lin_index = WHERE((statdata.exptime[domeflats] GT 1.0) AND (statdata.exptime[domeflats] NE 10.0)) 
    
    IF long_lin_index[0] NE -1 THEN BEGIN
    
      ; Do these ONLY IF WHERE does not return a -1
      ; Define the arrays to receive the single exposure linearity data
      long_signal = FLTARR(n_elements(long_lin_index))
      long_noise  = FLTARR(n_elements(long_lin_index))
      
      ; Now fill those arrays by selecting from the domeflats indices the 
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
    p2 = PLOT(short_signal,short_noise,'2D',NAME='Averaged exp. (DN)',SYM_FILLED=1,/OVERPLOT)
          
    ; Check the 10.0 second lamp levels
    ; Find the lamp level check frames
    lamp_lev_index = WHERE(statData.exptime EQ 10.0)
    lamp_signal = statData.median[lamp_lev_index]
    lamp_noise  = statData.sigma[lamp_lev_index]
    lamp_count = INDGEN(n_elements(lamp_lev_index),START=1,/ULONG) 
    
    ; If long exposres that are not averaged (1.5 seconds) or more are present 
    ; plot these
      
    ; stop  
      
    IF long_lin_index[0] NE -1 THEN BEGIN
      
      ; Overplot the single exposure raw signal and noise PTC values
      p3 = PLOT(long_signal, long_noise,'2D',NAME='Single exp. (DN)',SYM_FILLED=0,/OVERPLOT)     
      
      ; Now add the legend through p4
      
      leg_p3 = LEGEND(TARGET=[p1,p2,p3], POSITION=[800,1000],/DATA,/AUTO_TEXT_COLOR)
    
    ENDIF ELSE BEGIN
      
      ; No single exposure data, so plot a legend without p3
      
      leg_p2 = LEGEND(TARGET=[p1,p2], POSITION=[800,1000],/DATA,/AUTO_TEXT_COLOR)
      
    ENDELSE
   
    ; What is the lamp power variation in percent?
    
    lamp_sig_delta_perc = 100.0*(mean(lamp_signal) - MIN(lamp_signal))/mean(lamp_signal)
    
    ; Plot the lamp levels
  
    plamp = ERRORPLOT(lamp_count,lamp_signal,lamp_noise,'-2H',SYM_COLOR='red', $ 
                       SYM_FILLED=1,XRANGE=[0,n_elements(lamp_count)+1],$
                       YRANGE=[0.9*MIN(lamp_signal),1.1*MAX(lamp_signal)], $ 
                       FONT_STYLE='Bold', FONT_NAME='Times', $
                       XTITLE='Lamp Check #',YTITLE='Lamp Signal (DN)', $
                       TITLE=yearDay+' Linearity Lamp Level Tests')
    perc_chg_text = 'Mean lamp signal change (%) =: '+STRTRIM(STRING(lamp_sig_delta_perc),1)
    
    t_percent = TEXT(5,5000,perc_chg_text,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')
    
    ; save the progress!
    
    SAVE,/VARIABLES,/VERBOSE,FILENAME=yearDay+'.lin_dat.sav'
  
  ENDIF
  
  ; Now get the darks, make a smaller structure and pass that and the bias numbers 
  ; for the dark count analysis
  
  ; find the darks
  
  the_darks = WHERE(statdata.obstype EQ 'DARK')
  
  ; If we have darks, extract these from the statdata structure and pass the 
  ; dark frames to the proceduer do_darks for analysis with the bias measurement. 
  
  IF the_darks[0] NE -1 THEN BEGIN
  
  
    darkData = {filename:statdata.filename[the_darks],statsec:statdata.statsec[the_darks], $ 
               npix:statdata.npix[the_darks],mean:statdata.mean[the_darks],median:statdata.median[the_darks], $ 
               sigma:statdata.sigma[the_darks],exptime:statdata.exptime[the_darks],$ 
               utcobs:statdata.utcobs[the_darks],obstype:statdata.obstype[the_darks]}
             
    do_darks,darkData
    
  ENDIF

  stop

END

;******************************

PRO do_darks,darkData

  ; Procedure to perform analysis of dark frames using a 
  ; Dark Transfer Curve (DTC) per Chapter 11 of Janesick 
  ; ¨Photon Transfer DN -> Lambda"
  
  ; RELATED:
  ;
  ; This procedure called by PRO plot_ptc IF darks are present on the night
  ; processed by plot_ptc
  
  ; Plot the DTC
  
  ; Extract the year and day from a filename

  COMMON SHARED,imagesDir,rawFilter,oneImstatLine

  yearDay = STRMID(darkData.FILENAME[0],0,7)

  
  ; set up the DTC plot dummy indices

  xs = DINDGEN(400000, START=1.0d)
  ys=xs
  
  ; Setup a theoretical Poisson noise
  
  small_xs = WHERE(xs LT 10000)
  ys_pois = SQRT(small_xs)
  
  ; Initialize the DTC plot
  
  ; How large should the Signal (X axis) range go?
  ; Use the lower X limit as the bias_x/10.0
  maxSig = MAX(darkData.median)
  IF maxSig LT 10000.0 THEN BEGIN
    maxXrange = 1e4
  ENDIF ELSE BEGIN
    maxXrange = 1e5
  ENDELSE
  
  ; Set the Y axis range
  maxSigma = MAX(darkData.sigma)
  maxYrange = 1.33*maxSigma
  
  p1 = PLOT(xs, ys_pois,'--', XLOG=1, YLOG=1, XRANGE=[1,maxXrange], $ 
    YRANGE=[1,maxYrange], $
    XTITLE='Signal (DN)', YTITLE='NOISE (DN)', FONT_SIZE=12, $
    FONT_STYLE='Bold', FONT_NAME='Times', $
    Title = yearDay+' DTC curve for PI Sophia CCD Camera', $
    XTHICK=2, YTHICK=2,NAME='$\sigma_{D_SHOT}$')
    
  p2 = PLOT(darkData.median,darkData.sigma,'2D',NAME='Dark (DN)',$ 
            SYM_FILLED=0,SYM_THICK=2,/OVERPLOT)
  
  ; Make a linear plot
  
;  p4 = PLOT(darkData.median,dk_nobias,'2D',NAME='Dark - Bias (DN)',SYM_FILLED=1, $
;    XTITLE='Signal (DN)', YTITLE='NOISE (DN)', FONT_SIZE=12, $
;    FONT_STYLE='Bold', FONT_NAME='Times', $
;    Title = yearDay+' Linear DTC curve for PI Sophia CCD Camera', $
;    XTHICK=2, YTHICK=2)

  ; Add a legend
    
  leg_p1 = LEGEND(TARGET=[p1,p2], POSITION=[50,100],/DATA)
    
  ; Assuming Equation 11.18 sigma_D_FPN = D*D_n holds use the 
  ; entire set of dark data for a nighto estiamte the percentage factor 
  ; D_n the Dark Current Quality Factor
  
  ; Use the MOMENT function to get a MEAN and VARIANCE of the dark 
  ; quality factor for this date.
  
  dQualFacts = MOMENT(darkData.sigma/darkData.median)
  
  print,' '
  print,'For '+yearDay+' the Dark Quality Factor $D_{N}$ is:'
  print, STRTRIM(STRING(dQualFacts[0]),1)+' +/- '+STRTRIM(STRING(SQRT(dQualFacts[1])),1)+' %.'
  print,' '
  
 ; Save the dark data and combine it with other nights for the dark DTC
  
  SAVE,/VARIABLES,/VERBOSE,FILENAME=yearDay+'.darkData.sav'
  
END

;******************************

PRO do_dtc

  ptc
  
  COMMON SHARED,imagesDir

  ; In here I will plot the readnoise signal and sigma in DN determined on 
  ; g23d319 from 200 bias frames using:
  ; biases=WHERE(statdata.obstype EQ 'BIAS')
  ; rd_noise_signal = MEAN(statdata.MEDIAN[biases])
  ; rd_noise_dn = MEAN(statdata.SIGMA[biases])
  
  ; Arrays are used for plotting
  rd_noise_signal = [0.3207,0.3207]
  rd_noise_sigma  = [3.9332,3.9332] 

  ; Get a list of all the g??d???.darkData.sav files in this 
  ; directory and count the number of saved files present. 
  
  darkData_save_files = FILE_SEARCH('g??d???.darkData.sav',COUNT=fileCount)
  
  ; How many entries are in the individual nights that have dark data?
  ; Add these up and make one large structure to hold everything.
  
  ; Initialize a variable to hold the count of the entries, and it is positive
  ; and make it long so it will definitely be large enough. 
  ; Also make an integer array to hold the number of entries per night
  
  entryTotal = 0UL
  theEntries = INTARR(fileCount) 
  
  ; Loop over the remaining darkData_save_files
  
  FOR i=0,fileCount-1 DO BEGIN
  
    ; Restore the file
    RESTORE,filename=darkData_save_files[i]
    theEntries[i] = n_elements(darkData.filename)
    entryTotal = entryTotal + theEntries[i]   
  
  ENDFOR

  ; Now make a single structure identical in definition to darkData to 
  ; hold all of the darkData from the individual nights
   
  allDarkData = {filename:STRARR(entryTotal),statsec:STRARR(entryTotal), $
    npix:ULONARR(entryTotal),median:FLTARR(entryTotal), $
    sigma:FLTARR(entryTotal),exptime:FLTARR(entryTotal),$
    utcobs:STRARR(entryTotal),obstype:STRARR(entryTotal)}
  
  ; Now fill in the large allDarkData structure with the data from the individual nights
  ; The indices first range is [0:theEntries[0]-1] and I will use this to 
  ; fill in the first set of data to simplify the counting
  
  RESTORE,filename=darkData_save_files[0]
  
  allDarkData.filename[0:theEntries[0]-1] = darkData.filename[0:theEntries[0]-1]
  allDarkData.statsec[0:theEntries[0]-1] = darkData.statsec[0:theEntries[0]-1]  
  allDarkData.npix[0:theEntries[0]-1] = darkData.npix[0:theEntries[0]-1]
  allDarkData.median[0:theEntries[0]-1] = darkData.median[0:theEntries[0]-1]
  allDarkData.sigma[0:theEntries[0]-1] = darkData.sigma[0:theEntries[0]-1]
  allDarkData.exptime[0:theEntries[0]-1] = darkData.exptime[0:theEntries[0]-1]
  allDarkData.utcobs[0:theEntries[0]-1] = darkData.utcobs[0:theEntries[0]-1]
  allDarkData.obstype[0:theEntries[0]-1] = darkData.obstype[0:theEntries[0]-1]
  
  ; Set the number to increment for the counting loop to increment correctly
  ; This recorsd that the first "theSkip" indices are already filled. 
  
  theSkip = theEntries[0]-1
  
  ; Now enter the data from the remaining dark data night save files
  
  FOR j=1,fileCount-1 DO BEGIN

    ; Restore a file
    RESTORE,filename=darkData_save_files[j]
    
    ; Loop over the files and fill the large structure allDarkData by
    ; using indices NOT loops.
    
    allDarkData.filename[theSkip + 1:theEntries[j] + theSkip] = darkData.filename[0:theEntries[j]-1]
    allDarkData.statsec[theSkip + 1:theEntries[j] + theSkip] = darkData.statsec[0:theEntries[j]-1]
    allDarkData.npix[theSkip + 1:theEntries[j] + theSkip] = darkData.npix[0:theEntries[j]-1]
    allDarkData.median[theSkip + 1:theEntries[j] + theSkip] = darkData.median[0:theEntries[j]-1]
    allDarkData.sigma[theSkip + 1:theEntries[j] + theSkip] = darkData.sigma[0:theEntries[j]-1]
    allDarkData.exptime[theSkip + 1:theEntries[j] + theSkip] = darkData.exptime[0:theEntries[j]-1]
    allDarkData.utcobs[theSkip + 1:theEntries[j] + theSkip] = darkData.utcobs[0:theEntries[j]-1]
    allDarkData.obstype[theSkip + 1:theEntries[j] + theSkip] = darkData.obstype[0:theEntries[j]-1]
    
    ; Now iterate the skip increment
    
    theSkip = theEntries[j] + theSkip

  ENDFOR
    
  ; Plot the combined DTC data from the several nights with a Shot Noise comparison
  
  ; set up the DTC plot dummy indices

  xs = DINDGEN(400000, START=1.0d)
  ys=xs

  ; Setup a theoretical Poisson noise

  small_xs = WHERE(xs LT 10000)
  ys_pois = SQRT(small_xs)

  ; How large should the Signal (X axis) range go?
  ; Use the lower X limit as the bias_x/10.0
  maxSig = MAX(darkData.median)
  IF maxSig LT 10000.0 THEN BEGIN
    maxXrange = 1e4
  ENDIF ELSE BEGIN
    maxXrange = 1e5
  ENDELSE

  ; Set a reasonable maximum for the Y axis range
  maxYrange = 1.33*MAX(allDarkData.sigma)
  
  p1 = PLOT(xs, ys_pois,'--', XLOG=1, YLOG=1, XRANGE=[0.1,maxXrange], $
    YRANGE=[1,maxYrange], $
    XTITLE='Signal (DN)', YTITLE='NOISE (DN)', FONT_SIZE=12, $
    FONT_STYLE='Bold', FONT_NAME='Times', $
    Title = 'Multi-night DTC curve for PI Sophia CCD Camera', $
    XTHICK=2, YTHICK=2,NAME='$\sigma_{D\_SHOT}$')

  p2 = PLOT(allDarkData.median,allDarkData.sigma,'2D', $ 
           NAME='$\sigma_{D\_FPN}$ (DN)',SYM_FILLED=0, $ 
           SYM_THICK=2,/OVERPLOT)
           
  ; Plot the read noise from g23d319
  
  p_rd = PLOT(rd_noise_signal,rd_noise_sigma,'sb2',NAME='$\sigma_{READ} = 3.93 DN$',$ 
             SYM_FILLED=1,/OVERPLOT)

  ; Let DCQF = Dark Current Quality Factor or D_N in Janesick's book 
  ; and using Equation 11.22 fill an array with the values for DCQF
  ; Only use data for which there is a decently high SNR
  
  goodSignal = WHERE(allDarkData.MEDIAN GE 1500.0)
  
  DCQF_arr = allDarkData.SIGMA[goodSignal]/allDarkData.MEDIAN[goodSignal]
  
  ; Estiamte the DCQF and the statistics using the MOMENT function
  
  DCQF = MOMENT(DCQF_arr)
  
  deScaled_pois = SQRT(allDarkData.sigma/DCQF[0])
  
  pChk = PLOT(allDarkData.median,deScaled_pois,'ro',NAME='$\sigma_{D\_CORR}$ (DN)', $ 
             SYM_THICK=2,/OVERPLOT)
  
  ; Add a legend

  leg_p1 = LEGEND(TARGET=[p1,p2,pChk,p_rd], POSITION=[40000,10],/DATA)
  
  ; Label the Dark Current Quality Factor

  D_N_text = 'D$_{N}$ = '+STRTRIM(STRING(DCQF[0],FORMAT='(F5.3)',/PRINT),1)

  D_N = TEXT(0.5,90,D_N_text,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')
  
  ; As the "de-scaling" looks like it works well for SIgnal >= 2900 DN 
  ; use this threshold for detemring K_adc and get the mean. This is 
  ; well behaved so uncertainty is negligible
  
  Sig_29hun = WHERE(allDarkData.median GE 2900.0)
  
  ; Determine K_adc (electrons/DN) using Equation 11.20
  
  K_29hun = allDarkData.median[goodSignal[Sig_29hun]]/(deScaled_pois[goodSignal[Sig_29hun]])^2
  
  K = MEAN(K_29hun)
  
  ; Label K_adc
  
  ; Label the Dark Current Quality Factor

  K_adc_text = 'K$_{ADC} (e^{-}/DN)$ = '+STRTRIM(STRING(K,FORMAT='(F5.3)',/PRINT),1)

  K_adc = TEXT(0.5,60,K_adc_text,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')
  
  ; Now that we know the ADC sensitivity (gain in amateur speak) plot up the 
  ; dark signal versus time
  
  D_etron = K*allDarkData.median
  D_etron_sig = SQRT(D_etron)
  
  ; Linear fit for the dark rate against time
  
  darkRate = LINFIT(allDarkData.exptime,D_etron,MEASURE_ERRORS=D_etron_sig,$ 
                   CHISQR=drate_chisq,COVAR=drate_covar,SIGMA=drate_sigma,$
                   YFIT=drat_yfit)
  
  ; Plot the counts in electrons versus time.
  
  maxTime = MAX(allDarkData.exptime)
  
  maxEtron = 1.5*MAX(D_etron)
  
  plDk = PLOT(allDarkData.exptime,D_etron,'o',XLOG=1, YLOG=1, $
             YRANGE=[1,maxEtron], XRANGE=[1.0,1e5],$ 
             XTITLE='Dark exposure time (sec)', $
             YTITLE='Dark signal (e$^{-}$)', FONT_SIZE=12, $
             FONT_STYLE='Bold', FONT_NAME='Times', $
             Title = 'Multi-night Dark rate curve for PI Sophia CCD Camera', $
             XTHICK=2, YTHICK=2,SYM_THICK=2,NAME='D$_{R}$')

  ; Label the Dark Rate

  D_rate_text = 'D$_{R} (e^{-}/sec)$ = '+STRTRIM(STRING(darkRate[1],FORMAT='(F5.3)',/PRINT),1)
  D_rate_text2= '                 $\pm$ '+STRTRIM(STRING(drate_sigma[1],FORMAT='(F5.3)',/PRINT),1) 
  D_rate = TEXT(2.5,1000,D_rate_text,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')
  D_rate2= TEXT(2.5, 680,D_rate_text2,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')

  ; Overplot a line for the dark rate
  
  xs_dk = DINDGEN(11930, START=70.0d)
  ys_dk=darkRate[0]+darkRate[1]*xs_dk
  
  plDk2 = PLOT(xs_dk,ys_dk,'b--',THICK=2,NAME='D$_{R}$ (Fit)',/OVERPLOT)
  
  ; Add a legend

  leg_dk = LEGEND(TARGET=[plDk,plDk2], POSITION=[50,10000],/DATA)
   
             
  STOP




END
