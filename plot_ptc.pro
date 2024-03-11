  ;+
  ; :Author: bob.zavala
  ;-
  
  ;**************************************************************
  ;
  ; CODE FOR EXPERIMENTAL PHOTON TRANSFER CURVE ANALYSIS
  ;
  ; Based on Chapters 1-8 and 11 of 
  ;       "Photon Transfer DN -> Lambda" by James R. Janesick
  ; 
  ; Using measured data of the signal in DN from a digital scientific 
  ; camera analyse the performance of the camera. Plots of the 
  ; standard deviation versus Signal level are the Photon Transfer 
  ; curve plots. These plots and data are used to determine camera 
  ; performance such as:
  ; 
  ; Full well :           DN
  ; Read noise  :         DN 
  ; Sensitivity (gain):   electrons per DN
  ; Dark current :        electrons per second
  ; Fixed Pattern Noise : dimensionless 
  ;      Quality factor
  ; Check linearity, identify sources of non-linearity
  ;      
  ; CAMERA(S) EVALUATED:
  ; 
  ; As of 08 Mar 2024 the #### Princeton Instruments (PI) Sophia
  ; 4k 4096B CCD serial # # performance is analysed by 
  ; this code.
  ; 
  ; RESULT VERIFICATION:
  ; 
  ; Table below shows results for linearity data obtained 
  ; on the # with the Sophia camera on # (UT date ####). 
  ; Camera set to -60 degrees C and 
  ; light levels provided by flat field lamps (set to X V) 
  ; and V filter. 
  ;
  ; Dark data were taken on multiple nights and results analysed in 
  ;   
  ;   images/dark-dtc 
  ; 
  ;               PI specification              PTC result
  ;               ----------------              ----------
  ; Dark current  0.8 e^-/sec @-60C             0.869 +/- 0.003 e^-/sec
  ; Gain (K)      1 e^-/DN                      1.024 e^-/DN
  ; read noise    4 e^- rms                     4.026 e^-
  ; Full well     350 kilo e^- (typical)        63174.0 DN = 64690 e^- 
  ;                                             (K as measured and set)
  ; 
  ; Teledyne e2v dark current figure of merit D_m (nanoAmp/cm^2 at 300K)
  ; predicted using documentation A1A-765136 Version 8, August, 2018
  ; = 1.8 nA/cm^2 and measured at 173K and "... some variation may 
  ; be seen between devies" (pg. 2 of documentation).
  ;   
  ;   Dark Transfer Curve (DTC) measured D_m: 1.177 nA/cm^2
  ; 
  ; Stability of flat field lamps on #:
  ;     10.0 second exposures in V lamps at X V did show occasional
  ;     signal drops of 4%. Evaluation of linearity and sources of 
  ;     any non-linearity requires a stable laboratory environment. 
  ;     Data taken prior to fully reading chapters in Janesick and 
  ;     FPN noise removal with image differencing requires image pairs
  ;     taken at signal levels > 4000.0 DN. (See section below 
  ;     entitled "CODE STILL TO BE WRITTEN")
  ;     
  ; RUNNING THE CODE:
  ; 
  ; The code expects an ASCII text file 'raw.signal.txt' of PTC data 
  ; extracted from FITS files. The columns are expected to be with 
  ; a first row containing this column header information:
  ; 
  ; filename statsec npix mean median stddev min max exptime utc-obs obstype
  ; 
  ; First row format is not important, as a "#" in the first row will cause
  ; this row to be ignored. 
  ; 
  ; filename: the FITS image filename
  ; statsec:  region in the images from where the statistics are drawn
  ; npix:     number of pixels in region statsec
  ; mean:     the mean pixel value in statsec
  ; median:   median pixel value in statsec
  ; stddev:   standard deviation of pixels in statsec
  ; min:      minimum pixel value in statsec
  ; max:      maximum pixel value in statsec
  ; exptime:  exposure time in seconds
  ; utc-obs:  UTC time of exposure as hh:mm:ss.sss
  ; obstype:  type of observation (BIAS, DARK, DOME_FLAT)
  ; 
  ; NOTE: For this PTC analysis the above ASCII data are expected to 
  ; have the overscan substracted. This was conveniently done prior 
  ; to reading in the ASCII data as the PI Sophia camera does not have 
  ; an overscan region. It was easier to compute the "offset" via a Zero
  ; image and subtract this AND THEN get the image statistics.  
  ;                      
  ; The code to perform the statistics is a separate script and may be 
  ; a user preference. There is also some Zero image analysis needed, 
  ; which is a separate task and depends on the presence/absence of an 
  ; overscan region.
  ; 
  ; CODE RUNNING STEPS:
  ; 
  ; Assuming the ASCII file described above is present in a directory 
  ; containing one night of data
  ; 
  ;         somePath/images/g##d###/
  ;         
  ;   1) Run plot_ptc: Select the raw.signal.txt file. Results will 
  ;   be plotted to the screen. If darks were taken do_darks will 
  ;   be called.  
  ;   
  ;   2) If you will analyse many nights of dark data together copy 
  ;   the *.darkData.sav file to:
  ;   
  ;         somePath/images/dark-dtc/ 
  ;         
  ;   or a suitable named directory.
  ;   
  ;   3) Dark Transfer Curve (DTC) of many nights in the dark data 
  ;   containing directory is performed with do_dtc. Results again 
  ;   plotted to the screen and plots amy be saved. 
  ;   
  ;   4) Finish PTC analysis and linearity analysis in a single 
  ;   night directory with linear_ptc procedure. 
  ;   
  ; CODE STILL TO BE WRITTEN:
  ; 
  ; At the time the linearity data were collected I did not thoroughly 
  ; understand the PTC analysis. SUfficient data is required to remove 
  ; the Fixed Pattern Niose (FPN) from all images for the entire dynamic 
  ; range of the CCD. Another round of data collection is required and 
  ; code written to complete the final linearity checks and tests for 
  ; any sources of non-linearity that are found.  
  ;     
  ;**************************************************************
  
        
  ;**************************************************************      
  ;      
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
  imagesDir = '/images/'
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

  ; Make an ASCIII template to read the statistics files created by 
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
     
    x_min_p1 = 0.1 
     
    p1 = PLOT(xs, ys, XLOG=1, YLOG=1, XRANGE=[x_min_p1,1e6], YRANGE=[1,3*MAX(statdata.sigma)], $
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
    
    ; If long exposures that are not averaged (1.5 seconds) or more are present 
    ; plot these
      
    ; stop  
      
   
    ; Determine and overplot the sigma_readnoise in DN
    
    the_biases = WHERE(statData.obstype EQ 'BIAS')
    biasSignal = [MEAN(statData.median[the_biases]),MEAN(statData.median[the_biases])]
    biasNoise  = [MEAN(statData.sigma[the_biases]),MEAN(statData.sigma[the_biases])]
    
    name_rd = '$\sigma_{READ} = '+STRTRIM(STRING(biasNoise[0],FORMAT='(F5.3)',/PRINT),1)+' DN$'   
    
   IF x_min_p1 GT biasSignal[0] THEN BEGIN
     bias_Signal_plot = [x_min_p1,x_min_p1]
   ENDIF ELSE BEGIN
     bias_Signal_plot = biasSignal
   ENDELSE
     
   p_rd = PLOT(bias_Signal_plot,biasNoise,'sb2',NAME=name_rd,SYM_FILLED=1,/OVERPLOT)

   IF long_lin_index[0] NE -1 THEN BEGIN

     ; Overplot the single exposure raw signal and noise PTC values
     p3 = PLOT(long_signal, long_noise,'2D',NAME='Single exp. (DN)',SYM_FILLED=0,/OVERPLOT)

     ; Now add the legend through p4

     leg_p3 = LEGEND(TARGET=[p1,p2,p3,p_rd], POSITION=[800,1000],/DATA,/AUTO_TEXT_COLOR)

   ENDIF ELSE BEGIN

     ; No single exposure data, so plot a legend without p3

     leg_p2 = LEGEND(TARGET=[p1,p2,p_rd], POSITION=[800,1000],/DATA,/AUTO_TEXT_COLOR)

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

  ; stop

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
  
  IF maxYrange GT 100.0 THEN BEGIN
    
    leg_p1 = LEGEND(TARGET=[p1,p2], POSITION=[50,100],/DATA)
  
  ENDIF ELSE BEGIN
     
    leg_p1 = LEGEND(TARGET=[p1,p2], POSITION=[3000,4],/DATA)
    
  ENDELSE
    
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

FUNCTION dfigm,ccd_temp,pix_a,dkrate

; A function used with PROCEDURE do_dtc to calculate 
; and RETURN the valure of the Dark Figure of Merit
; in nanoAmperes per centimeters^{2}. 
; 
; INPUT:
;
; ccd_temp  = The CCD oeprating temperature in degrees Kelvin
; pix_a     = Area of a pixel in square centimeters
; dkrate    = Fitted dark rate in electrons per second

  ; Set Boltzmann's constant in eV per degree Kelvin
  
  k_Boltz = 8.62e-5

  ; Calculate the energy band gap (eV) using Equation 11.17  
  ; from Janesick's "Photon Transfer: DN -> Lambda"
  
  e_gap = 1.1557 - 7.021e-4*(ccd_temp)^2/(1108.0 + ccd_temp)
  
  ; Now get the Dark Figure of Merit using Equation
  ; 11.21 from Janesick
  
  denom = (2.55e15*pix_a*(ccd_temp)^1.5)*exp(-1.0*e_gap/(2*k_Boltz*ccd_temp))
  
  ;print,' '
  ;print,' The denominator =: ',STRTRIM(STRING(denom),1)
  ;print,' '
  
  dark_f_m = dkrate/denom
  
  print,''
  print,'The Dark Figure of Merit is: ',STRTRIM(STRING(dark_f_m),1),' nA per cm^2'
  print,'at 300 degrees Kelvin. '
  print,' '
  
  RETURN,dark_f_m

END

;******************************

PRO do_dtc

  ptc
  
  COMMON SHARED,imagesDir

  ; In here I will plot the readnoise signal and sigma in DN determined on 
  ; # from 200 bias frames using:
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
  ; This records that the first "theSkip" indices are already filled. 
  
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

  D_N = TEXT(0.5,430,D_N_text,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')
  
  ; As the "de-scaling" looks like it works well for Signal >= 2900 DN 
  ; use this threshold for detemring K_adc and get the mean. This is 
  ; well behaved so uncertainty is negligible
  
  Sig_29hun = WHERE(allDarkData.median GE 2900.0)
  
  ; Determine K_adc (electrons/DN) using Equation 11.20
  
  K_29hun = allDarkData.median[goodSignal[Sig_29hun]]/(deScaled_pois[goodSignal[Sig_29hun]])^2
  
  K = MEAN(K_29hun)
  
  ; Label K_adc
  
  ; Label the Dark Current Quality Factor

  K_adc_text = 'K$_{ADC}$ = '+STRTRIM(STRING(K,FORMAT='(F5.3)',/PRINT),1)+' $(e^{-}/DN)$'

  K_adc = TEXT(0.5,300,K_adc_text,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')
  
  ; Now that we know the ADC sensitivity (gain in amateur speak) fit for the 
  ; dark signal versus time
  
  D_etron = K*allDarkData.median
  D_etron_sig = SQRT(D_etron)
  
  ; Linear fit for the dark rate against time for exposure times 
  ; greater than 200 seconds. These are well and scatter increases
  ; below that exposure time
  
  times_to_fit = WHERE(allDarkData.exptime GE 200.0)
  
  darkRate = LINFIT(allDarkData.exptime[times_to_fit],D_etron[times_to_fit], $ 
                   MEASURE_ERRORS=D_etron_sig[times_to_fit],$ 
                   CHISQR=drate_chisq,COVAR=drate_covar,SIGMA=drate_sigma,$
                   YFIT=drat_yfit)
  
  ; Now to calculate the camera Dark figure of merit in nano Amperes cm^(-2)
  ; Right now this is the Sophia camera, remind the user and press on.
  msg_text = 'Parameters will be set for the Sophia CCD Dark Figure of Merit'
  df_msg = DIALOG_MESSAGE(msg_text,/INFORMATION,TITLE='D_F MESSAGE',/CENTER)
  DELVAR,df_msg

  ; Sophia temperature = -60C = 213 K
  ccd_temp = 213.0  ; deg K
  pix_a = 2.25e-6   ; cm^2 from 15 micron pixels

  dark_fig_merit = dfigm(ccd_temp,pix_a,darkRate[1])

  ; Label the Dark Figure of Merit

  D_fig_merit_text = 'D$_{M}$ = '+STRTRIM(STRING(dark_fig_merit,FORMAT='(F5.3)',/PRINT),1)+' $(nA/cm^{2})$'

  D_fig_merit = TEXT(0.5,200,D_fig_merit_text,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')
  
  ; Label the read noise in electrons
  
  sigma_rd_etrons = rd_noise_sigma[0]*K
  
  sigma_rd_text = '$\sigma_{READ} = $'+STRTRIM(STRING(sigma_rd_etrons,FORMAT='(F5.3)',/PRINT),1)+' $e^{-}$'
  
  sigma_rd = TEXT(0.5,140,sigma_rd_text,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')
  
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

  D_rate_text = 'D$_{R}$ = '+STRTRIM(STRING(darkRate[1],FORMAT='(F5.3)',/PRINT),1)
  D_rate_text2= '    $\pm$ '+STRTRIM(STRING(drate_sigma[1],FORMAT='(F5.3)',/PRINT),1)+' $(e^{-}/sec)$' 
  D_rate = TEXT(2.5,1000,D_rate_text,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')
  D_rate2= TEXT(2.5, 670,D_rate_text2,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')

  ; Overplot a line for the dark rate
  
  xs_dk = DINDGEN(11700, START=200.0d)
  ys_dk=darkRate[0]+darkRate[1]*xs_dk
  
  plDk2 = PLOT(xs_dk,ys_dk,'b--',THICK=2,NAME='D$_{R}$ (Fit)',/OVERPLOT)
  
  ; Add a legend

  leg_dk = LEGEND(TARGET=[plDk,plDk2], POSITION=[50,9000],/DATA)

END

;******************************

FUNCTION map_fpn_diff, shotReadData,fpn_diffData,yearDay

COMMON SHARED,imagesDir,rawFilter,oneImstatLine

; MAPPING function used with class LIST theMults 
; IN DEVELOPMENT
; 
; INPUT:
;
; z  = Integer array containing the location indices 
;      in the structure statData where the same exposure 
;      time image statistics are to used to compute the 
;      FPN differences
;
; shotReadData = structure holding signal, total noise, 
;            in DN and filename for the images that will 
;            be used to produce the shot+read and shot noise
;            curves. Also holds the readnoise. Data selected 
;            from shotData structure. 
;            
; fpn_diffData = a structure containing the IMSTAT 
;            results for the FPN difference images. 
; 
; OUTPUT:
; 
; shotReadData = Strucutre described above is returned 
;            with data entered for delta_noise, sigma_(SHOT)
;            and sigma_(SHOT+READ)
; 
 
  print,' '
  print,' In FUNCTION map_fpn_diff'
  print,' '
  
  ; Print the elements in a mult in theMults
  ; print,z,':'
  ; print,' '

  ; Create a LIST of the image number pairs that go into 
  ; creating the differences for subtracting out the FPN noise.
  ; I do this within the function as passing a LIST within a List::Map
  ; method iterates over the elements in the List. Recall this only 
  ; passes one element in a list to the function. If fpn_DiffData.filename
  ; contains n elements theDiffPairs is a two column LIST indexed as:
  ; 
  ; theDiffPairs[n,2] as n_columns = 2
  
  theDiffPairs = STRSPLIT(STRMID(fpn_diffData.FILENAME[*],8,7),'-',/EXTRACT)
  
  ; Loop through theDiffPairs and fil in the values for the required 
  ; entries in shotReadData
  
  FOR k = 0, n_elements(theDiffPairs)-1 DO BEGIN
    furstOne = yearDay+'.'+theDiffPairs[k,0]+'.fits'
    secundOne = yearDay+'.'+theDiffPairs[k,1]+'.fits'
    furst_loc = WHERE(shotReadData.filename EQ furstOne)
    secund_loc = WHERE(shotReadData.filename EQ secundOne)
    shotReadData.sig_r_shot[furst_loc] = fpn_diffData.sigma[k]/SQRT(2)
    shotReadData.sig_r_shot[secund_loc] = shotReadData.sig_r_shot[furst_loc] 
    shotReadData.sig_shot[furst_loc] = SQRT(shotReadData.sig_r_shot[furst_loc]^2 - shotReadData.readnoise^2)
    shotReadData.sig_shot[secund_loc] = SQRT(shotReadData.sig_r_shot[secund_loc]^2 - shotReadData.readnoise^2)
  ENDFOR


  ; STOP 
  ; Return 1 for Success
  ; When ready, RETURN,shotReadData
  RETURN,shotReadData

END

;******************************

PRO linear_ptc

ptc 

COMMON SHARED,imagesDir,rawFilter,oneImstatLine

; Procedure to conduct linarity PTC analysis as in Chapter 5 of 
; Janesick's "Photon Transfer DN -> Lambda"

; USES
; 
; IDL SAVE file created for a given g??d??? night by the 
; procedure plot_ptc and previously computed dark current results 
; in the ../images/dark-dtc directory. The path assusmes the user 
; is in an images/g??d???/ directory.

  ; Select the IDL save file with the linearity data
  statData = DIALOG_PICKFILE(FILTER='*lin_dat*')
  RESTORE,FILENAME=statData

  ; I know for g23d314 by visual inspection that the four PTC 
  ; pairs (signal and noise) before saturation are
  ;  on the linear FPN curve and Eqn. 5.8 
  ; of Janesick gives P_n the FPN quality factor. 
  
  P_n_values = long_noise[12:15]/long_signal[12:15]
  
  P_n = MEAN(P_n_values)
  
  ; Is this g23d314? If so I do not have the pairs of equal light exposure 
  ; images across the entire PTC curve. Some manual handling is necessary.
  
  IF yearDay EQ 'g23d314' THEN BEGIN
 
    ; To determine the sigma_(SHOT+READ) pairs of exposures at the same 
    ; light levels are needed. I have these at the low light level, shorter 
    ; exposure times for #. At the longer exposures I do not as 
    ; these were single exposures. An approximation on # is possible 
    ; for images ????314 and ????316 as the exposure times were 150.0
    ; and 160.0 seconds 6.3% different in exposure times. This was made into 
    ; the image test-FPN-diff.fits with the statistice recorded to:
    ; test-FPN-diff.imstat.txt.
    ; 
    ; Average the signal levels and record this to:
    ; test_fpn_diff variables manually. Note that the SIGNAL for 
    ; plotting is the AVERAGE signal of the two images input to the 
    ; difference image. 
    ;
    ; Reduce the noise in the difference image by a factor of SQRT(2)
    ; as the difference image over-estimates the noise.   
    
    average_input_sig = (63174.0+59407.0)/2.0
    test_fpn_diff_noise = 301.9/SQRT(2)
    
    ; Calculate this "high end" signal and noise of sigma_(READ+SHOT)
    ; using Equation 5.11 of Janesick
    
    test_sigma_rd_shot = SQRT(test_fpn_diff_noise^2 - biasnoise[0]^2)
    
  ENDIF 

  ; Now do the diferences for the equal exposure times. These are SO FAR
  ; 0.2, 0.4, 0.8, 1.0 and 10.0 second exposures. Recall that the 
  ; 0.1 second exposures are too sensitive to the variable lamp fluxes.
  
  ; the locations of the exposures
  theTwos = WHERE(statData.exptime EQ 0.2)
  theFours = WHERE(statData.exptime EQ 0.4)
  theEights = WHERE(statData.exptime EQ 0.8)
  theOnes = WHERE(statData.exptime EQ 1.0)
  theTens = WHERE(statData.exptime EQ 10.0)
  
  ; All together in order
  
  allToDiff = [theTwos,theFours,theEights,theOnes,theTens]
  
  ; Place these arrays containing index locations into an array because 
  ; I will go through a class LIST called theMults for "the multiples"
  
  theMults = LIST(theTwos,theFours,theEights,theOnes,theTens)
  
  ; Read in the FPN difference IMSTAT data.
  ; *NOTE* the use of a generic IMSTAT 7 column template
  ;        and I expect to always name give the FPN 
  ;        difference files the same names as they aer
  ;        in unique directories. 
  
  fpnDiffStatFile = 'all_fpnDiffs.imstat.txt'
  RESTORE,FILENAME=imagesDir+'imstat.template.sav'
  fpn_diffData = READ_ASCII(fpnDiffStatFile,TEMPLATE=imstat_templ)
  
  ; Make a structure to hold the sigma_(SHOT+READ) and sigma_(SHOT)
  ; results obtained using the difference data
 
  array_size = 2*n_elements(fpn_diffData.FILENAME)
 
  shotReadData = {readNoise:biasnoise[0], filename:STRARR(array_size), $ 
                 im_num:STRARR(array_size), signal:FLTARR(array_size), $
                 total_noise:FLTARR(array_size), sig_fpn:FLTARR(array_size), $ 
                 sig_r_shot:FLTARR(array_size), sig_shot:FLTARR(array_size)}
                 
  shotReadData.filename[*] = statData.FILENAME[allToDiff]
  shotReadData.signal[*] = statData.median[allToDiff]
  shotReadData.total_noise[*] = statData.sigma[allToDiff]
  shotReadData.im_num[*] = STRMID(statData.FILENAME[allToDiff],8,3)

  ; Use the function map_fpn_diff to subtract the FPN noise, fill in 
  ; sigma(read + shot) and sigma(shot)

  shotReadData = map_fpn_diff(shotReadData,fpn_diffData,yearDay)

  ; Locate the full well in DN
  
  donde_fullWell = WHERE(statDAta.sigma[domeflats[long_lin_index]] $ 
                  EQ MAX(statDAta.sigma[domeflats[long_lin_index]]))
                  
  fullWell_DN = statData.median[domeflats[long_lin_index[donde_fullWell]]]

  ; Plot the PTC, replotting with the code from plot_ptc
  
  x_min_p1 = 0.1

  p1 = PLOT(xs, ys, XLOG=1, YLOG=1, XRANGE=[x_min_p1,1e6], YRANGE=[1,3*MAX(statdata.sigma)], $
    XTITLE='Signal (DN)', YTITLE='NOISE (DN)', FONT_SIZE=12, $
    FONT_STYLE='Bold', FONT_NAME='Times', $
    Title = yearDay+' PTC curve for PI Sophia CCD Camera', $
    XTHICK=2, YTHICK=2, /NODATA)

  ; Begin plotting the raw signal and noise
  p2 = PLOT(short_signal,short_noise,'2D',NAME='Averaged exp. (DN)',SYM_FILLED=1,/OVERPLOT)

  p_rd = PLOT(bias_Signal_plot,biasNoise,'sr2',NAME=name_rd,SYM_FILLED=1,/OVERPLOT)

  ; Overplot the single exposure raw signal and noise PTC values
  p3 = PLOT(long_signal, long_noise,'2D',NAME='Single exp. (DN)',SYM_FILLED=0,/OVERPLOT)

  p_shot = PLOT(shotReadData.signal,shotReadData.sig_shot,'2gX',NAME='$\sigma_{SHOT}$ (DN)', /OVERPLOT)  

  ; Overplot a Slope = 1 FPN line for a Pn = 0.01
  
  Pn_x = [100.0,300000.0] 
  Pn_y = [1.0,3000.0]
  
  ; Label the FPN Quality Factor
  
  p_FPN = PLOT(Pn_x,Pn_y,'--',NAME='$\sigma_{FPN}$',/OVERPLOT)
  P_n_text = 'P$_{n}$ = 0.01'
  P_n_show = TEXT(300,2,P_n_text,/DATA,FONT_SIZE=12,FONT_STYLE='BOLD')
  
  ; Label the Full Well in DN
  
  fullWell_text = 'S$_{FW} = $'+STRTRIM(STRING(fullWell_DN,FORMAT='(F7.1)',/PRINT),1)+' (DN)'
  fW_text_show = TEXT(8000.0,1100.0,fullWell_text,/DATA,FONT_SIZE=10,FONT_STYLE='BOLD')
  
  ; Overplot a Slope = 1/2 shot noise
  
  shot_line_x = [1,1e6]
  shot_line_y = [1,1000.0]
  
  p_shot_line = PLOT(shot_line_x,shot_line_y,NAME='$\sigma_{shot}$',/OVERPLOT)
  
  ; Now add the legend through p4

  leg_p_shot = LEGEND(TARGET=[p1,p2,p3,p_rd,p_shot,p_FPN,p_shot_line], $ 
               POSITION=[400,1000],/DATA,/AUTO_TEXT_COLOR)


  ;STOP
END