; Examples worked out from "Photon Transfer: DN -> Lambda" by James R. Janesick
;

;**********
; Chapter 2
;**********

PRO ex_2_3

; Example 2.3 of J. Janesick "Photon Transfer DN -> Lambda"

COMPILE_OPT IDL2, hidden

; prepare ASCII template to read the data

; templ = ASCII_TEMPLATE('ex_2.3.txt')

; SAVE,templ,filename='ex_2_3.templ.sav'

; Read the data

restore,filename='ex_2_3.templ.sav'

exam_data = READ_ASCII('ex_2.3.txt',TEMPLATE=templ)

; Set the optical active thickness T_epi in microns

T_epi_1 = 10.0

T_epi_2 = 25.0

; Now calculate the interacting QE with and without AR coating

QE_non_ar_1 = (1 - exam_data.reflCoeff)*(1.0 - exp(-1.0*T_epi_1/exam_data.aLen))

QE_with_ar_1 = exam_data.ar_factor*(1.0 - exp(-1.0*T_epi_1/exam_data.aLen))

QE_non_ar_2 = (1 - exam_data.reflCoeff)*(1.0 - exp(-1.0*T_epi_2/exam_data.aLen))

QE_with_ar_2 = exam_data.ar_factor*(1.0 - exp(-1.0*T_epi_2/exam_data.aLen))


; Plot the results with non-AR first in blue

pl10_1 = PLOT(exam_data.lam,QE_non_ar_1,'b-',xtitle='Wavelength (nm)', $ 
           ytitle = 'Interacting QE',xrange=[300,1000],yrange=[0,1.0], $
           FONT_SIZE=12)
      
pl10_2  = PLOT(exam_data.lam,QE_non_ar_2,'b-.',/OVERPLOT)    

; Now plot the AR coated

pl10_3 = PLOT(exam_data.lam,QE_with_ar_1,'r-',/OVERPLOT)

pl10_4  = PLOT(exam_data.lam,QE_with_ar_2,'r-.',/OVERPLOT)

           
END

;**********
; Chapter 3
;**********

FUNCTION my_gauss,x,the_mean,sig
  
  ; My own personal Gaussian distribution calculated for input x
  ; with mean the_mean and standard deviation sig
  
  ; The return value is normalized

  A = 1/(sig*SQRT(2*!pi))
  B = exp(-0.5*((x-the_mean)/sig)^2)
  
  the_gauss = A*B
  norm_factor = MAX(the_gauss)
  
  RETURN, the_gauss/norm_factor

END

; Example 3.2 of J. Janesick "Photon Transfer DN -> Lambda"

; Example 3.2

PRO ex_3_2

COMPILE_OPT IDL2, hidden

; How many discrete Poisson events do we want to simulate?
; Set an integer to this number and then create an increasing 
; floating point array to hold these input integers

print,' '
print,'**********'
print,' Begin Example 3.2 Photon Transfer DN -> Lambda '
print,'**********'
print,' '

do_poisson = 11

poiss_arr = FINDGEN(do_poisson)

; Now make an array to hold the Gaussian niose distribution on top
; of the discrete Poisson probability distribution

rand_noise = [0.1,0.2,0.3,0.4,0.5]


; signal_count is an integer to record the number of 
; electrons from 0 to signal_count that we will plot in this procedure
; and it equals the number of elements in the rand_noise array

signal_count = n_elements(rand_noise)

; Assuming the quantum yield gain eta is 1 electron/interacting photon
; what is the average interacting photon per pixel?

P_i = 1.0

; Now calculate the Poisson probabilities for numbers of photons interacting
; with a pixel

pois_prob = ((P_i)^poiss_arr/FACTORIAL(poiss_arr))*exp(-1*P_i)


; g1=my_gauss(3.0,3.0,1.7321)

; print,' '
; print,'g1 = ',g1
; print,' '

; Make an array of x values that will be input to the normalized
; Gaussian to produce the random noise values that are applied to 
; the Poisson noise values.
; This normalized Gaussian spans 0.0 <= x <= 2.0 in increments 
; of 0.1 of x.

xs = FINDGEN(21,INCREMENT=0.1,START=0.0)

; Create a single structure to hold the noise results to plot for one 
; random noise level

one_noise = {nnn:FLTARR(10*(signal_count) + 11)}

; Now replicate this structure signal_count + 1 times to hold all the 
; noise results as many_noise[i].nnn

many_noise=REPLICATE(one_noise,signal_count)

; Loop over random noise array. 

FOR i = 0, n_elements(rand_noise) - 1 DO BEGIN
    
    ; Now loop over the signal count. Note that signal_count
    ; begins at 0 so I do not need a n_elements( ) minus 1
    ; in the FOR loop.
    FOR j = 0, signal_count - 1 DO BEGIN
       
    ; Create the Gaussian random noise values
    g_noise = pois_prob[j]*my_gauss(xs,P_i,rand_noise[i])
    
    ; here's how far we shift over in x-coordinates
    the_shift = (j * 10)
    
    ; Enter the random noise numbers into the structure
    
    ;print,'i = ',STRTRIM(STRING(i),1)
    ;print,'j = ',STRTRIM(STRING(j),1)
    ;print,'coords for xs are: [',STRTRIM(STRING(0+the_shift),1),':',STRTRIM(STRING(20+the_shift),1),']'
    ;print,' '
    
    many_noise[i].nnn[0+the_shift:20+the_shift] = g_noise + many_noise[i].nnn[0+the_shift:20+the_shift]
     
    ENDFOR
    
ENDFOR

title_0 = STRMID(STRING(rand_noise[0]),5,3)
title_1 = STRMID(STRING(rand_noise[1]),5,3)
title_2 = STRMID(STRING(rand_noise[2]),5,3)
title_3 = STRMID(STRING(rand_noise[3]),5,3)
title_4 = STRMID(STRING(rand_noise[4]),5,3)

; Plot the results as per Figure of the text. 

; Create a floating point array to nicely display he electron signal
; from 0 -> 6.0 electrons

sig_etrons = many_noise[0].nnn/10

;inx = FINDGEN(10)
;iny = FINDGEN(10)

;bb0 = BARPLOT(inx, iny,XRANGE=[-1,5],YRANGE=[0,0.4], YTITLE = 'Occurences', $ 
;         XTITLE='Signal, $e^{-}$', $
;         title='Random noise = '+title_0+' e$^{-}$',YMINOR=0,XMINOR=0,FONT_SIZE=12,/NODATA)

;b0 = BARPLOT(many_noise[0].nnn,FILL_COLOR='blue',/CURRENT,AXIS_STYLE=0)
; b0.XSHOWTEXT=1
; b0.XTICKN=['0','1','2','3','4','5','6']
     ; XTICKN=['0','1','2','3','4','5','6'])

b0 = BARPLOT(many_noise[0].nnn,FILL_COLOR='blue',title='Random noise = '+title_1+' e$^{-}$', $
            YTITLE = 'Occurences',XTITLE='Signal, $e^{-} \times$ 10',YMINOR=0,XMINOR=0, $ 
            FONT_SIZE=12)
                 
b1 = BARPLOT(many_noise[1].nnn,FILL_COLOR='blue',title='Random noise = '+title_1+' e$^{-}$', $
     YTITLE = 'Occurences',XTITLE='Signal, $e^{-} \times$ 10',YMINOR=0,XMINOR=0,FONT_SIZE=12)
     
b2 = BARPLOT(many_noise[2].nnn,FILL_COLOR='blue',title='Random noise = '+title_2+' e$^{-}$', $
     YTITLE = 'Occurences',XTITLE='Signal, $e^{-} \times$ 10',YMINOR=0,XMINOR=0,FONT_SIZE=12)

b3 = BARPLOT(many_noise[3].nnn,FILL_COLOR='blue',title='Random noise = '+title_3+' e$^{-}$', $
     YTITLE = 'Occurences',XTITLE='Signal, $e^{-} \times$ 10',YMINOR=0,XMINOR=0,FONT_SIZE=12)
     
b4 = BARPLOT(many_noise[4].nnn,FILL_COLOR='blue',title='Random noise = '+title_4+' e$^{-}$', $
     YTITLE = 'Occurences',XTITLE='Signal, $e^{-} \times$ 10',YMINOR=0,XMINOR=0,FONT_SIZE=12)

print,' '
print,'**********'
print,' Example 3.2 Photon Transfer DN -> Lambda complete'
print,'**********'
print,' '

;stop

END

;**********

;**********
; Chapter 5
;**********

;****************************************************

FUNCTION my_wien,input_x,input_y,turnover

; A function I created to model in a coarse manner the turnover that 
; occurs when pixels reach full well. The turnover is too sharp, but is 
; good enough for my purposes. I use the Wien Law as a basis. 
; 
; THE PARAMETERS ARE:
;
; input_x:  an array of input signal values (DN or electrons)
; input_y:  an array of the pixel values in DN or electrons for the 
;           Y axis
; turnover: the pixel full well value at which the turnover occurs
; 
; RETURNS:
; 
; wien_values: an array of pixel values calculated using the Wien Law 
;              as a proxy
; 

  wien_values = 2.57*input_y*(exp(-1.0*input_x/turnover))
  
  RETURN, wien_values 

END 

;****************************************************

PRO ex_5_1

  COMPILE_OPT IDL2, hidden

  ; Generate Photon Transfer Curves (PTC's) for sigma_total,  
  ; sigma_read, sigma_shot and sigma_FPN given the parameters 
  ; listed in Example 5.1 of Janesick's  
  ; "Photon Transfer DN -> Lambda".
  
  ; Slopes are for a log-log plot.
  
  ; sigma_read = constant = 3.33 DN
  
  ; maximum signal at saturation (ordinate) = 233,000 = 10^5.3674
  
  ; maximum sigma_FPN (slope = 1) = P_n*(Full well) = 
  ;                               = 0.02*233,000 = 4,666 = 10^3.669
  
  ; maximum sigma_shot (slope = 1/2) = SQRT(233,000) = 482.7
  ;                                  = 10^2.6837
                        
  ; Set up a blank log-log plot, then overplot the data            
  
  xs = DINDGEN(400000, START = 1.0d)
  ; let ys all = 1.0
  ys = xs/xs 
  
  p1 = PLOT(xs, ys, XLOG=1, YLOG=1, XRANGE=[1,1e6], YRANGE=[1,1e4], $ 
           XTITLE='Signal (DN)', YTITLE='NOISE (DN)', FONT_SIZE=12, $ 
           FONT_STYLE='Bold', FONT_NAME='Times', Title = 'Example 5.1 DN $\rightarrow$ $\lambda$', $
           XTHICK=2, YTHICK=2, /NODATA)
           
  ; Plot the read noise
  
  sigma_read = 3.33*ys         ; 
  
  p2 = PLOT(xs, sigma_read, NAME='$\sigma_{RN}(DN)$', THICK=2, /OVERPLOT)
 
  ; Plot the shot noise
  ; Minimum signal is 1 electron generates 1.5 DN
  
  K_adc = 1.5 

  sigma_shot = SQRT(xs/K_adc)  ; Eqn. 5.7 = Eqn. 4.20 
  
  ; Modify for full well
  
  ; fw = WHERE(xs GE 233000)
  fw = WHERE(xs GE 220000)
  
  ; This gives a wrong shape
  ; sigma_shot[fw] = sigma_shot[fw] - SQRT((xs[fw]-233000.0)/K_adc) 

  ; Better, but too extreme of a turnover
  ; sigma_shot[fw] = sigma_shot[fw] - exp(xs[fw] - 233000)
  
  ; Try a function based on Planck Blackbody law
  ; sigma_shot[fw] = sigma_shot[fw]/(exp(xs[fw]/233000) - 1)
  
  ; Try a function based on Wien's law
  ; sigma_shot[fw] = 2.57*sigma_shot[fw]*(exp(-1.0*xs[fw]/233000))
  sigma_shot[fw] = my_wien(xs[fw],sigma_shot[fw],233000)

  p3 = PLOT(xs, sigma_shot, NAME='$\sigma_{SHOT}(DN)$', THICK=2, LINESTYLE=2, $ 
           /OVERPLOT)

 ; Plot the FPN noise

 P_n = 0.02
 
 sigma_fpn = P_n*xs            ; Eqn. 5.8 = Eqn. 3.12

 ; Allow for full well FPN
 sigma_fpn[fw] = my_wien(xs[fw],sigma_fpn[fw],233000)

 p4 = PLOT(xs, sigma_fpn, NAME='$\sigma_{FPN}(DN)$', THICK=2, LINESTYLE=1, $ 
          /OVERPLOT)

 ; sigma_total is the sum of Read Noise, Shot Noise and FPN Noise
 
 ; sigma_total = sigma_read + sigma_shot + sigma_fpn NO NO NO Quadrature!!!
 
 ; Eqn. 5.9
 sigma_total = SQRT(sigma_read^2 + sigma_shot^2 + sigma_fpn^2)
 
 p5 = PLOT(xs, sigma_total, NAME='$\sigma_{TOTAL}(DN)$', THICK=2, LINESTYLE=4, $
          /OVERPLOT)
          
 ; Add the LEGEND to the plot
 
 leg = LEGEND(TARGET=[p2,p3,p4,p5], POSITION=[40,4000], /DATA)
 
 ; Now re-make the plots in terms of electrons
 
 xs_etrons = xs*K_adc
 
 p11 = PLOT(xs, ys, XLOG=1, YLOG=1, XRANGE=[1,1e6], YRANGE=[1,1e4], $
   XTITLE='Signal (e$^-$)', YTITLE='NOISE (e$^-$)', FONT_SIZE=12, $
   FONT_STYLE='Bold', FONT_NAME='Times', Title = 'Example 5.1 DN $\rightarrow$ $\lambda$', $
   XTHICK=2, YTHICK=2, /NODATA)

 ; Read noise in electrons
 
 sigma_read_etrons = 3.33*ys*K_adc         ;

 p22 = PLOT(xs_etrons, sigma_read_etrons, NAME='$\sigma_{RN}(e^{-})$', THICK=2, /OVERPLOT)
 
 ; Plot the Shot noise in electrons
 
 sigma_shot_etrons = sigma_shot*K_adc
 
 p33 = PLOT(xs_etrons, sigma_shot_etrons, NAME='$\sigma_{SHOT}(e^{-})$', THICK=2, $ 
           LINESTYLE=2, /OVERPLOT)
           
 ; Plot the FPN noise in electrons
 
 sigma_fpn_etrons = sigma_fpn*K_adc
 
 p44 = PLOT(xs_etrons, sigma_fpn_etrons, NAME='$\sigma_{FPN}(e^{-})$', THICK=2, $ 
          LINESTYLE=1, /OVERPLOT)
          
 ; Plot sigma_total in electrons 
 
 sigma_total_etrons = sigma_total*K_adc
 
 p55 = PLOT(xs_etrons, sigma_total_etrons, NAME='$\sigma_{TOTAL}(e^{-})$', THICK=2, $ 
           LINESTYLE=4, /OVERPLOT)
 
 ; Add the LEGEND to the plot

 leg_etrons = LEGEND(TARGET=[p22,p33,p44,p55], POSITION=[40,4000], /DATA)
           
END  

;****************************************************

PRO ex_5_2

  COMPILE_OPT IDL2, hidden

  ; Using the information provided and the data file in 
  ; this directory: ex_5.2.txt generate PTC for sigma_total, 
  ; sigma_read+shot, and sigma_fpn in DN. Determine K_ADC, 
  ; P_n, sigma_read, and S_fw.

  ; Read in the ascii data file
  ; Initialize files and a template
  
  txtFile = 'ex_5.2.txt'
  
  txt_dat = 'ex_5.2.txt.dat'
  
  ; Restore the ASCII file template

  templ_file = 'ex_5_2.templ.sav'

  IF FILE_TEST(templ_file) THEN BEGIN  
  
    restore, filename=templ_file, /VERBOSE
    print,' '
    print,'********************'
    print,' '
  ;  WAIT, 5
  
  ENDIF ELSE BEGIN
    
    templ_file_err = DIALOG_MESSAGE('I did not find the file: '+templ_file, $
      TITLE='File not Found', /CENTER, /ERROR )
      
    GOTO, DONE
    
  ENDELSE
  
  ; Test if the data ascii file is present. If it is here, 
  ; read it in using the template. If the text file is not
  ; present an error dialog box will appear.
  
  IF FILE_TEST(txtFile) THEN BEGIN
    
      IF ~FILE_TEST(txt_dat) THEN BEGIN

      txt_dat = READ_ASCII(txtFile, TEMPLATE=the_template, /VERBOSE)
        
      ENDIF
      
  ENDIF ELSE BEGIN
      
      txtFile_err = DIALOG_MESSAGE('I did not find the file: '+txtFile, $
                    TITLE='File not Found', /CENTER, /ERROR )
                    
      GOTO, DONE
    
  ENDELSE

  ; Set up a blank log-log plot, then overplot the data

  p1 = PLOT(txt_dat.raw, txt_dat.total_noise, XLOG=1, YLOG=1, XRANGE=[1,0.5e6], YRANGE=[0.01,1e4], $
    XTITLE='Signal (DN)', YTITLE='NOISE (DN)', FONT_SIZE=12, $
    FONT_STYLE='Bold', FONT_NAME='Times', Title = 'Example 5.2 DN $\rightarrow$ $\lambda$', $
    XTHICK=2, YTHICK=2, /NODATA)
    
  ; Plot the signal and sigma_signal from the txt_dat structure
  
  p2 = PLOT(txt_dat.signal, txt_dat.total_noise, "b4D", NAME='$\sigma_{TOT}(DN)$', /OVERPLOT)
  
  ; Plot the signal and sigma_fpn and determine the P_n
  
  p3 = PLOT(txt_dat.signal, txt_dat.fpn_sigma, "b4+", NAME='$\sigma_{FPN}(DN)$', /OVERPLOT)
   
  ; Query user on number of points to exclude from P_n averaging
  
  num_exclude = 0
  
  READ, num_exclude, PROMPT='Enter number of FPN points to exclude at or after full well: '
  
  ; Get an average P_n and sigma_P_N
  
  many_P_n = txt_dat.fpn_sigma/txt_dat.SIGNAL
  
  avg_P_N = MOMENT(many_P_n[0:n_elements(many_P_N)-(1 + num_exclude)])
  
  ; P_n from one data pair:
  
  one_P_n = txt_dat.fpn_sigma[17]/txt_dat.SIGNAL[17]
  
  print,' '
  print,' Average P_n = ',STRTRIM(STRING(avg_P_n[0]), 1),' +/- ', $
          STRTRIM(STRING(avg_P_n[1]), 1)
  print,' '
  print,' P_n from one FPN data pair = ',STRTRIM(STRING(one_P_n), 1)
  print,' '

  ; Plot the shot noise PTC
  
  p4 = PLOT(txt_dat.signal, txt_dat.shot_noise, "r4*", NAME='$\sigma_{shot}(DN)$', /OVERPLOT)
  
  ; Report the K_adc using element #17
  
  K_adc = txt_dat.SIGNAL[17]/txt_dat.SHOT_NOISE[17]^2
  
  print,' '
  print,' K_adc from one Shot Noise data pair = ',STRTRIM(STRING(K_adc), 1)
  print,' '
  
  ; Calculate K_adc from a linear fit to the FPN line. Knowing slope and intercept
  ; get a y value from an x value and K_adc = x/y^2. Note that LINFIT 
  ; will be passed a keyword to return a vector of calculated y-fit values.
  
  y_to_fit = txt_dat.shot_noise[15:n_elements(txt_dat.fpn_sigma)-(1 + num_exclude)]
  x_to_fit = txt_dat.SIGNAL[15:n_elements(txt_dat.fpn_sigma)-(1 + num_exclude)]
  
  ; stop
  
  print,' '
  print,' Fitting some of the Shot noise data for K_adc ...'
  print,' '
  
  FPN_fit = LINFIT(x_to_fit, y_to_fit, CHISQR=chisq, COVAR=cvar, PROB=probal, SIGMA=linfit_sigma, $ 
                   YFIT=y_fitted)
  
    print,' Shot noise slope        = ',STRTRIM(STRING(FPN_fit[1]), 1)
    
    ; K_adc is found when the y value for noise = 1 DN
    
  ; print,' Shot noise  y-intercept = ',STRTRIM(STRING(FPN_fit[0]), 1)
  
  ; Get the average of the fitted y values to the signal values for a K_adc
  ; estimate
  
  mean_K_adc = MEAN(x_to_fit/(y_fitted^2))
  ; print,' '
  print,' Average K_adc from the fit is: = ',STRTRIM(STRING(mean_K_adc), 1)
  print,' '
  
  ;stop

  ; Calculate sigma_read usign Equations 5.10 and 5.11 of Janesick's book
  ; excluding first the saturated points from sigma_sqrd_shot_read
  
  sigma_sqrd_shot_read = txt_dat.total_noise[0:n_elements(txt_dat.total_noise) - num_exclude]^2 - $ 
                         txt_dat.fpn_sigma[0:n_elements(txt_dat.total_noise) - num_exclude]^2
  
  ; Next, make a check in the existing arrays within the SQRT to that I do not 
  ; NaN if the difference is < 0
  
  good_ones = WHERE(sigma_sqrd_shot_read GE txt_dat.shot_noise^2)
  
  sigma_read = SQRT(sigma_sqrd_shot_read[good_ones] - txt_dat.shot_noise[good_ones]^2)
  
  sigma_read_mom = MOMENT(sigma_read)
  
  print,' '
  print,' Average sigma_read = ',STRTRIM(STRING(sigma_read_mom[0]), 1),' +/- ', $
    STRTRIM(STRING(SQRT(sigma_read_mom[1])), 1)
  print,' '
  
  
  ; Graphically sigma_read looks to be about 5 DN
  
  rd_n_xs = [1,3.5e05]
  rd_n_ys = [5.0,5.0]
  
  p5 = PLOT(rd_n_xs, rd_n_ys, "k:", NAME='$\sigma_{RN}(DN)$', THICK=2, /OVERPLOT)

  ; Plot a vertical line for the full well
  
  fw_xs = [41800,41800]
  fw_ys = [0.01,10000]
  
  p6 = PLOT(fw_xs, fw_ys, "k-.", NAME='FW', THICK=2, /OVERPLOT)
  
  
  ; Add the LEGEND to the plot

  leg = LEGEND(TARGET=[p2,p3,p4,p5, p6], POSITION=[40,4000], /DATA)

  ; I get to the DONE: if I cannot find the template file or the ASCII 
  ; data file
  DONE:
  
  END

;****************************************************

PRO ex_5_3

  COMPILE_OPT IDL2, hidden

  ; Generate Photon Transfer Curves (PTC's) for sigma_total,
  ; sigma_read, sigma_shot and sigma_FPN given the parameters
  ; listed in Example 5.1 of Janesick's This time with eta's Quantum
  ; Yields greater than one

  ; "Photon Transfer DN -> Lambda".

  ; Slopes are for a log-log plot.

  ; sigma_read = constant = 3.33 DN

  ; maximum signal at saturation (ordinate) = 233,000 = 10^5.3674

  ; maximum sigma_FPN (slope = 1) = P_n*(Full well) =
  ;                               = 0.02*233,000 = 4,666 = 10^3.669

  ; maximum sigma_shot (slope = 1/2) = SQRT(233,000) = 482.7
  ;                                  = 10^2.6837

  ; Set up a blank log-log plot, then overplot the data

  xs = DINDGEN(400000, START = 0.01d)
  ; let ys all = 1.0
  ys = xs/xs

  ; Array for the Quantum Yields eta's

  etas = [1.0, 3.0, 10.0, 25.0, 50.0, 100.0]

  ; Make a sort-of 'Identity' array to expand out the read noise and
  ; FPN noise column vectors into 6 by 400000 arrays to match the dimensions
  ; of the shot noise column vector expanded by the etas

  ones_array = FLTARR(n_elements(etas)) + 1.0d

  p1 = PLOT(xs, ys, XLOG=1, YLOG=1, XRANGE=[1,1e6], YRANGE=[1,1e4], $
    XTITLE='Signal (DN)', YTITLE='NOISE (DN)', FONT_SIZE=14, $
    FONT_STYLE='Bold', FONT_NAME='Times', $ 
    Title = 'Example 5.3 DN $\rightarrow$ $\lambda \n \eta_{i} \ge 1$', $
    XTHICK=0.5, YTHICK=0.5, XTICKLEN=1.0, YTICKLEN=1.0, /NODATA)

  ; Plot the read noise

  sigma_read = 3.33*ys         ;

  ; p2 = PLOT(xs, sigma_read, NAME='$\sigma_{RN}(DN)$', THICK=2, /OVERPLOT)

  ; Plot the shot noise
  ; Minimum signal is 1 electron generates 1.5 DN

  K_adc = 1.5 ; for eta = 1.0

  sigma_shot = SQRT(xs/K_adc)  ; Eqn. 5.7 = Eqn. 4.20

  ; Modify for full well

  ; fw = WHERE(xs GE 233000)
  fw = WHERE(xs GE 220000)

  ; Use the modified Wien's law for the saturation tail
  sigma_shot[fw] = my_wien(xs[fw],sigma_shot[fw],233000)

  P_n = 0.02

  sigma_fpn = P_n*xs            ; Eqn. 5.8 = Eqn. 3.12

  ; Allow for full well FPN
  sigma_fpn[fw] = my_wien(xs[fw],sigma_fpn[fw],233000)

  ; sigma_total is the quadrature sum of Read Noise, Shot Noise and FPN Noise

  ; Eqn. 5.9
  sigma_total = SQRT(ones_array#sigma_read^2 + etas#sigma_shot^2 + ones_array#sigma_fpn^2)

  ; stop

  p2 = PLOT(xs, sigma_total[0,*], NAME='$\eta = 1.0$',   THICK=2, /OVERPLOT)
  p3 = PLOT(xs, sigma_total[1,*], NAME='$\eta = 3.0$',   THICK=2, /OVERPLOT)
  p4 = PLOT(xs, sigma_total[2,*], NAME='$\eta = 10.0$',  THICK=2, /OVERPLOT)
  p5 = PLOT(xs, sigma_total[3,*], NAME='$\eta = 25.0$',  THICK=2, /OVERPLOT)
  p6 = PLOT(xs, sigma_total[4,*], NAME='$\eta = 50.0$',  THICK=2, /OVERPLOT)
  p7 = PLOT(xs, sigma_total[5,*], NAME='$\eta = 100.0$', THICK=2, /OVERPLOT)
  
  ; Plot the shot noise only
  
  eta_sigma_shot = SQRT(etas)#sigma_shot
  
  p22 = PLOT(xs, eta_sigma_shot[0,*], NAME='$\eta_{SHOT} = 1.0$', $ 
            LINESTYLE=2, THICK=2, /OVERPLOT)
  p32 = PLOT(xs, eta_sigma_shot[1,*], NAME='$\eta_{SHOT} = 3.0$', $ 
            LINESTYLE=2, THICK=2, /OVERPLOT)
  p42 = PLOT(xs, eta_sigma_shot[2,*], NAME='$\eta_{SHOT} = 10.0$', $
            LINESTYLE=2, THICK=2, /OVERPLOT)
  p52 = PLOT(xs, eta_sigma_shot[3,*], NAME='$\eta_{SHOT} = 25.0$', $
            LINESTYLE=2, THICK=2, /OVERPLOT)
  p62 = PLOT(xs, eta_sigma_shot[4,*], NAME='$\eta_{SHOT} = 50.0$', $
            LINESTYLE=2, THICK=2, /OVERPLOT)
  p72 = PLOT(xs, eta_sigma_shot[5,*], NAME='$\eta_{SHOT} = 100.0$', $ 
            LINESTYLE=2, THICK=2, /OVERPLOT)
  
  ; Determine Signal levels in electrons when sigma_shot = sigma_FPN
  ; using Eqn. 3.14
  
  noise_equiv=etas/(P_n^2)
  
  print,' '
  print,' Shot noise = FPN noise at a signal level of:'
  print,' For eta = ',STRTRIM(STRING(etas[0]),1),' a Signal of ', $
          STRTRIM(STRING(noise_equiv[0]),1),' electrons.'
  print,' For eta = ',STRTRIM(STRING(etas[1]),1),' a Signal of ', $
          STRTRIM(STRING(noise_equiv[1]),1),' electrons.'
  print,' For eta = ',STRTRIM(STRING(etas[2]),1),' a Signal of ', $
          STRTRIM(STRING(noise_equiv[2]),1),' electrons.'
  print,' For eta = ',STRTRIM(STRING(etas[3]),1),' a Signal of ', $
          STRTRIM(STRING(noise_equiv[3]),1),' electrons.'
  print,' For eta = ',STRTRIM(STRING(etas[4]),1),' a Signal of ', $
          STRTRIM(STRING(noise_equiv[4]),1),' electrons.'
  print,' For eta = ',STRTRIM(STRING(etas[5]),1),' a Signal of ', $
          STRTRIM(STRING(noise_equiv[5]),1),' electrons.'
  print,' '
  
  ; Make a re-scaled plot so one can determine the individual
  ; quantum yields eta_i using Eqn. 4.5
  
  ; stop
  
  p12 = PLOT(xs, ys, XLOG=1, YLOG=1, XRANGE=[0.01,100], YRANGE=[0.1,100.0], $
    XTITLE='Signal (DN)', YTITLE='NOISE (DN)', FONT_SIZE=14, $
    FONT_STYLE='Bold', FONT_NAME='Times', $ 
    Title = 'Example 5.3 DN $\rightarrow$ $\lambda \n \eta_{i} \ge 1$', $
    XTHICK=0.5, YTHICK=0.5, XTICKLEN=1.0, YTICKLEN=1.0, /NODATA)

  ; Now re-execute the plotting commands for the sigma_total data for 
  ; quantum yields >= 1.0

  p23 = PLOT(xs, sigma_total[0,*], NAME='$\eta = 1.0$',   THICK=2, /OVERPLOT)
  p33 = PLOT(xs, sigma_total[1,*], NAME='$\eta = 3.0$',   THICK=2, /OVERPLOT)
  p43 = PLOT(xs, sigma_total[2,*], NAME='$\eta = 10.0$',  THICK=2, /OVERPLOT)
  p53 = PLOT(xs, sigma_total[3,*], NAME='$\eta = 25.0$',  THICK=2, /OVERPLOT)
  p63 = PLOT(xs, sigma_total[4,*], NAME='$\eta = 50.0$',  THICK=2, /OVERPLOT)
  p73 = PLOT(xs, sigma_total[5,*], NAME='$\eta = 100.0$', THICK=2, /OVERPLOT)
  
  ; Now plot the shot noise only lines so that the eta_i may be determined 
  ; graphically for the different values of quantum yield. 

  p24 = PLOT(xs, eta_sigma_shot[0,*], NAME='$\eta_{SHOT} = 1.0$', $
    LINESTYLE=2, THICK=2, /OVERPLOT)
  p34 = PLOT(xs, eta_sigma_shot[1,*], NAME='$\eta_{SHOT} = 3.0$', $
    LINESTYLE=2, THICK=2, /OVERPLOT)
  p44 = PLOT(xs, eta_sigma_shot[2,*], NAME='$\eta_{SHOT} = 10.0$', $
    LINESTYLE=2, THICK=2, /OVERPLOT)
  p54 = PLOT(xs, eta_sigma_shot[3,*], NAME='$\eta_{SHOT} = 25.0$', $
    LINESTYLE=2, THICK=2, /OVERPLOT)
  p64 = PLOT(xs, eta_sigma_shot[4,*], NAME='$\eta_{SHOT} = 50.0$', $
    LINESTYLE=2, THICK=2, /OVERPLOT)
  p74 = PLOT(xs, eta_sigma_shot[5,*], NAME='$\eta_{SHOT} = 100.0$', $
    LINESTYLE=2, THICK=2, /OVERPLOT)

  
  ; Add LEGENDs to the plot
  
  labels1 = ['0.015','0.03','0.05','0.15','0.5','1.5 = $K_{ADC}(e^{-}/DN)$']
  xloc_lbl = [0.03,0.04,0.1,0.3,0.8,1.6]
  yloc_lbl = [1.6, 1.2, 1.3, 1.4,1.3,1.1]

  text1 = TEXT(xloc_lbl, yloc_lbl, labels1, /DATA, FILL_BACKGROUND=1, $
              COLOR='red', FONT_STYLE='bold', FONT_SIZE=10)
  text2 = TEXT(0.015, 2.2, '$K_{ADC}(P_{I}/DN)$',/DATA,FILL_BACKGROUND=1, $
              COLOR='red', FONT_STYLE='bold', FONT_SIZE=10)

  
END

;****************************************************

PRO ex_5_4

  COMPILE_OPT IDL2, hidden

  ; "Photon Transfer DN -> Lambda".

  ; Using a measured relation given in the example, 
  ; generate a PTC for the data and plot a PTC and an 
  ; image lag factor I_lag as a function of signal.

  ; Set up a blank log-log plot, then overplot the data

  xs = DINDGEN(40000, START = 0.01d)
  ys = xs
  
  ; Set the ADC Sensitivity K_adc

  K_adc = 19.0 ; electrons/DN

  ; Initiate the plot window

  p1 = PLOT(xs, ys, XLOG=1, YLOG=1, XRANGE=[100,1e4], YRANGE=[10,1e3], $
    XTITLE='Signal (DN)', YTITLE='NOISE (DN)', FONT_SIZE=14, $
    FONT_STYLE='Bold', FONT_NAME='Times', $
    Title = 'Example 5.4 DN $\rightarrow$ $\lambda \n Image Lag (I_{LAG}) Part 1$', $
    XTHICK=0.5, YTHICK=0.5, XTICKLEN=1.0, YTICKLEN=1.0, /NODATA)
  
  ; Create the vectors for the Ideal shot noise sigma_shot  
  ; measured shot noise sigma_m_shot

  ; A normalization factor is included to match the observed 
  ; Fig. 5.19. This is SQRT(K_adc). I call it out separately here 
  ; to remind myself of Equation 4.20
  
  sigma_shot = SQRT(K_adc)*SQRT(xs)

  sigma_m_shot = sigma_shot*(1.0 - 0.29*exp(-1.0*(xs/1579.0)))

  ; Calculate the Image lag
  
  imLag = sigma_m_shot/sigma_shot

  ; Plot the ideal shot noise sigma_shot
  
  p2 = PLOT(xs, sigma_shot, '-', THICK=2, /OVERPLOT, NAME='Ideal $\sigma_{SHOT}$')
  
  ; Plot the measured shot noise sigma_m_shot
  
  p3 = PLOT(xs, sigma_m_shot, '--', THICK=2, /OVERPLOT, NAME='$\sigma_{M\_SHOT}$')
  
  leg = LEGEND(POSITION=[500, 500], TARGET=[p2, p3], /DATA, FONT_NAME='Times', $ 
              FONT_STYLE='bold')
              
  pl_imLag = PLOT(xs, imLag, XLOG=1, XRANGE=[100,1e4], YRANGE=[0.65,1.05], $
                 XTITLE='Signal (DN)', YTITLE='Image lag ($\sigma_{M\_SHOT}/\sigma_{SHOT}$)', $
                 Title = 'Example 5.4 DN $\rightarrow$ $\lambda \n Image Lag (I_{LAG}) Part 2$', $
                 XTHICK=0.5, YTHICK=0.5, XTICKLEN=1.0, YTICKLEN=1.0, THICK=2.0, FONT_SIZE=14, $ 
                 FONT_NAME='Times')

END

;****************************************************

PRO ex_5_5

  COMPILE_OPT IDL2, hidden

  ; Use the Variance PTC method with the numbers from Example 5_1
  ; to plot noise variance as a function of signal on a 
  ; linear scale. Determine the K_adc and the read noise
  ; from the graph
  
  ; Plot the total variance versus signal for comparison to 
  ; the shot + read_noise variance
  
  ; Set up a blank plot, then overplot the data

  xs = DINDGEN(300000, START = 0.0d)
  ; let ys all = 1.0
  ys = xs/xs

  p1 = PLOT(xs, ys, XRANGE=[1,3e5], YRANGE=[1,1.6e5], $
    XTITLE='Variance (DN)', YTITLE='NOISE (DN)', FONT_SIZE=12, $
    FONT_STYLE='Bold', FONT_NAME='Times', $ 
    Title = 'Example 5.5 DN $\rightarrow$ $\lambda \n$ Variance PTC', $
    XTHICK=2, YTHICK=2, /NODATA)

  ; Set the read noise

  sigma_read = 3.33*ys
  sigma_read[0] = 3.33
  
  var_read = sigma_read^2
  
  K_adc = 1.5 
  
  var_shot = xs/K_adc
  
  P_n = 0.02

  var_fpn = (P_n*xs)^2
  
  tot_var = var_read + var_shot +var_fpn
  
  var_read_shot = var_read + var_shot
  
  ; Allow for saturation with a quadratic as shown in handwritten 
  ; notes for this example
  
  ; fw = WHERE(xs GE 233000)
  fw = WHERE(xs GE 200000)
  
  var_shot_read_sat = (-1.0e-5)*(xs[fw] - 233000.0d)^2 + 143849.0d
  
  var_read_shot[fw] = var_shot_read_sat
  
  ; Plot the two variance curves
  
  p2 = PLOT(xs, var_read_shot, '-', THICK=2, NAME='$\sigma^{2}_{READ+SHOT}(DN)$', $
           FONT_NAME='Times', /OVERPLOT)

  p3 = PLOT(xs, tot_var,'-:', THICK=2, NAME='$\sigma^{2}_{TOTAL}(DN)$', $
         FONT_NAME='Times', /OVERPLOT)
         
  leg = LEGEND(POSITION=[1.1e5, 1.3e5], TARGET=[p2, p3], /DATA, FONT_NAME='Times', $
         FONT_STYLE='bold')

  ; Indicate the full well

  x_fw = [233000.0d, 233000.0d]
  y_fw = [0.0d, 160000.0d]
  
  p4 = PLOT(x_fw, y_fw, '__', THICK=2, /OVERPLOT)
  
  fw_text = TEXT(1.8e5, 1.1e05, 'Full Well $\rightarrow$',/DATA,FILL_BACKGROUND=1, $
    COLOR='red', FONT_STYLE='bold', FONT_SIZE=10)


  ; Zoom in on low signal regime for read noise
  
  p12 = PLOT(xs, ys, XRANGE=[0,100], YRANGE=[0,100], $
    XTITLE='Variance (DN)', YTITLE='NOISE (DN)', FONT_SIZE=12, $
    FONT_STYLE='Bold', FONT_NAME='Times', $
    Title = 'Example 5.5 DN $\rightarrow$ $\lambda \n$ Variance PTC', $
    XTHICK=2, YTHICK=2, /NODATA)
  
  p22 = PLOT(xs, var_read_shot, '-', THICK=2, NAME='$\sigma^{2}_{READ+SHOT}(DN)$', $
    FONT_NAME='Times', /OVERPLOT)

  p32 = PLOT(xs, tot_var,'-:', THICK=2, NAME='$\sigma^{2}_{TOTAL}(DN)$', $
    FONT_NAME='Times', /OVERPLOT)
    
  p42 = PLOT(xs[0:30], var_read[0:30], 'b-.', THICK=2, NAME='$\sigma^{2}_{READ}(DN)$', $
    FONT_NAME='Times', /OVERPLOT)

  p5 = PLOT(xs, var_shot, 'r__', THICK=2, NAME='$\sigma^{2}_{SHOT}(DN)$', $  
    FONT_NAME='Times', /OVERPLOT)

  ;stop

  ; Insert a legend
  
  leg = LEGEND(POSITION=[35, 90], TARGET=[p22, p32, p42, p5], /DATA, FONT_NAME='Times', $
    FONT_STYLE='bold')

END

;****************************************************

PRO ex_5_6

  COMPILE_OPT IDL2, hidden

  ; Use a Variance PTC and non-zero Quantum yields to 
  ; demonstrate a calculation of K_adc for the case of 
  ; eta = 1 and 10. This is a simple example that doesn't mention 
  ; read noise so I will set that to 3.33 DNas usual
  
  ; Set up a blank plot, then overplot the data

  xs = DINDGEN(300000, START = 0.0d)
  ; let ys all = 1.0
  ys = xs/xs

  p1 = PLOT(xs, ys, XRANGE=[0,4e3], YRANGE=[0,500], $
    XTITLE='Variance (DN)', YTITLE='NOISE (DN)', FONT_SIZE=12, $
    FONT_STYLE='Bold', FONT_NAME='Times', $
    Title = 'Example 5.6 DN $\rightarrow$ $\lambda \n$ Variance PTC $\eta \ge 1.0$', $
    XTHICK=2, YTHICK=2, /NODATA)

  ; Set the read noise

  sigma_read = 3.33*ys
  sigma_read[0] = 3.33

  ; Make a small array to hold the two etas
  
  etas = [1.0,10.0]

  ; Make a sort-of 'Identity' array to expand out the read noise 
  ; vectors into arrays to match the dimensions
  ; of the shot noise column vector expanded by the etas

  ones_array = FLTARR(n_elements(etas)) + 1.0d

  K_adc = 100.0d ; electrons/DN given in problem
  
  ; Create the array of variances for shot noise
  
  var_shot = etas#(xs/K_adc)
  
  var_total = etas#sigma_read + var_shot
  
  ; Run the plots

  ;stop

  pl2 = PLOT(xs, var_total[0,*], '__', THICK=2, Name='$\eta_{i} = 1$', /OVERPLOT)
  
  pl3 = PLOT(xs, var_total[1,*], '-', THICK=2, Name='$\eta_{i} = 10$', /OVERPLOT)
  
  ; Add a legend

  leg = LEGEND(POSITION=[1200, 400], TARGET=[pl2, pl3], /DATA, FONT_NAME='Times', $
    FONT_STYLE='bold')
    
  ; Get the slopes, get the eta
    
  steep_slope = (var_total[1,2000] - var_total[1,1000])/(xs[2000] - xs[1000])
  
  shallow_slope = (var_total[0,2000] - var_total[0,1000])/(xs[2000] - xs[1000])  
  
  big_eta = steep_slope/shallow_slope
  
  ; Instead of STRTRIM(STRING(big_eta,1) use an object
  
  ; convert big_eta to a string
  
  a = STRING(big_eta)
  
  out_string_cl = ' Quantum efficiency eta = '+STRMID(a.Trim(),0,4)+' electrons per photon.'
  
  print,' '
  print,out_string_cl
  print,' '
  
  out_string = '$\leftarrow\ \eta\ = $'+STRMID(a.Trim(),0,4)+' e$^{-}/\gamma$'
  
  eta_text = TEXT(2000, 190, out_string,/DATA, $ 
                 FILL_BACKGROUND=1, COLOR='red', FONT_STYLE='bold', FONT_SIZE=12)
                 
  stop
    
END

;**********

;**********
; Chapter 7
;**********

;****************************************************

PRO ex_7_1

  COMPILE_OPT IDL2

  ; Make plots versus time (t: sec) of S (sense node electrons),
  ; sigma_shot(Sense node e-trons), sigma_fpn(Sense node e-trons),
  ; S(DN), sigma_shot(DN), and sigma_fpn(DN). Generate PTC for 
  ; for sigma_shot(DN), and sigma_fpn(DN). Graphically validate 
  ; P_n and determine extreme values for K_adc. Plot K_adc as a 
  ; function of S (sense node e-trons) and S(DN). Assume P_n = 0.03. 
  
  ; make an array for times in seconds
  ts = DINDGEN(10001, START = 0.0d)

  ; S = t (given)
  S = ts
  
  ; sigma_shot assumes Poisson noise
  sigma_shot = SQRT(S)
  
  ; sigma_fpn is linear and multiplied by P_n
  ; and graphically P_n is in the same location on the Dn and e-tron plots 
  P_n = 0.03d
  sigma_fpn = P_n*S
  
  ; Make the plots of the sense-node electrons versus time
  
  pl_S = PLOT(ts, S, '-', NAME='S(e$^{-}_{SN}$)', XRANGE=[1,41000], YRANGE=[1,41000], $
             XLOG=1, YLOG=1, XTITLE='Time (sec)', YTITLE='Sense node $e^{-}/$DN', $
             FONT_SIZE=14, FONT_STYLE='Bold', FONT_NAME='Times', $
             TITLE='Example 7-1 part A $\n$ Sense Node e$^{-}$', XTHICK=0.5, YTHICK=0.5, $ 
             THICK=2)

  pl_S2 = PLOT(ts, sigma_shot, '-:', NAME='$\sigma_{SHOT}$', THICK=2, /OVERPLOT)
  
  pl_S3 = PLOT(ts, sigma_fpn, '--', NAME='$\sigma_{FPN}$', THICK=2, /OVERPLOT)

  leg = LEGEND(POSITION=[15, 1500], TARGET=[pl_S, pl_S2, pl_S3], /DATA, $ 
              FONT_NAME='Times', FONT_STYLE='bold', FONT_SIZE=12)
  
  p_n_string = 'P$_{n}$ = 1/33'
  
  p_n_text = TEXT(120, 1.5, p_n_string, /DATA, FONT_SIZE=12, FILL_COLOR='white', $ 
                 FILL_BACKGROUND=1, COLOR='blue', FONT_STYLE='bold')
                 
  ; Now plot in DN to se the non-linear effects of the modeled Source Follower 
  ; (SF) amplifier

  K_adc = 1.5^(1.0 + 0.001*ts)

  S_DN = ts/(K_adc)
  
  sigma_shot_dn = SQRT(S_DN/K_adc) ; Eqn. 5.7
  
  sigma_fpn_dn = P_n*S_DN
  
  pl_S_DN = PLOT(ts, S_DN, '-', NAME='S($_{DN}$)', XRANGE=[1,41000], YRANGE=[1,1000], $
    XLOG=1, YLOG=1, XTITLE='Time (sec)', YTITLE='DN', XTICKLEN=0, YTICKLEN=0, $
    FONT_SIZE=14, FONT_STYLE='Bold', FONT_NAME='Times', $
    TITLE='Example 7-1 part B $\n$ Source Follower (SF) non-linearity', XTHICK=0.5, $
    YTHICK=0.5, THICK=2)
    
  pl_S_DN_2 = PLOT(ts, sigma_shot_dn, '-:', NAME='$\sigma_{SHOT}(DN)$', THICK=2, $ 
                  /OVERPLOT)
                  
  pl_S_DN_3 = PLOT(ts, sigma_fpn_dn, '--', NAME='$\sigma_{FPN}(DN)$', THICK=2, $
                  /OVERPLOT) 
  
  leg = LEGEND(POSITION=[25, 150], TARGET=[pl_S_DN, pl_S_DN_2, pl_S_DN_3], /DATA, $
                  FONT_NAME='Times', FONT_STYLE='bold', FONT_SIZE=12)
  
  ; Plots versus DN
  
  pl_DN = PLOT(S_DN, sigma_shot_dn, 'r-', NAME='$\sigma_{SHOT}(DN)$', THICK=2, $ 
              XRANGE=[1,1000], YRANGE=[0.1,50], XLOG=1, YLOG=1, XTITLE='S(DN)', $ 
              YTITLE='$\sigma$(DN)', FONT_SIZE=14, FONT_STYLE='Bold', $ 
              XTICKLEN=1, YTICKLEN=1, $
              FONT_NAME='Times', TITLE='Example 7-1 part C $\n$ Determine $K_{ADC}(e^{-}/DN)$')
              
  pl_DN_2 = PLOT(S_DN, sigma_fpn_dn, 'r-:', NAME='$\sigma_{FPN}(DN)$', THICK=2, $ 
                /OVERPLOT) 
                
  ; Seeing the sigma_shot curve the min and max values of K_adc from the intercepts
  ; with x-intercepts = 1 are 1.5 and about 90 as the slope 1/2 with 
  ; K_adc = 90 meets the end of the sigma_shot_dn curve at ~[103,1.2]
  
  min_k_sig_dn = SQRT(S_DN/1.5d)
  max_k_sig_dn = SQRT(S_DN/90.0d)
  
  pl_low_k = PLOT(S_DN, min_k_sig_dn, 'b:', NAME='K$_{ADC}$(low) = 1.5', $ 
                 THICK=2, /OVERPLOT)
                 
  pl_high_k = PLOT(S_DN, max_k_sig_dn, 'b-.', NAME='K$_{ADC}$(high) = 90.0', $
                  THICK=2, /OVERPLOT)
   
  leg = LEGEND(POSITION=[20,40], TARGET=[pl_DN, pl_DN_2, pl_low_k, pl_high_k], $ 
              /DATA, FONT_NAME='Times', FONT_STYLE='bold', FONT_SIZE=12)
                
  p_n_string = 'P$_{n}$ = 1/33'

  p_n_text = TEXT(10, 1.1, p_n_string, /DATA, FONT_SIZE=12, FILL_COLOR='white', $
                 FILL_BACKGROUND=1, COLOR='blue', FONT_STYLE='bold')
               
  ; Calculate and plot K_adc (e-tron/DN) versus DN and electrons
  
  K_adc_DN = S_DN/(sigma_shot_dn)^2
 
  pl_K_dn = PLOT(S_DN, K_adc_dn, 'r-', THICK=2, $ 
               XRANGE=[1,1000], YRANGE=[1,500], XLOG=1, YLOG=1, XTITLE='S(DN)', $ 
               YTITLE='K$_{ADC}(e^{-}/DN)$', FONT_SIZE=14, FONT_STYLE='Bold', $ 
               FONT_NAME='Times', $ 
               TITLE='Example 7-1 part D $\n$ V/V non-linearity [SF amplifier]')
               
  ; Knowing the sensitivity we can work back with photon transfer and plot 
  ; K_adc_dn versus electrons from the sense node
  
  non_lin_e = S_DN*K_adc_DN
  
  pl_K_e =  PLOT(non_lin_e, K_adc_dn, 'r-', THICK=2, $
            XRANGE=[1,10000], YRANGE=[1,500], XLOG=1, YLOG=1, XTITLE='S(e$^{-}$)', $
            YTITLE='K$_{ADC}(e^{-}/DN)$', FONT_SIZE=14, FONT_STYLE='Bold', $
            FONT_NAME='Times', $
            TITLE='Example 7-1 part E $\n$ V/V non-linearity [SN electrons]')
            
END

;****************************************************

PRO ex_7_2

  COMPILE_OPT IDL2

  ; Sense node V/e-tron non-linearity.  

  ; Reference voltage at Drain of reset FET transistor. 
  V_ref = 3.1
  
  ; Capacitance constant
  k1 = 10.909e-15 ; femtoFarad = 10^-15 Farad
  ; Charge of electron constant
  q = 1.6e-19 ; Coulombs

  ; Array for the Sense node voltage range
  V_sn = DINDGEN(222, INCREMENT=0.01, START=0.90)
  ; Sense node capacitance
  C_sn = k1/V_sn
  
  ; Arrays for vertical lines for Voltage extremes on Cap
  Vh   = [3.1, 3.1]
  Vh_C = [0.0, 12.0]
  
  Vl   = [0.9,0.9]
  Vl_C = Vh_C
  
  pl_C_sn = PLOT(V_sn, C_sn/1e-15, XRANGE=[0.7,3.3], YRANGE=[2.0,13.0], $
                XTITLE='Sense node Voltage (V$_{SN}$, Volts)', $
                YTITLE='C$_{SN} (femtoFarads)$', $ 
                Title='Exercise 7.2 V/V non-linearity $\n$ Sense Node Capacitance', $
                THICK=2, FONT_NAME='Times', FONT_STYLE='Bold', FONT_SIZE=14)
  pl_h = PLOT(Vh, Vh_C, '--', THICK=2, /OVERPLOT)
  pl_l = PLOT(Vl, Vl_C, '--', THICK=2, /OVERPLOT)
  ; Add text to indicate direction of curve for charging
  
  cap_string = '$\leftarrow$ Capacitor charging'
  reset = 'Reset'
  disch = 'Discharged'

  cap_text = TEXT(1.7, 8.0, cap_string, FONT_NAME='TIMES', FONT_STYLE='Bold', $ 
                 FONT_SIZE=14, /DATA)

  diode_reset = TEXT(3.0,12.1, reset, FONT_NAME='TIMES', FONT_STYLE='Bold', $
                FONT_SIZE=12, /DATA)
  diode_disch = TEXT(0.8,12.2, disch, FONT_NAME='TIMES', FONT_STYLE='Bold', $
                FONT_SIZE=12, /DATA)

  ; Plot S(electrons) and sigma_shot(electrons) versus V_SN
  ; Use the solution of the differential equation for S(e^-) 
  ; and sigma_shot(e^-)
   
  S_etron    = (k1/q)*ALOG(V_ref/V_sn)
  sigma_shot = SQRT(S_etron)
  
  pl_e = PLOT(V_sn, S_etron, XRANGE=[0.7,3.3], YLOG=1, YRANGE=[7,3e5], $
    XTITLE='Sense node Voltage (V$_{SN}$, Volts)', $
    YTITLE='Sense node Electrons', $
    Title='Exercise 7.2 V/V non-linearity $\n$ Sense Node Signal and noise', $
    THICK=2, FONT_NAME='Times', FONT_STYLE='Bold', FONT_SIZE=14, NAME='S$_{SN}$')
    
  pl_sig = PLOT(V_sn, sigma_shot, '-:', THICK=2, NAME='$\sigma_{SHOT}$', /OVERPLOT)
              
  leg_e = LEGEND(POSITION=[1.7, 10000], TARGET=[pl_e, pl_sig], $ 
                /DATA, FONT_NAME='Times', FONT_STYLE='bold', FONT_SIZE=12) 

  diode_reset = TEXT(3.0,105000, reset, FONT_NAME='TIMES', FONT_STYLE='Bold', $
                FONT_SIZE=12, /DATA)
  diode_disch = TEXT(0.8,105000, disch, FONT_NAME='TIMES', FONT_STYLE='Bold', $
                FONT_SIZE=12, /DATA)

              pl_h = PLOT(Vh, [5,100000], '--', THICK=2, /OVERPLOT)
              pl_l = PLOT(Vl, [5,100000], '--', THICK=2, /OVERPLOT)

  ; Part C Plot S(e/V_sn) and N(e/SN) versus V_sn
  
  ; S_{SN}(e/V_sn) = S(e)/S(V_sn) AND S(V_sn) = V_ref - V_sn
  
  S_V_sn = V_ref -1.0*V_sn
  
  S_sn = S_etron/S_V_sn
  
  ; N_{SN}(e/V_sn) = C_sn/q
  
  N_sn = C_sn/q
  
  aName = '$S_{SN}(e^{-}/V_{SN})$
  
  ; Now the plots, using C for part C
  
  pl_S_V_sn = PLOT(V_sn, S_sn, XRANGE=[0.7,3.3], YLOG=1, YRANGE=[1e4,1.5e5], $ 
    XTITLE='Sense node Voltage (V$_{SN}$, Volts)', NAME='$S_{SN}(e^{-}$ $V^{-1}_{SN})$', $
    YTITLE='Signal and Noise Sen. (e$^{-}/V)$', $ 
    Title='Exercise 7.2 V/V non-linearity $\n$ Sense Node Signal & Noise Sensitivities', $
    THICK=2, FONT_NAME='Times', FONT_STYLE='Bold', FONT_SIZE=14, $
    XTICKLEN=1, YTICKLEN=1)
    
  pl_N_sn = PLOT(V_sn, N_sn, '-:', THICK=2, NAME='$N_{SN}(e^{-}$ $V^{-1}_{SN})$', /OVERPLOT)
     
  leg = LEGEND(TARGET=[pl_N_sn, pl_S_V_sn], POSITION=[2.5, 1.1e5], /DATA, $
              FONT_NAME='Times', FONT_STYLE='bold', FONT_SIZE=12) 

  diode_reset = TEXT(3.0,105000, reset, FONT_NAME='TIMES', FONT_STYLE='Bold', $
    FONT_SIZE=12, FONT_COLOR='red', /DATA)
  diode_disch = TEXT(0.8,105000, disch, FONT_NAME='TIMES', FONT_STYLE='Bold', $
    FONT_SIZE=12, FONT_COLOR='red', /DATA)

  pl_h = PLOT(Vh, [5,100000], 'r--', THICK=4, /OVERPLOT)
  pl_l = PLOT(Vl, [5,100000], 'r--', THICK=4, /OVERPLOT)

  
  ; Plot Sense node Gain (microVolts per electron)
  ; Inverse of Sensitivities are Gains (A)
  
  inv_S_sn = 1e6/S_sn
  inv_N_sn = 1e6/N_sn
  
  pl_micro = PLOT(V_sn, inv_S_sn, XRANGE=[0.7,3.3], YRANGE=[10,50], $
      XTITLE='Sense node Voltage (V$_{SN}$, Volts)', NAME='$S_{SN}$', $
      YTITLE='Signal and Noise Gain ($\mu$V/e$^{-})$', $
      Title='Exercise 7.2 V/V non-linearity $\n$ Sense Node Gains', $
      THICK=2, FONT_NAME='Times', FONT_STYLE='Bold', FONT_SIZE=14)
          
  pl_m2 = PLOT(V_sn, inv_N_sn, '-:', THICK=2, NAME='$N_{SN}$', /OVERPLOT)
          
  leg = LEGEND(TARGET=[pl_micro, pl_m2], POSITION=[1.8,45], /DATA, $
    FONT_NAME='Times', FONT_STYLE='bold', FONT_SIZE=12)
  
  diode_reset = TEXT(3.0,46, reset, FONT_NAME='TIMES', FONT_STYLE='Bold', $
    FONT_SIZE=12, /DATA)
  diode_disch = TEXT(0.8,46, disch, FONT_NAME='TIMES', FONT_STYLE='Bold', $
    FONT_SIZE=12, /DATA)

  pl_h = PLOT(Vh, [10,46], '--', THICK=2, /OVERPLOT)
  pl_l = PLOT(Vl, [10,46], '--', THICK=2, /OVERPLOT)
  
  
  ; Now use sensitivities S_sn and N_sn to plot Signal in Volts output 
  ; from the sense node [Signal(V_sn)] and the shot noise in volts 
  ; for the signal from the sense node [sigma_shot(V_sn)]
  
  ; Calculate the Sense node signal voltage S(V_sn) = S(e-trons)/S_sn(e/V_sn)
  ; and sigma_shot(V_sn) = sigma_shot(e-trons)/N_sn(e/V_sn)
  S_V_sn          = S_etron/S_sn   
  sigma_shot_V_sn = sigma_shot/N_sn
  
  ; Now plot the results to part 4
  
  pl_4 = PLOT(V_sn, S_V_sn, XRANGE=[0.7,3.3], YRANGE=[0.0005,8], $ 
       XTITLE='Sense node Voltage (V$_{SN}$, Volts)', NAME='$S(V_{SN})$', $
       YTITLE='$S(V_{SN}) and \sigma(V_{SN}$) (Volts)', XLOG=0,YLOG=1, $  
       TITLE='Exercise 7.2 V/V non-linearity $\n$ Sense Node Signal and Noise Voltages', $     
       THICK=2, FONT_NAME='Times', FONT_STYLE='Bold', FONT_SIZE=14)
     
  pl_42 = PLOT(V_sn, sigma_shot_V_sn, '-:', THICK=2, NAME='$\sigma_{SHOT}(V_{SN})$', $
       /OVERPLOT)
       
  leg=LEGEND(TARGET=[pl_4, pl_42], POSITION=[1.8,0.5], /DATA, $
       FONT_NAME='Times', FONT_STYLE='bold', FONT_SIZE=12)
       
   diode_reset = TEXT(3.0,3.5, reset, FONT_NAME='TIMES', FONT_STYLE='Bold', $
     FONT_SIZE=12, /DATA)
   diode_disch = TEXT(0.8,3.5, disch, FONT_NAME='TIMES', FONT_STYLE='Bold', $
     FONT_SIZE=12, /DATA)

   pl_h = PLOT(Vh, [0.0005,3], '--', THICK=2, /OVERPLOT)
   pl_l = PLOT(Vl, [0.0005,3], '--', THICK=2, /OVERPLOT)
   
   ; Results for part 5
   
   ; Gain A(DN/V_sn) = A_sf*A_cds*A_adc = 1.76e4 DN/V 
   ; A_adc for shorthand notation
   A_adc = 1.76e4
   
   ; S_adc(e/DN) = S_sn(e/V_sn)/A(DN/V_sn)
   S_adc = S_sn/A_adc
   
   ; N_adc(e/DN) = N_sn(e/V_sn)/A(DN/V_sn)
   N_adc = N_sn/A_adc
   
   ; K_adc(e/ADC) = S(DN)/[sigma_shot(DN)]^2 = A_adc(DN/V)*S(V_sn)/[A_adc(DN/V)*sigma_shot(V_sn)]^2
   K_adc = A_adc*S_V_sn/(A_adc*sigma_shot_V_sn)^2
   
   ; stop
      
   ; Plot away! This is versus Signal in electrons, S_etron
   
   pl_5 = PLOT(S_etron, K_adc, XRANGE=[0,90000.0], YRANGE=[1,5], $ 
        XTITLE='S($e^{-})$', NAME='$K_{ADC}$', $
        YTITLE='Sensitivity (e$^{-}/DN$)', XLOG=0, YLOG=0, $ 
        TITLE='Exercise 7.2 V/V non-linearity $\n$ ADC Sensitivities', $
        THICK=2, FONT_NAME='Times', FONT_STYLE='Bold', FONT_SIZE=14, $ 
        XTICKLEN=1, YTICKLEN=1, XMINOR=1, YMINOR=1)
        
   pl_52 = PLOT(S_etron, S_adc, '-:', THICK=2, NAME='S$_{ADC}$', $
        /OVERPLOT)
        
   pl_53 = PLOT(S_etron, N_adc, '__', THICK=2, NAME='N$_{ADC}$', $
        /OVERPLOT)
  
   leg = LEGEND(TARGET=[pl_5, pl_53, pl_52], POSITION=[25000,4], /DATA, $
        FONT_NAME='Times', FONT_STYLE='bold', FONT_SIZE=12)
        
   pl_fw = PLOT([84326.0,84326.0], [1,4.6], 'r--', THICK=2, /OVERPLOT)
   
   full_well = TEXT(78000.0, 4.7, 'Full Well', FONT_NAME='TIMES', FONT_STYLE='Bold', $
     FONT_SIZE=12, FONT_COLOR='red', /DATA)
  
   ; Calculate and plot S(DN) the signal in DN and sigma(DN) the noise in DN 
   ; versus S signal in electrons. Note this is linear and does not include 
   ; non-linear S_adc and N_adc effects. At this stage in Photon Transfer 
   ; the calculation of the transfer is S(V_SN) -> S(DN)
   
   S_DN = A_adc*S_V_sn
   
   sigma_shot_DN = A_adc*sigma_shot_V_sn
   
   ; Plot 'em
   
   pl_SDN = PLOT(S_etron, S_DN, 'b', XRANGE=[150,1.1e5], YRANGE=[10,45000], $
     XTITLE='S($e^{-})$', NAME='$S(DN)$', $
     YTITLE='DN', XLOG=1, YLOG=1, $
     TITLE='Exercise 7.2 V/V non-linearity $\n$ Linear DN Signal and Noise', $
     THICK=3, FONT_NAME='Times', FONT_STYLE='Bold', FONT_SIZE=14, $
     XTICKLEN=1, YTICKLEN=1,XMINOR=4, YMINOR=4)
   
   pl_sigDN = PLOT(S_etron, sigma_shot_DN, 'b-:', THICK=3, NAME='$\sigma(DN)$', $
     /OVERPLOT)
     
   leg = LEGEND(TARGET=[pl_SDN, pl_sigDN], POSITION=[1000, 10000], /DATA, $
     FONT_NAME='Times', FONT_STYLE='bold', FONT_SIZE=12)
     
   ; Plot False signal (electrons) based on sensitivity K_adc(e_tron/DN) 
   ; versus true signal (electrons)
   
   false_S_etron = S_DN*K_adc 
   
   ; Make the plot showing false and true signal
   
   pl_falseSig = PLOT(S_etron, false_S_etron, 'b', XRANGE=[150,1.0e5], $ 
     YRANGE=[0.1,3.5e5],  XTITLE='S($e^{-})$', NAME='$S_{False}$', $
     YTITLE='False signal (Non-linear $-$ linear)', $
     TITLE='Exercise 7.2 V/V non-linearity $\n$ False & True Signal ($e^{-}$)', $
     THICK=3, FONT_NAME='Times', FONT_STYLE='Bold', FONT_SIZE=14, $
     XTICKLEN=1, YTICKLEN=1,XMINOR=1, YMINOR=1, XLOG=0, YLOG=0)
   
   pl_trueSig = PLOT(S_etron, S_etron, 'b--', THICK=3, NAME='$S_{True}$', $
     /OVERPLOT)
     
   ; Overplot the Full well
   pl_fw = PLOT([84326.0,84326.0], [0,3.3e5], 'r--', THICK=2, /OVERPLOT)  

   full_well = TEXT(84310.0, 3.33e5, 'Full Well', FONT_NAME='TIMES', FONT_STYLE='Bold', $
     FONT_SIZE=12, FONT_COLOR='red', /DATA)
     
   leg = LEGEND(TARGET=[pl_falseSig, pl_trueSig], POSITION=[4e4, 2e5], /DATA, $
     FONT_NAME='Times', FONT_STYLE='bold', FONT_SIZE=12, COLOR='blue')
     
   ;; LAST ONE FOR 7-2! Non-linearity 
   ; Non-lin for K_adc(e-/DN)
   
   nl_K_adc = 100*(K_adc - MIN(K_adc))/MIN(K_adc)
   
   pl_nl = PLOT(S_etron, nl_K_adc, 'b', XRANGE=[150,1.0e5], $ 
    YRANGE=[0,100], XTITLE='S($e^{-})$', NAME='$NL%(K_{ADC})$', $
    YTITLE='Non-linearity (%)', $
    TITLE='Exercise 7.2 V/V non-linearity $\n$ Percent Non-linearity', $
    THICK=3, FONT_NAME='Times', FONT_STYLE='Bold', FONT_SIZE=14, $
    XTICKLEN=1, YTICKLEN=1,XMINOR=1, YMINOR=1, XLOG=1, YLOG=0)
    
   ; Non-lin for S_sn(e/DN). I call this S_adc here. 
    
   nl_S_adc = 100*(S_adc - MIN(S_adc))/MIN(S_adc)
    
   pl_nl2 = PLOT(S_etron, nl_S_adc, 'b--', THICK=3, NAME='$NL%(S_{SN})$', $
     /OVERPLOT)
     
   ; Non-lin for N_sn(e/DN). I call this N_adc here.

   nl_N_adc = 100*(N_adc - MIN(N_adc))/MIN(N_adc)

   pl_nl3 = PLOT(S_etron, nl_N_adc, 'b-:', THICK=3, NAME='$NL%(N_{SN})$', $
     /OVERPLOT)
     
   leg = LEGEND(TARGET=[pl_nl, pl_nl3, pl_nl2], POSITION=[5.5e3, 80], /DATA, $
     FONT_NAME='Times', FONT_STYLE='bold', FONT_SIZE=12, COLOR='blue')
   
   ; Overplot the Full well
   pl_fw = PLOT([84326.0,84326.0], [0,90], 'r--', THICK=2, /OVERPLOT)

   full_well = TEXT(65000.0, 92, 'Full$\n$Well', FONT_NAME='TIMES', FONT_STYLE='Bold', $
     FONT_SIZE=12, FONT_COLOR='red', /DATA)
   
     
END

;**********

;**********
; Chapter 8
;**********

;****************************************************

PRO ex_8_2

; Generate a set of PTC's after flat fielding for the PTC's in 
; Example 5-1. Assume Q_FF = 10^2, 10^3, 10^4, 10^5, and 10^6
; electrons and N_FF = 1. Calculate the flat field level 
; where the corrected noise is greater than the raw noise. 
; Assume the FPN Quality Factor P_N = 0.02 or 2%. 

; USES FUNCTION my_wien


  xs = DINDGEN(400000, START = 1.0d)
  ; let ys all = 1.0
  ys = xs/xs
  
  p1 = PLOT(xs, ys, XLOG=1, YLOG=1, XRANGE=[1,1e6], YRANGE=[1,1e4], $
    XTITLE='SIGNAL (e$^{-}$)', YTITLE='NOISE (e$^{-}$)', FONT_SIZE=12, $
    FONT_STYLE='Bold', FONT_NAME='Times', $ 
    Title = 'Example 8.2.1 $\sigma_{COR} and Q_{FF} = 100$ with $\Delta$', $
    XTHICK=2, YTHICK=2, /NODATA)

  ; Plot the read noise
  
  ; sigma_read = 3.33*ys
  ; From Fig. 8.4 Readnoise = 5 electrons         
  sigma_read = 5.0*ys
  
  p2 = PLOT(xs, sigma_read, NAME='$\sigma_{RN}$', THICK=2, /OVERPLOT)

  ; Plot the shot noise
  ; Minimum signal is 1 electron generates 1.0 DN, Fig. 8.4 is in electrons

  K_adc = 1.0

  sigma_shot = SQRT(xs/K_adc)  ; Eqn. 5.7 = Eqn. 4.20

  ; Modify for full well

  ; fw = WHERE(xs GE 233000)
  fw = WHERE(xs GE 220000)

  ; Call my_wien for the full well values
  sigma_shot[fw] = my_wien(xs[fw],sigma_shot[fw],233000)


  p3 = PLOT(xs, sigma_shot, NAME='$\sigma_{SHOT}(DN)$', THICK=2, LINESTYLE=2, $
          /OVERPLOT)

  ; Plot the Flat-fielded corrected noise sigma_COR

  ; set the FPN Quality Factor P_n
  P_n = 0.02

  ; Set the number of Flat Field images
  N_ff = 1

  ; Create an array to hold the Flat Field Quality Factor
  ; Q_ff values, recall Q_ff = (S_ff)*(N_ff)
  
  Q_ff = [1e02, 1e03, 1e04, 1e05, 1e06]*N_ff
  
  sigma_cor_1 = SQRT(sigma_read^2 + (xs*(1 + xs/Q_ff[0])))
  
  sigma_cor_2 = SQRT(sigma_read^2 + (xs*(1 + xs/Q_ff[1])))
  sigma_cor_3 = SQRT(sigma_read^2 + (xs*(1 + xs/Q_ff[2])))
  sigma_cor_4 = SQRT(sigma_read^2 + (xs*(1 + xs/Q_ff[3])))
  sigma_cor_5 = SQRT(sigma_read^2 + (xs*(1 + xs/Q_ff[4])))
  
  ; Apply the full well turnover to all the sigma_cor's
  sigma_cor_1[fw] = my_wien(xs[fw],sigma_cor_1[fw],233000)
  sigma_cor_2[fw] = my_wien(xs[fw],sigma_cor_2[fw],233000)
  sigma_cor_3[fw] = my_wien(xs[fw],sigma_cor_3[fw],233000)
  sigma_cor_4[fw] = my_wien(xs[fw],sigma_cor_4[fw],233000)
  sigma_cor_5[fw] = my_wien(xs[fw],sigma_cor_5[fw],233000)
  

  diff_1 = sigma_cor_1 - sigma_shot
  
  min_diff_1 = MIN(diff_1, min_loc)
  
  print,'**********'
  print,' Minimum difference sigma_cor Q_ff = 100 minus sigma_shot: ' 
  print,' = '+STRTRIM(STRING(min_diff_1),1)+' at Signal = '+ $ 
          STRTRIM(STRING(xs[min_loc]),1)+' electrons.'
  print,'**********'
  
  
  p4 = PLOT(xs, sigma_cor_1, NAME='$\sigma_{COR}(Q_{ff}=100)$',THICK=2, $ 
            LINESTYLE=1, /OVERPLOT)

  ; Plot the delta sigma_cor_1 - sigma_shot
  
  p5 = PLOT(xs, diff_1, NAME='$\Delta(\sigma_{COR} - \sigma_{SHOT})$', THICK=2, $ 
           COLOR='red', /OVERPLOT)
  
  ; Plot the total noise for Q_ff = 100
  
  ; Eqn. 5.9
  sigma_total = SQRT(sigma_read^2 + sigma_shot^2 + sigma_cor_1^2)

  p6 = PLOT(xs, sigma_total, NAME='$\sigma_{TOTAL}(^{-})$', THICK=2, LINESTYLE=4, $
    /OVERPLOT)

  ; Plot the curves for Q_ff > 100
  p7 = PLOT(xs, sigma_cor_2, NAME='$\sigma_{COR}(Q_{ff}=1000)$',THICK=2, $
    LINESTYLE=1, /OVERPLOT)
  p8 = PLOT(xs, sigma_cor_3, NAME='$\sigma_{COR}(Q_{ff}=10^{4})$',THICK=2, $
    LINESTYLE=1, /OVERPLOT)
  p9 = PLOT(xs, sigma_cor_4, NAME='$\sigma_{COR}(Q_{ff}=10^{5})$',THICK=2, $
    LINESTYLE=1, /OVERPLOT)
  p10 = PLOT(xs, sigma_cor_5, NAME='$\sigma_{COR}(Q_{ff}=10^{6})$',THICK=2, $
    LINESTYLE=1, /OVERPLOT)

  ; Add the LEGEND to the plot

  leg = LEGEND(TARGET=[p2,p3,p4,p5,p6], POSITION=[200,4000], /DATA)

  ; Now re-produce the Exercise solution, Figure 8.4 in the text

  pa = PLOT(xs, ys, XLOG=1, YLOG=1, XRANGE=[1,1e6], YRANGE=[1,1e4], $
    XTITLE='SIGNAL (e$^{-}$)', YTITLE='NOISE (e$^{-}$)', FONT_SIZE=12, $
    FONT_STYLE='Bold', FONT_NAME='Times', $
    Title = 'Example 8.2.2 $\sigma_{COR}$ Solution', $
    XTHICK=2, YTHICK=2, /NODATA)

  ; Overplot readnoise and shot noise
  pb = PLOT(xs, sigma_read, NAME='$\sigma_{RN}$', THICK=2, /OVERPLOT)
  pc = PLOT(xs, sigma_shot, NAME='$\sigma_{SHOT}(DN)$', THICK=2, /OVERPLOT)

  ; Calculate and plot the FPN noise sigma_fpn = P_n*signal
  sigma_fpn = P_n*xs
  ; Allow for the full well in the FPN curve
  sigma_fpn[fw] = my_wien(xs[fw],sigma_fpn[fw],233000)
  
  pd = PLOT(xs, sigma_fpn, NAME='$\sigma_{FPN}$', THICK=2, /OVERPLOT)
  
  ; The total, ideal noise
  
  sigma_tot_ideal = SQRT(sigma_read^2 + sigma_shot^2 + sigma_fpn^2)
  
  pe = PLOT(xs, sigma_tot_ideal,'r', NAME='$\sigma_{TOT}$(Ideal)', THICK=2, $
           /OVERPLOT)
           
  ; Plot the curves for Q_ff >= 100
  p4 = PLOT(xs, sigma_cor_1, NAME='$\sigma_{COR}(Q_{ff}=100)$',THICK=2, $
    LINESTYLE=1, /OVERPLOT)
  p7 = PLOT(xs, sigma_cor_2, NAME='$\sigma_{COR}(Q_{ff}=1000)$',THICK=2, $
    LINESTYLE=3, /OVERPLOT)
  p8 = PLOT(xs, sigma_cor_3, NAME='$\sigma_{COR}(Q_{ff}=10^{4})$',THICK=2, $
    LINESTYLE=4, /OVERPLOT)
  p9 = PLOT(xs, sigma_cor_4, NAME='$\sigma_{COR}(Q_{ff}=10^{5})$',THICK=2, $
    LINESTYLE=5, /OVERPLOT)
  p10 = PLOT(xs, sigma_cor_5, 'g', NAME='$\sigma_{COR}(Q_{ff}=10^{6})$',THICK=2, $
    LINESTYLE=1, /OVERPLOT)
    
  ; Some TEXT annotations
  
  t_rn   = TEXT(40000, 6.0, '$\sigma_{RN}$', /DATA, FONT_SIZE=10)
  t_fpn  = TEXT(900, 10.5, '$\sigma_{FPN}$', /DATA, FONT_SIZE=10)
  t_shot = TEXT(40000, 140, '$\sigma_{SHOT}$', /DATA, FONT_SIZE=10)
  
  ; Add a LEGEND for the Q_ff's to this plot

  leg = LEGEND(TARGET=[p4,p7,p8,p9,p10], POSITION=[200,4000], /DATA)


END