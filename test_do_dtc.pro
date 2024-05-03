; An IDL batch file to run a regression test on the Photon Transfer Code
; plot_ptc.pro

; TO RUN THIS TEST:
; 1) Start idl at the terminal 
; 2) at the IDL prompt type the following and Enter
; 
; IDL> @test_do_dtc

; THe cretor (Bob Zavala) typically opened the code in the IDL Developement 
; Environment and compiled it after the directories were changed as needed.
; Within the directory:
; 
; /mnt/nofs/projects/solarSystemEphem/images/dark-dtc/testCode
; 
;   this batch file may be run if the IDLDE is not used. 

; Compile the PTC code
.r /mnt/nofs/projects/solarSystemEphem/pro/plot_ptc.pro

; CD to the test directory
CD, '/mnt/nofs/projects/solarSystemEphem/images/dark-dtc/testCode'

; run the test
do_dtc

; Output information to user
print,' '
print,'********************'
print,' Test of procedure do_dtc completed.'
print,' Compare the two plots that were created to the plots in this directory:'
print,' '
print,'     multiNight_test_do_dtc.png'
print,'     multiNight_test_do_dtc_darkRate.png'
print,' '
print,' The plots just created and the comparisons should agree.'
print,'********************'
print,' '
