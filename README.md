# photonTransfer

Code in various languages that I wrote that is helpful for Photon Transfer Curve (PTC) analysis of actual experimental data. I have some different code written for simulating PTC analysis. 

For learning about PTC see:

 "Photon Transfer DN -> Lambda" by James R. Janesick published by the SPIE Vol. No. PM170 (2007) 
ISBN 9780819467225

Brief descriptions:

lincheck.cl = An IRAF Command Language (CL) script

columnSep.sh = A shell script called by lincheck.cl

plot_ptc.pro = IDL code containing procedures to perform the PTC analysis. This is only a barebones beginning at this stage. 

output-lincheck.txt = Output captured from running lincheck.cl on Community IRAF 2.17.1 on some example data. 

examples.pro = IDL code used for solutions to examples in Janesick's book "Photon Transfer DN -> Lambda."
