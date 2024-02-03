#! /bin/sh

## Shell script columnSep.sh
## Bob Zavala me fecit

## A shell script that is called with the Community
## IRAF CL script lincheck.cl to add one whitespace
## between the image filename and the pixel statistics
## section in the raw.signal.txt file. The raw.signal.txt
## is the output of lincheck.cl and is used in Photon
## Transfer Curve (PTC) analysis. 

## REQUIRES:
##   1) The IRAF CL script lincheck.cl is called within IRAF,
##   and this script resides in the same directory as lincheck.cl.
##
##   2) This script must be executable by the user. 

## INPUT:
##   1) File called 'raw.signal.txt' and this filename
##   is written within the script.

## OUTPUT:
##   2) File called 'raw.signal.txt', same name as the input
##   This script adds a single column of whitespace so there
##   seems little need to create a new filename and clog up
##   the directory. 

## A script using the ed editor to separate the 
## filename from the statistics section from the 
## imstatistics output so that these will easily be
## read as individual columns becasue they will be
## separated by white space.

echo " "
echo "  Using ed to edit all lines in the raw.signal.txt file.\n"
echo "  This will change 'fits[' -> 'fits [' \n"
echo "  as needed to separate columns in the file.\n"
ed -s raw.signal.txt <<EOF
!echo "  In the ed editor\n"
,s/fits/fits /g
,wq
EOF
echo "  File raw.signal.txt file imstat columns separated.\n"
