# CL script lincheck.cl 
# Bob Zavala me fecit

# An IRAF CL (command language) script for processing a night
# of CCD Photon Transfer Curve (PTC) linearity data taken at
# the 61-inch telescope. Illumination provided by flat-field lamps
# on the dome and these flats and biases are present in the directory.

# REQUIRES:
#   1) Script calls a Bash shell columnSep.sh and expects this script
#   to be in the same sirectory as lincheck.cl
#
#   2) The CCD Database in IRAF is setup. Do you need help?
#   Bob Zavala can help with this, or check the IRAF Community
#   online documentation. Once setup CCD, CMOS and IR reducitions
#   are streamlined. 

# RUNS ON:
#   Runs on Community IRAF V2.17.1 with a run time on a NOFS
#   VM of 36 seconds for 348 Princeton Sophia 4096 x 4096
#   CCD images.
#
#   Tested running within the directory where the PTC images reside. 

# INPUT
#   1) The images in a single directory in the usual NOFS naming
#   convention g??.???.??? are the only input. No argument is passed
#   as the script looks for filenames with this pattern. 

# OUTPUT
#   A file called 'raw.signal.txt' containing 10 columns useful for
#   subsequent PTC analysis. The columns are:
#
#	1)  Image filename
#	2)  Image section used for statistics
#	3)  Mean pixel value in that section (DN)
#	4)  Median pixel value in that section (DN)
#	5)  Std. deviation in that section (DN)
#	6)  Minimum pixel value in that section (DN)
#	7)  Maximum pixel value in that section (DN)
#	8)  Exposure time in seconds
#	9)  UTC time hh:mm:ss.sss of the esposure start
#	10) OBSTYPE header keyword (e.g. BIAS, Dome_Flat) 

# At start tell the user the data and time. Report this again
# at the end so if interested a quick eyeball time estimate of the
# runtime is available. 
print "\n"
date

# Remind the user of the name of the script ;)
print "\n"
print"Running script lincheck.cl\n"

# Check to see if the ccdred IRAF package is loaded.
# If not, load it now.

if (deftask ("ccdred"))
   print ("Package ccdred loaded, continuing.\n")
else
   print ("I need to load the ccdred package.\n")
   noao
   imred
   ccdred

# Run the ccdlsit task as this is useful to have around
# even if not needed. THis part requires the CD Database
# mentioned above in "preliminaries".
ccdlist ("g23d???.???.fits", ccdtype = "", names=no,
long=no,ccdproc="", > "ccdlist.txt")

print "File ccdlist.list created.\n"

# The object type e.g. Bias, Dark or Flat is required
# for the subsequent analysis.
hedit ("g23d???.???.fits",
"OBJECT", ".", add=no, addonly=no, delete=no, verify=no,
show=yes, update=no, > "object.list")

print "File object.list created.\n"

# This copies a one-line file in the parent directory
# so a "column header" is added. I could add this with an
# editing command, but I like this method. 
cp -v /mnt/nofs/projects/solarSystemEphem/images/obstype.txt .

# Becuase of compatibility issues with header keywords
# as the NOFS headers evolved thw OBSTYPE keyword is needed. 
hselect ("g23d???.???.fits",
"OBSTYPE", "yes", missing="INDEF", >> "obstype.txt")

print "File obstype.list created.\n"

# Get the image statistice required for the PTC analysis.
# The region is chosen for the Sophia CCD. A different
# camera is likely to have a different region.
### In future I should change this so the region is an argument.
### or maybe a camera and a lookup table is provided for
### different cameras. 
imstat ("g23d???.???.fits[900:3200,1200:3300]",
fields="image,npix,mean,midpt,stddev,min,max", lower=INDEF, upper=INDEF,
nclip=0, lsigma=3., usigma=3., binwidth=0.1, format=yes, cache=no,
> "raw.signal.imstat.txt")

print "File raw.signal.imstat.txt created.\n"

# Another column header line copy.
cp -v /mnt/nofs/projects/solarSystemEphem/images/exptimes.txt .

# As you might imagine, exposure times are also needed
# for the PTC analysis. 
hselect ("g23d???.???.fits",
"EXPTIME", "yes", missing="INDEF", >> "exptimes.txt")

print "File exptimes.txt created.\n"

# The UTC time of observation is a useful parameter
# to extract for the PTC analysis. Date is already
# coded in the filename. 
cp -v /mnt/nofs/projects/solarSystemEphem/images/utcobs.txt .

hselect ("g23d???.???.fits",
"UTC-OBS", "yes", missing="INDEF", >> "utcobs.txt")

print "File utcobs.txt created.\n"

# Use the existing Unix 'paste' command to stitch everything
# together into a file that will be sued for the PTC analysis.
!paste raw.signal.imstat.txt exptimes.txt utcobs.txt obstype.txt > raw.signal.txt

# Call the Shell script that uses the tried and true ed editor
# to add one needed whitespace so we get the columns we want
# that will be read in for the subsequent analysis.
!./columnSep.sh

# Miller time. 
print "\n"
print "Script lincheck.cl finished\n"
print "File raw.signal.txt ready for PTC analysis\n"
print "\n"

# Again pass the system date and time so the eyeball runtime estimate may
# be made. If you really want a more precise estiamte to be made, feel
# free to improve upon this flawless piece of code :) 
date

# End cleanly from the script
bye

