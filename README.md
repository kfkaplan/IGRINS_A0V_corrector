# IGRINS_A0V_corrector
Program for IGRINS spectra to correct differences in flexure and airmass between the A0V Standard Star and science target.

# Contact
Questions, comments, want to contribute?
Email me at kkaplan@usra.edu

# What does it do?
This program is designed to do two things to correct the A0V spectrum to better match the science target spectrum to improve your telluric corrections and/or relative flux calibration.  The first thing it does is apply a 0th order correction to the flexure in the detector's x direction.  The second thing is does is apply a power-law correction to account for the difference in airmass betweem when the A0V star was observed and when the science target was observed.  

# Requirements
- Python 3.7 or above (I reccomend installing with anaconda)
- Astropy
  -  If you have anaconda on the command line type `conda install astropy`
  -  If you have python installed some other way, on the command line type `pip install -U astropy`

# Installation
Download or clone the git repo to your machine.  Copy `a0v_corrector.py` to the IGRINS pipeline `plp/output/YYYYMMDD` directory (or directories) for the night(s) you have data you want to apply the correction to.

# Set up
Identify the frame numbers for the A0V star and science target you want to correct between
