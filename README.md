# IGRINS_A0V_corrector
Program for IGRINS spectra to correct differences in flexure and airmass between the A0V Standard Star and science target.

# Contact
Questions, comments, want to contribute?
Email me at kkaplan@usra.edu

# What does it do?
This program is designed to correct the A0V spectrum flexure and airmass match the science target, to improve telluric correction and/or relative flux calibration.
- It applies a 0th order correction pixel shift in the detector's x-direction to correct for the difference flexure in the detector's x direction between when.
- It also applies a power-law correction to account for the difference in airmass betweem when the A0V star was observed and when the science target was observed. 
- Both corrections are fit iteratively
- The correction is applied to the A0V spectrum and then saved as corrected fits files

# Requirements
- Python 3.7 or above (I reccomend installing with anaconda)
- Astropy
  -  If you have anaconda on the command line type `conda install astropy`
  -  If you have python installed some other way, on the command line type `pip install -U astropy`

# Installation
Download or clone the git repo to your machine.  

# Set up and run
- Identify night and the frame numbers for the A0V star and science target you want to correct between
- Copy `a0v_corrector.py` to the IGRINS pipeline `plp/output/YYYYMMDD` directory (or directories) for the night(s) you have data you want to apply the correction to.
- Open that copy of `a0v_corrector.py` in your favorite text editor
- Modify the file paths to match the the `.spec.fits` and `.spec_flattened.fits` file for the A0V star spectrum and the `.spec.fits` file for the science spectrum
