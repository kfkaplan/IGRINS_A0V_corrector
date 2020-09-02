# The IGRINS A0V Corrector
This program is designed to correct the A0V spectrum flexure and airmass to match the science target in order to improve telluric correction.
- It applies a 0th order correction pixel shift in the detector's x-direction to correct for the difference flexure in the detector's x direction between when.
- It also applies a power-law correction to account for the difference in airmass betweem when the A0V star was observed and when the science target was observed. 
- Both corrections are fit iteratively
- The correction is applied to the A0V spectrum and then saved as corrected fits files

Questions, comments, want to contribute?
Email me at kkaplan@usra.edu

# Requirements
- Python 3.7 or above
- Astropy
  -  If you have anaconda on the command line type `conda install astropy`
  -  If you have python installed some other way, on the command line type `pip install -U astropy`
- IGRINS pipeline output `.spec.fits` and `.spec_flattened.fits` files for your A0V standard star and the `.spec.fits` file for your science target

# Installation
Download or clone the git repo to your machine.  

# Set up and run
- Identify night and the frame numbers for the A0V star and science target you want to correct between.
- Copy `a0v_corrector.py` to the IGRINS pipeline `plp/output/YYYYMMDD` directory (or directories) for the night(s) you have data you want to apply the correction to.
- Open that copy of `a0v_corrector.py` in your favorite text editor.
- Modify the file paths to point the the `.spec.fits` and `.spec_flattened.fits` file for the A0V star spectrum and the `.spec.fits` file for the science spectrum. Note that the H and K bands are treated as seperate spectra so you will need to do this process twice, once for each band. In the example below the night will be 20200901, the A0V frame number will be 0093, the science target frame number will be 0105, and we will use the H-band spectra (files start with `SDCH_`, for the K-band they start with `SDCK_`).
 ```
  #--------------------------------------------User Inputs (modify these)-----------------------------------------------------------------------------------------------
a0v_fits_path = 'SDCH_20200901_0093.spec.fits' #Path to A0V .spec.fits file
flattened_a0v_fits_path = 'SDCH_20200901_0093.spec_flattened.fits' #Path to A0V .spec_flattened.fits file
sci_fits_path = 'SDCH_20200901_0093.spec.fits' #Path to science target .spec.fits file
  ```
  - You can manually enter your corrections or let the program try to automatically determine the corrections.  Either way, you need to define the input corrections  `input_flexure` and `input_power` which will be fixed if manual or initial guesses if automatic.  Then you need to define if you want to let the program automatically find the corrections or not.  If this is your first time using this program and you want to let the program run automatically, just leave the values alone.
```
input_flexure = 0.0 #Initial guess (or manually set) flexure correction
input_power = 1.0 #Initial guess (or manually set) powerlaw (airmass) correction 
automate_fit = True #Set True for automatic fitting, set False for a simple manual correction (user must supply manual correction in the two variables above)
```
- You can modify the parameters for the automated fitting.  In general these are best left unchanged, unless you know what you are doing.  Please contact me with any questions about these parameters at kkaplan@usra.edu.
```
#-------------------------------Automation Inputs (modify these if you know what you are doing)------------------------------------------------------------------------
n_iterations = 5 #Number of iterations
skip_first_orders = 4 #Skip how many first and last orders?
skip_last_orders = 4
flexure_array = np.arange(-2.0, 2.0 + 0.2, 0.2) #Define the range and step size in pixels to cross-correlate the flexure between the a0v and science spectra
power_array = np.arange(0., 2.0 + 0.2, 0.2) #Define the range of powers and step size to corss-correlate the telluric line depths between the a0v and science spectra
min_x, max_x = [300, 1750] #Cuts to sides pixels of detector in x-direction
min_a0v, max_a0v = [0.2, 0.85] #Cuts to continuum normalized A0V spectrum (we are interested only in the telluric line edges)
```
- To run simply call the python script from the command line
```
python a0v_corrector.py
```
- Next repeat the above process for the K-band by changing SDCH to SDCK in the file paths.
```
  #--------------------------------------------User Inputs (modify these)-----------------------------------------------------------------------------------------------
a0v_fits_path = 'SDCK_20200901_0093.spec.fits' #Path to A0V .spec.fits file
flattened_a0v_fits_path = 'SDCK_20200901_0093.spec_flattened.fits' #Path to A0V .spec_flattened.fits file
sci_fits_path = 'SDCK_20200901_0093.spec.fits' #Path to science target .spec.fits file
```
- Run the python script again to process the K-band
```
python a0v_corrector.py
```
# Outputs

The A0V star spectrum corrected for flexure and airmass it outputted as `.spec.corrected.fits` and `.spec_flattened.corrected.fits` files.  These files are indentical to the original files except the spectra in the first extension have had the corrections applied to them.  You will want to back up your original `.spec.fits` and `.spec_flattened.fits` files and use the `.spec.corrected.fits` and `.spec_flattened.corrected.fits` files instead.  There are no outputs for the science spectrum since the A0V spectrum's flexure and airmass is corrected to match the science spectrum.

The `.pdf` files show the correction itself.  If you run the automated correction, the first two pages will be the chisq plots for the flexure and power (airmass).  The rest of the pages show each echelle order of the science spectrum telluric corrected by the A0V spectrum before and after the flexure and airmass corrections have been applied so you can see how good the corrections are.

For the example above,  `a0v_corrector.py` would output the following files
```
SDCH_20200901_0093.spec.corrected.pdf
SDCH_20200901_0093.spec.corrected.fits
SDCH_20200901_0093.spec_flattened.corrected.fits
SDCK_20200901_0093.spec.corrected.pdf
SDCK_20200901_0093.spec.corrected.fits
SDCJ_20200901_0093.spec_flattened.corrected.fits
```
