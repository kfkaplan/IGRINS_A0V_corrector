# This script takes the outputs of the IGRINS PLP and corrects for flexure and differneces in airmass between the standard A0V
# stars and the science target.  For automatif fitting, it assumes the science target has enough continuum to see and fit telluric lines.


#--------------------------------------------User Inputs (modify these)-----------------------------------------------------------------------------------------------
a0v_fits_path =  #Path to A0V .spec.fits file
flattened_a0v_fits_path = #Path to A0V .spec_flattened.fits file
sci_fits_path = #Path to science target .spec.fits file
input_flexure = 0.0
input_power = 1.0
automate_fit = True
#-------------------------------Automation Inputs (modify these if you know what you are doing)------------------------------------------------------------------------
n_iterations = 5 #Number of iterations
skip_first_orders = 4 #Skip how many first and last orders?
skip_last_orders = 4
flexure_array = np.arange(-2.0, 2.0 + 0.2, 0.2) #Define the range and step size in pixels to cross-correlate the flexure between the a0v and science spectra
power_array = np.arange(0., 2.0 + 0.2, 0.2) #Define the range of powers and step size to corss-correlate the telluric line depths between the a0v and science spectra
min_x, max_x = [300, 1750] #Cuts to sides pixels of detector in x-direction
min_a0v, max_a0v = [0.2, 0.85] #Cuts to continuum normalized A0V spectrum (we are interested only in the telluric line edges)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Import necessary libraries
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.ndimage import binary_dilation, binary_erosion
import numpy as np
#import os
import copy
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages  #For outputting a pdf with multiple pages (or one page)
#import pdb
try:  #Try to import bottleneck library, this greatly speeds up things such as nanmedian, nanmax, and nanmin
	import bottleneck as bn #Library to speed up some numpy routines
except ImportError:
	print("Bottleneck library not installed.  Code will still run but might be slower.  You can try to bottleneck with 'conda install bottleneck' or 'pip install bottleneck' for a speed up.")
	import numpy as bn
import matplotlib.gridspec as grd
np.seterr(all='ignore') #Turn off those annoying runtime warnings




#Define some useful variables
n_flexure = len(flexure_array)
n_power = len(power_array)

#Do a robost running median filter that ignores nan values and outliers, returns result in 1D
def robust_median_filter(flux, size=375):
	if size%2 == 0: size = size+1 #Make even results odd
	half_sizes = np.array([-(size-1)/2, ((size-1)/2)+1], dtype='int')		
	if np.ndim(flux) == 2: #For 2D spectrum
		ny, nx = np.shape(flux) #Calculate npix in x and y
	else: #Else for 1D spectrum
		nx = len(flux) #Calculate npix
	median_result = np.zeros(nx) #Create array that will store the smoothed median spectrum
	if np.ndim(flux) == 2: #Run this loop for 2D
		for i in range(nx): #This loop does the running of the median down the spectrum each pixel
			x_left, x_right = i + half_sizes
			if x_left < 0:
				x_left = 0
			elif x_right > nx:
				x_right = nx
			median_result[i] = bn.nanmedian(flux[:,x_left:x_right]) #Calculate median between x_left and x_right for a given pixel
	else: #Run this loop for 1D
		for i in range(nx): #This loop does the running of the median down the spectrum each pixel
			x_left, x_right = i + half_sizes
			if x_left < 0:
				x_left = 0
			elif x_right > nx:
				x_right = nx
			median_result[i] = bn.nanmedian(flux[x_left:x_right])  #Calculate median between x_left and x_right for a given pixel
	return median_result


def normalize_continuum(input_flux_array, size=375): #Normalize spectrum to continuum using robust running median filter
		median_result_1d = robust_median_filter(input_flux_array, size=size) #Take a robust running median along the trace, this is the found continuum
		normalized_flux = input_flux_array / median_result_1d #Normalize continuum 
		return(normalized_flux) #Replace this order's flux array with one that has been continuum normalized

#Definition takes a high resolution spectrum and rebins it (via interpolation and integration) onto a smaller grid
#while conserving flux, based on Chad Bender's "srebin" which is in turn based on "TERRASPEC_CALCULATE_REBIN"
def srebin(oldWave, newWave, oldFlux, kind='linear'):
	nPix = len(newWave) #Number of pixels in new binned spectrum
	newFlux = np.zeros(len(newWave)) #Set up array to store rebinned fluxes
	interpObj = interp1d(oldWave, oldFlux, kind=kind, bounds_error=False) #Create a 1D linear interpolation object for finding the flux density at any given wavelength
	#wavebindiffs = newWave[1:] - newWave[:-1] #Calculate difference in wavelengths between each pixel on the new wavelength grid
	wavebindiffs = np.diff(newWave) #Calculate difference in wavelengths between each pixel on the new wavelength grid
	wavebindiffs = np.hstack([wavebindiffs, wavebindiffs[-1]]) #Reflect last difference so that wavebindiffs is the same size as newWave
	wavebinleft =  newWave - 0.5*wavebindiffs #Get left side wavelengths for each bin
	wavebinright = newWave + 0.5*wavebindiffs #get right side wavelengths for each bin
	fluxbinleft  = interpObj(wavebinleft)
	fluxbinright = interpObj(wavebinright)
	for i in range(nPix): #Loop through each pixel on the new wavelength grid
		useOldWaves = (oldWave >= wavebinleft[i]) & (oldWave <= wavebinright[i]) #Find old wavelength points that are inside the new bin
		nPoints = bn.nansum(useOldWaves)
		wavePoints = np.zeros(nPoints+2)
		fluxPoints = np.zeros(nPoints+2)
		wavePoints[0] = wavebinleft[i]
		wavePoints[1:-1] = oldWave[useOldWaves]
		wavePoints[-1] = wavebinright[i]
		fluxPoints[0] = fluxbinleft[i]
		fluxPoints[1:-1] = oldFlux[useOldWaves]
		fluxPoints[-1] = fluxbinright[i]
		newFlux[i] =  0.5 * bn.nansum((fluxPoints[:-1]+fluxPoints[1:])*np.diff(wavePoints)) / wavebindiffs[i]
	return newFlux



hdul_a0v = fits.open(a0v_fits_filename) #Open fits files
hdul_flattened_a0v = fits.open(flattened_a0v_fits_filename)
hdul_sci = fits.open(sci_fits_filename)
data_a0v = hdul_a0v[0].data #Grab data from each fits file (for IGRINS, the data are stored in the zeroith extension)
data_a0v_flattened =  hdul_flattened_a0v[0].data
data_sci = hdul_sci[0].data
n_orders, n_x = np.shape(data_a0v) #Get the number of orders and pixels in x
x = np.arange(n_x) #Create an array to store the pixel x values (across the detector) for each order
x_cuts = (x > min_x) & (x < max_x) #Set cuts in x to trim the left and right edges of the detector


#Iterate between flexure and airmass (power) corrections, narrowing the range for each iterations
best_fit_flexure = input_flexure #Define the default best fit flexure and power which will be iterated upon
best_fit_power = input_power

if automate_fit:
	all_iterations_flexure = []
	all_iterations_chisq_flexure = []
	all_iterations_power = []
	all_iterations_chisq_power = []
	for iteration in range(n_iterations):  #Loop through the iterations
		print('Iteration '+str(iteration+1)+' of '+str(n_iterations))
		print('Best fit flexure from last iteration = ', best_fit_flexure)
		print('Best fit power from last iteration= ', best_fit_power)
		print('Flexure range to search = '+str(flexure_array[0])+' to '+str(flexure_array[-1]))
		print('Power range to search = '+str(power_array[0])+' to '+str(power_array[-1]))
		#Calculate flexure difference using cross-correlation
		chisq_flexure = np.zeros(n_flexure)	 #Define array to store the chisq for the flexure
		for order in range(skip_first_orders, n_orders-skip_last_orders): #Loop through all orders in specified order range
			sci_order = data_sci[order,:] #Grab subset of science and A0V data for this order
			a0v_flattened_order = data_a0v_flattened[order,:]
			smoothed_sci_order = normalize_continuum(sci_order, size=75) #Contiunuum normalize the science data to place the telluric lines in it on the same flux scale as the flattened A0V data
			cuts = np.isfinite(a0v_flattened_order) & np.isfinite(sci_order) &  binary_dilation((a0v_flattened_order > min_a0v) & (a0v_flattened_order < max_a0v)) & x_cuts #Set cuts for which pixels to use
			if np.any(cuts == True): #Error catch to ignore orders with no useful telluric lines
				for i in range(n_flexure): #Loop to check all flexure values for this order and iteration
					chisq_flexure[i] += bn.nansum((srebin(x, x+flexure_array[i], a0v_flattened_order)[cuts]**best_fit_power - smoothed_sci_order[cuts])**2) #Calculate chisq for flexure using the difference between the A0V and a continuum normalized science spectrum
		best_fit_flexure = flexure_array[chisq_flexure == bn.nanmin(chisq_flexure)][0]
		#Calculate airmass (powerlaw) difference using cross-correlation
		chisq_power = np.zeros(n_power) #Define array to store chisq for the power
		for order in range(skip_first_orders, n_orders-skip_last_orders): #Loop through all orders in specified order range
			sci_order = data_sci[order,:] #Grab subset of science and A0V data for this order
			a0v_flattened_order = data_a0v_flattened[order,:]
			for i in range(n_power): #Loop to check all power values for this order and iteration
				ratio = sci_order / srebin(x, x+best_fit_flexure, a0v_flattened_order)**power_array[i] #Divide science data by A0V data for a given power
				smoothed_ratio = robust_median_filter(ratio, size=75) #Smooth the ratio
				cuts = np.isfinite(ratio) & np.isfinite(smoothed_ratio) &  binary_dilation((a0v_flattened_order > min_a0v) &  (a0v_flattened_order < max_a0v)) & x_cuts #Set cuts for which pixels to use
				if np.any(cuts == True): #Error catch
					chisq_power[i] += bn.nansum((ratio[cuts] - smoothed_ratio[cuts])**2) #Calculate chicsq for the power by comparing the unsmoothed ratio of the science data to the power/airmass corrected A0V data to a smoothed version of said ratio, a good correction would not vary much from the smoothed correction
		best_fit_power = power_array[chisq_power == bn.nanmin(chisq_power)][0]
		#Store chisq values for all iterations
		all_iterations_flexure.append(flexure_array)
		all_iterations_chisq_flexure.append(chisq_flexure)
		all_iterations_power.append(power_array)
		all_iterations_chisq_power.append(chisq_power)
		#Recalculate parameter ranges for next iteration
		n_flexure = len(flexure_array)
		n_power = len(power_array)
		delta_flexure = (flexure_array[1] - flexure_array[0]) / 4.0
		delta_power = (power_array[1] - power_array[0]) / 4.0
		min_flexure = best_fit_flexure - (delta_flexure*(n_flexure)/2.0)
		max_flexure = best_fit_flexure + (delta_flexure*(n_flexure)/2.0)
		min_power = best_fit_power - (delta_power*(n_power)/2.0)
		max_power = best_fit_power + (delta_power*(n_power)/2.0)
		flexure_array = np.arange(min_flexure, max_flexure, delta_flexure)
		power_array = np.arange(min_power, max_power, delta_power)
		n_flexure = len(flexure_array)
		n_power = len(power_array)


#Correct data and save back to fits files
corrected_data_a0v = np.zeros(np.shape(data_sci))
corrected_data_a0v_flattened = np.zeros(np.shape(data_sci))
for order in range(n_orders): #Apply best fit corrections to A0V data
	corrected_data_a0v[order,:] = srebin(x, x+best_fit_flexure, data_a0v[order,:])**best_fit_power
	corrected_data_a0v_flattened[order,:] = srebin(x, x+best_fit_flexure, data_a0v_flattened[order,:])**best_fit_power
hdul_a0v[0].data = corrected_data_a0v #Put the corrected data back into the fits files and save the otuputted fits file as *.corrected.fits
hdul_flattened_a0v[0].data = corrected_data_a0v_flattened
hdul_a0v.writeto(a0v_fits_filename[:-5]+'.corrected.fits', overwrite=True)
hdul_flattened_a0v.writeto(flattened_a0v_fits_filename[:-5]+'.corrected.fits', overwrite=True)
hdul_a0v.close()
hdul_flattened_a0v.close()

#Plot telluric corrected science target before and after the found best fit flexure and airmass corrections for the A0V standard
with PdfPages(a0v_fits_filename[:-5]+'.corrected.pdf') as pdf:
	if automate_fit:
		#Plot chisq results for flexure and power (airmass) corrections
		pyplot.figure() #Plot flexure chisq
		for i in range(n_iterations):
			pyplot.plot(all_iterations_flexure[i], all_iterations_chisq_flexure[i], label='Iteration '+str(i+1))
		pyplot.legend()
		pyplot.xlabel('Flexure (x pixels)')
		pyplot.ylabel('Chisq')
		pdf.savefig()
		pyplot.close()
		pyplot.figure() #Plot power chisq
		for i in range(n_iterations):
			pyplot.plot(all_iterations_power[i], all_iterations_chisq_power[i], label='Iteration '+str(i+1))
		pyplot.legend()
		pyplot.xlabel('Power')
		pyplot.ylabel('Chisq')
		pdf.savefig()
		pyplot.close()
	#Plot uncorrected and corrected science data for comparison
	uncorrected_color = 'red'
	corrected_color = 'blue'
	for ordernum in range(n_orders):
		a0v_order = data_a0v_flattened[ordernum,:]
		x = np.arange(len(a0v_order))
		a0v_order_corrected = srebin(x, x+best_fit_flexure, a0v_order**best_fit_power,)
		sci_order = data_sci[ordernum,:]
		ratio = sci_order/a0v_order
		ratio_corrected = sci_order/a0v_order_corrected
		masked_ratio = copy.deepcopy(ratio)
		masked_ratio_corrected = copy.deepcopy(ratio_corrected)
		if (ordernum >= skip_first_orders) and (ordernum < n_orders - skip_last_orders):
			mask = binary_dilation((a0v_order > min_a0v) & (a0v_order < max_a0v)) & x_cuts
		else:
			mask = a0v_order == 100.0 #mask everything
		masked_ratio[~mask] = np.nan
		masked_ratio_corrected[~mask] = np.nan
		pyplot.clf()
		fig, axs = pyplot.subplots(4)
		fig.fontsize=8
		max_x = np.max(x)
		xi = (x >= 0) & (x < max_x/4)
		axs[0].tick_params(labelsize=6)
		axs[0].plot(x[xi], ratio[xi], linestyle='dotted', linewidth=0.5, color=uncorrected_color)
		axs[0].plot(x[xi], ratio_corrected[xi], linestyle='dotted', linewidth=0.5, color=corrected_color)
		axs[0].plot(x[xi], masked_ratio[xi], linewidth=0.5, color=uncorrected_color)
		axs[0].plot(x[xi], masked_ratio_corrected[xi], linewidth=0.5, color=corrected_color)
		xi = (x >= max_x/4) & (x < 2*max_x/4)
		axs[1].tick_params(labelsize=6)
		axs[1].plot(x[xi], ratio[xi], linestyle='dotted', linewidth=0.5, color=uncorrected_color)
		axs[1].plot(x[xi], ratio_corrected[xi], linestyle='dotted', linewidth=0.5, color=corrected_color)
		axs[1].plot(x[xi], masked_ratio[xi], linewidth=0.5, color=uncorrected_color)
		axs[1].plot(x[xi], masked_ratio_corrected[xi], linewidth=0.5, color=corrected_color)
		xi = (x >= 2*max_x/4) & (x < 3*max_x/4)
		axs[2].tick_params(labelsize=6)
		axs[2].plot(x[xi], ratio[xi],  linestyle='dotted', linewidth=0.5, color=uncorrected_color)
		axs[2].plot(x[xi], ratio_corrected[xi], linestyle='dotted', linewidth=0.5, color=corrected_color)
		axs[2].plot(x[xi], masked_ratio[xi], linewidth=0.5, color=uncorrected_color)
		axs[2].plot(x[xi], masked_ratio_corrected[xi], linewidth=0.5, color=corrected_color)
		pyplot.suptitle('Order = '+str(ordernum))
		xi = (x >= 3*max_x/4) & (x < 4*max_x/4)
		axs[3].tick_params(labelsize=6)
		axs[3].plot(x[xi], ratio[xi], linestyle='dotted', linewidth=0.5, color=uncorrected_color)
		axs[3].plot(x[xi], ratio_corrected[xi],  linestyle='dotted', linewidth=0.5, color=corrected_color)
		axs[3].plot(x[xi], masked_ratio[xi], label='Before', linewidth=0.5, color=uncorrected_color)
		axs[3].plot(x[xi], masked_ratio_corrected[xi], label='After', linewidth=0.5, color=corrected_color)
		pyplot.suptitle('Order = ' + str(ordernum) + '  Flexure = %7.4f' % best_fit_flexure + '   Power = %7.4f' % best_fit_power)
		pyplot.legend()
		maxy = bn.nanmedian(sci_order)*3.0
		miny = -0.1 * maxy
		pyplot.ylim([miny, maxy])
		pdf.savefig()
		pyplot.close()

print('best_fit_flexure = ', best_fit_flexure)
print('best_fit_powerlaw = ', best_fit_power)
