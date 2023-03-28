import fitsio
import numpy as np
import os
from astropy.io import fits
from scipy import stats as st

n = ["04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26"]
camera = "r"
date = "20220908"

for n in n:

	filename = '/media/consuelo/M2/projects/alopeke/' + date + '/N' + date + 'A00' + n + camera + '.fits'
	outfile = '/media/consuelo/M2/projects/alopeke/data_reduced/'+ date +'/N' + date + 'A00' + n + camera + '_reduced.fits'

	fitsfile = fits.open(filename)
	images = fitsfile[0].data
	images = images.astype(float)

	master_bias = fitsio.read('/media/consuelo/M2/projects/alopeke/reductions/masterbias_' + camera + '_cropped.fits')     #same for two dates
	master_flat = fitsio.read('/media/consuelo/M2/projects/alopeke/reductions/masterflat_1' + camera + '_cropped.fits')      #idem

	for i in range(5000):
		images[i] = (images[i] - master_bias) / master_flat

	images = np.delete(images,0,0)                 # cube without the first image
	hdu = fits.PrimaryHDU(images)
	hdulist = fits.HDUList([hdu])
	hdulist.writeto(outfile, overwrite=True)

	print(n)
