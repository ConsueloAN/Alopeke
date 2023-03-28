from astropy.io import fits
from astropy.table import Table, vstack
import sep
import numpy as np
import UTILS as ut

date = '20220824'
n = ["04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27"]
camera = "r"
tipo = "reduced"     #change filename

x_FRB_r = [129.3]
y_FRB_r = [124.7]
x_FRB_b = [119.6]
y_FRB_b = [128.2]

radio = 2
FWHM = 2*radio

sizebox = 20

# Centres are estimated in star_centres file in gdrive:

# star1 is the brighter
x_star1_r = 70.5
y_star1_r = 201.0
x_star2_r = 78.9
y_star2_r = 241.0
x_star3_r = 219.7
y_star3_r = 241.3
x_bkg_r = [66.5]
y_bkg_r = [124.0]

x_star1_b = 178.0
y_star1_b = 204.5
x_star2_b = 170.1
y_star2_b = 244.6
x_star3_b = 29.7
y_star3_b = 248.0
x_bkg_b = [182.0]
y_bkg_b = [127.5]


for n in n:
	filename = '/media/consuelo/M2/projects/alopeke/data_reduced/'+ date +'/N'+ date +'A00' + n + camera + '_reduced.fits'
	fitsfile = fits.open(filename)
	#UTC = fitsfile[0].header['UTC']
	data = fitsfile[0].data
	data = data.astype(float)

	ind = []
	flux_1FWHM = []
	flux_2FWHM = []
	fluxerr_1FWHM = []
	fluxerr_2FWHM = []

	flux_star1_1FWHM_recov = []
	flux_star1_2FWHM_recov = []
	x_star1_recov = []
	y_star1_recov = []
	SN_star1_recov = []

	flux_star2_1FWHM_recov = []
	flux_star2_2FWHM_recov = []
	x_star2_recov = []
	y_star2_recov = []
	SN_star2_recov = []

	flux_star3_1FWHM_recov = []
	flux_star3_2FWHM_recov = []
	x_star3_recov = []
	y_star3_recov = []
	SN_star3_recov = []	

	flux_bkg_1FWHM_recov = []
	flux_bkg_2FWHM_recov = []
	x_bkg_recov = []
	y_bkg_recov = []
	
	for i in range(4999):
		ind += [i+1]
		image = data[i]
		bkg = sep.Background(image)
		image_sub = image - bkg

		if camera == "r":

			flux_2FWHM_, fluxerr_2FWHM_, flag_2FWHM_ = sep.sum_circle(image_sub, x_FRB_r, y_FRB_r, 2*FWHM, err=bkg.globalrms, gain=1.0)		
			flux_1FWHM_, fluxerr_1FWHM_, flag_1FWHM_ = sep.sum_circle(image_sub, x_FRB_r, y_FRB_r, FWHM, err=bkg.globalrms, gain=1.0)
		
			flux_1FWHM += [flux_1FWHM_[0]]
			flux_2FWHM += [flux_2FWHM_[0]]
			fluxerr_1FWHM += [fluxerr_1FWHM_[0]]
			fluxerr_2FWHM += [fluxerr_2FWHM_[0]]

			mask1 = ut.create_mask(x_star1_r, y_star1_r, sizebox)
			flux_star1_1FWHM_recov_, x_star1_recov_, y_star1_recov_, SN_star1_recov_ = ut.bright_star_parameters(image_sub, mask1, FWHM)
			flux_star1_2FWHM_recov_, x_star1_recov_, y_star1_recov_, SN_star1_recov_ = ut.bright_star_parameters(image_sub, mask1, 2*FWHM)

			flux_star1_1FWHM_recov += [flux_star1_1FWHM_recov_[0]]
			flux_star1_2FWHM_recov += [flux_star1_2FWHM_recov_[0]]
			x_star1_recov += [x_star1_recov_[0]]
			y_star1_recov += [y_star1_recov_[0]]
			SN_star1_recov += [SN_star1_recov_[0]]

			mask2 = ut.create_mask(x_star2_r, y_star2_r, sizebox)
			flux_star2_1FWHM_recov_, x_star2_recov_, y_star2_recov_, SN_star2_recov_ = ut.bright_star_parameters(image_sub, mask2, FWHM)
			flux_star2_2FWHM_recov_, x_star2_recov_, y_star2_recov_, SN_star2_recov_ = ut.bright_star_parameters(image_sub, mask2, 2*FWHM)

			flux_star2_1FWHM_recov += [flux_star2_1FWHM_recov_[0]]
			flux_star2_2FWHM_recov += [flux_star2_2FWHM_recov_[0]]
			x_star2_recov += [x_star2_recov_[0]]
			y_star2_recov += [y_star2_recov_[0]]
			SN_star2_recov += [SN_star2_recov_[0]]

			mask3 = ut.create_mask(x_star3_r, y_star3_r, sizebox)
			flux_star3_1FWHM_recov_, x_star3_recov_, y_star3_recov_, SN_star3_recov_ = ut.bright_star_parameters(image_sub, mask3, FWHM)
			flux_star3_2FWHM_recov_, x_star3_recov_, y_star3_recov_, SN_star3_recov_ = ut.bright_star_parameters(image_sub, mask3, 2*FWHM)

			flux_star3_1FWHM_recov += [flux_star3_1FWHM_recov_[0]]
			flux_star3_2FWHM_recov += [flux_star3_2FWHM_recov_[0]]			
			x_star3_recov += [x_star3_recov_[0]]
			y_star3_recov += [y_star3_recov_[0]]
			SN_star3_recov += [SN_star3_recov_[0]]

			flux_bkg_2FWHM_, fluxerr_bkg_2FWHM_, flag_bkg_2FWHM_ = sep.sum_circle(image_sub, x_bkg_r, y_bkg_r, 2*FWHM, err=bkg.globalrms, gain=1.0)		
			flux_bkg_1FWHM_, fluxerr_bkg_1FWHM_, flag_bkg_1FWHM_ = sep.sum_circle(image_sub, x_bkg_r, y_bkg_r, FWHM, err=bkg.globalrms, gain=1.0)
		
			flux_bkg_1FWHM_recov += [flux_bkg_1FWHM_[0]]
			flux_bkg_2FWHM_recov += [flux_bkg_2FWHM_[0]]

		else:

			flux_2FWHM_, fluxerr_2FWHM_, flag_2FWHM_ = sep.sum_circle(image_sub, x_FRB_b, y_FRB_b, 2*FWHM, err=bkg.globalrms, gain=1.0)		
			flux_1FWHM_, fluxerr_1FWHM_, flag_1FWHM_ = sep.sum_circle(image_sub, x_FRB_b, y_FRB_b, FWHM, err=bkg.globalrms, gain=1.0)
		
			flux_1FWHM += [flux_1FWHM_[0]]
			flux_2FWHM += [flux_2FWHM_[0]]
			fluxerr_1FWHM += [fluxerr_1FWHM_[0]]
			fluxerr_2FWHM += [fluxerr_2FWHM_[0]]

			mask1 = ut.create_mask(x_star1_b, y_star1_b, sizebox)
			flux_star1_1FWHM_recov_, x_star1_recov_, y_star1_recov_, SN_star1_recov_ = ut.bright_star_parameters(image_sub, mask1, FWHM)
			flux_star1_2FWHM_recov_, x_star1_recov_, y_star1_recov_, SN_star1_recov_ = ut.bright_star_parameters(image_sub, mask1, 2*FWHM)

			flux_star1_1FWHM_recov += [flux_star1_1FWHM_recov_[0]]
			flux_star1_2FWHM_recov += [flux_star1_2FWHM_recov_[0]]
			x_star1_recov += [x_star1_recov_[0]]
			y_star1_recov += [y_star1_recov_[0]]
			SN_star1_recov += [SN_star1_recov_[0]]	

			mask2 = ut.create_mask(x_star2_b, y_star2_b, sizebox)
			flux_star2_1FWHM_recov_, x_star2_recov_, y_star2_recov_, SN_star2_recov_ = ut.bright_star_parameters(image_sub, mask2, FWHM)
			flux_star2_2FWHM_recov_, x_star2_recov_, y_star2_recov_, SN_star2_recov_ = ut.bright_star_parameters(image_sub, mask2, 2*FWHM)

			flux_star2_1FWHM_recov += [flux_star2_1FWHM_recov_[0]]
			flux_star2_2FWHM_recov += [flux_star2_2FWHM_recov_[0]]
			x_star2_recov += [x_star2_recov_[0]]
			y_star2_recov += [y_star2_recov_[0]]
			SN_star2_recov += [SN_star2_recov_[0]]	

			mask3 = ut.create_mask(x_star3_b, y_star3_b, sizebox)
			flux_star3_1FWHM_recov_, x_star3_recov_, y_star3_recov_, SN_star3_recov_ = ut.bright_star_parameters(image_sub, mask3, FWHM)
			flux_star3_2FWHM_recov_, x_star3_recov_, y_star3_recov_, SN_star3_recov_ = ut.bright_star_parameters(image_sub, mask3, 2*FWHM)

			flux_star3_1FWHM_recov += [flux_star3_1FWHM_recov_[0]]
			flux_star3_2FWHM_recov += [flux_star3_2FWHM_recov_[0]]
			x_star3_recov += [x_star3_recov_[0]]
			y_star3_recov += [y_star3_recov_[0]]
			SN_star3_recov += [SN_star3_recov_[0]]

			flux_bkg_2FWHM_, fluxerr_bkg_2FWHM_, flag_bkg_2FWHM_ = sep.sum_circle(image_sub, x_bkg_b, y_bkg_b, 2*FWHM, err=bkg.globalrms, gain=1.0)		
			flux_bkg_1FWHM_, fluxerr_bkg_1FWHM_, flag_bkg_1FWHM_ = sep.sum_circle(image_sub, x_bkg_b, y_bkg_b, FWHM, err=bkg.globalrms, gain=1.0)
		
			flux_bkg_1FWHM_recov += [flux_bkg_1FWHM_[0]]
			flux_bkg_2FWHM_recov += [flux_bkg_2FWHM_[0]]

	tab = Table()
	tab['frame'] = ind
	#tab['UTC'] = [UTC]
	tab['flux_1FWHM'] = flux_1FWHM
	tab['fluxerr_1FWHM'] = fluxerr_1FWHM
	tab['flux_2FWHM'] = flux_2FWHM
	tab['fluxerr_2FWHM'] = fluxerr_2FWHM

	tab['flux_star1_1FWHM'] = flux_star1_1FWHM_recov
	tab['flux_star1_2FWHM'] = flux_star1_2FWHM_recov
	tab['x_star1'] = x_star1_recov
	tab['y_star1'] = y_star1_recov
	tab['SN_star1'] = SN_star1_recov

	tab['flux_star2_1FWHM'] = flux_star2_1FWHM_recov
	tab['flux_star2_2FWHM'] = flux_star2_2FWHM_recov
	tab['x_star2'] = x_star2_recov
	tab['y_star2'] = y_star2_recov
	tab['SN_star2'] = SN_star2_recov

	tab['flux_star3_1FWHM'] = flux_star3_1FWHM_recov
	tab['flux_star3_2FWHM'] = flux_star3_2FWHM_recov
	tab['x_star3'] = x_star3_recov
	tab['y_star3'] = y_star3_recov
	tab['SN_star3'] = SN_star3_recov

	tab['flux_bkg_1FWHM'] = flux_bkg_1FWHM_recov
	tab['flux_bkg_2FWHM'] = flux_bkg_2FWHM_recov


	#print(tab)
	print(n)

	tab.write('/media/consuelo/M2/projects/alopeke/count_tables/'+ date +'/count_table_' + tipo + '_' + n + '_' + camera + '.fits', overwrite=True)
