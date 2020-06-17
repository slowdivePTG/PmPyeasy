#! /usr/bin/python

#aperture photometry with iraf/apphot

import numpy as np
from pyraf import iraf
from iraf import digiphot,apphot
import os
from astropy.io import fits, ascii
from astropy.table import Table

    
def apphot_phot_iraf(image, centroid_algorithm='centroid', coo_list='default', ap_list='default',coo_list_format='daophot', fwhmpsf=3, sigma=10, app=8, skyin=15, wsky=15, sky_sigma_down=3, sky_sigma_up=3, readnoise=5, epadu=1, datamin = 'INDEF', datamax = 'INDEF', emission='yes', zeropt=25, single_point=True, outkeys=None):
    '''    
    INPUTS:
	image:
	centroid_algorithm: 'none', 'centroid', 'gauss' or 'ofilter' 
	coo_list: the input coordinate list(s); default: image.coo.?
	ap_list: output photometry file(s); default: image.mag.?
	fwhmpsf: FWHM of the PSF in scale unit
	sigma: standard deiation of background in counts
	app:  aperture radius(radii) in scale unit
	skyin: inner sky annulus radius
	wsky: the width of sky annulus 
	sky_sigma_down: lower k-sigma rejection limit in sky sigma
	sky_sigma_up: upper k-sigma rejection limit in sky sigma
	readnoise: CCD readout noise in electrons
	epadu: gain in electron per count
	datamin: minimum good data value
	datamax: maximum good data value

    OUTPUTS:
	The full photometry will be saved in aperture photometry file ap_list, and for each star measured the following record is written :

	image  xinit  yinit  id  coords  lid
	xcenter  ycenter  xshift  yshift  xerr  yerr  cier cerror
	msky  stdev  sskew  nsky  nsrej  sier  serror
	itime  xairmass  ifilter  otime
	rapert  sum  area  mag  merr  pier  perr

	Image and coords are the name of the image and coordinate file respectively.
	Id and lid are the sequence numbers of stars in the output and coordinate files respectively.
	Sier and serror are the sky fitting error code and accompanying error message respectively.
	Msky, stdev and sskew are the best estimate of the sky value (per pixel), standard deviation and skew respectively.
	Nsky and nsrej are the number of sky pixels and the number of sky pixels rejected respectively.
	Itime is the exposure time, xairmass is self-evident, ifilter is an id string identifying the filter used in the observations, and otime is a string containing the time of the observation in whatever units the user has set up.
	Rapert, sum, area, and flux are the radius of the aperture in scale units, the total number of counts including sky in the aperture, the area of the aperture in square pixels, and the total number of counts excluding sky in the aperture.
	Mag and merr are the magnitude and error in the magnitude in the aperture (see below).
	
	        flux = sum - area * msky
       		mag = zmag - 2.5 * log10 (flux) + 2.5 * log10 (itime)
        	merr = 1.0857 * err / flux
         	err = sqrt (flux / epadu + area * stdev**2 + area**2 * stdev**2 / nsky)
    '''

    apphot.phot.unlearn()

    apphot.datapars.unlearn()
    apphot.datapars.scale = 1.0
    apphot.datapars.fwhmpsf = fwhmpsf
    apphot.datapars.emission = emission
    apphot.datapars.sigma = sigma

    apphot.datapars.datamin = datamin
    apphot.datapars.datamax = datamax

    apphot.datapars.noise = 'poisson'

    apphot.datapars.ccdread = ''
    apphot.datapars.readnoise = readnoise
    apphot.datapars.itime = 1.0
    apphot.datapars.epadu = epadu
    apphot.datapars.xairmass = 'INDEF'
    apphot.datapars.ifilter = 'INDEF'
    apphot.datapars.otime = 'INDEF'
    
    # iraf.digiphot.apphot.centerpars : 
    apphot.centerpars.calgorithm = centroid_algorithm
    apphot.centerpars.cbox = 10.0
    apphot.centerpars.cthreshold = 0.0
    apphot.centerpars.minsnratio = 1.0
    apphot.centerpars.cmaxiter = 10.0
    apphot.centerpars.maxshift = 2.0
    apphot.centerpars.clean = False
    apphot.centerpars.rclean = 1.0
    apphot.centerpars.rclip = 2.0
    apphot.centerpars.kclean = 3.0
    apphot.centerpars.mkcenter = False

    # iraf.digiphot.apphot.fitskypars : 
    apphot.fitskypars.unlearn()
    apphot.fitskypars.salgorithm = 'median'
    apphot.fitskypars.annulus = skyin
    apphot.fitskypars.dannulus = wsky
    apphot.fitskypars.skyvalue = 0.0
    apphot.fitskypars.smaxiter = 10.0
    apphot.fitskypars.sloclip = 0.0
    apphot.fitskypars.shiclip = 0.0
    apphot.fitskypars.snreject = 50.0
    apphot.fitskypars.sloreject = sky_sigma_down
    apphot.fitskypars.shireject = sky_sigma_up
    apphot.fitskypars.khist = 3.0
    apphot.fitskypars.binsize = 0.1
    apphot.fitskypars.smooth = False
    apphot.fitskypars.rgrow = 0.0
    apphot.fitskypars.mksky = False

    # iraf.digiphot.apphot.photpars : 
    apphot.photpars.unlearn()
    apphot.photpars.weighting = 'constant'
    apphot.photpars.apertures = app
    apphot.photpars.zmag = zeropt
    apphot.photpars.mkapert = False
          
    photparams = {
        'radplot':False,
        }
       
    if os.path.exists(ap_list): 
        os.remove(ap_list) 

    # run photometry using the newly created coxyfile for providing input coordinates
    apphot.phot(image=image, skyfile='', output=ap_list, coords=coo_list, verify='no',interactive='no',verbose=True, Stdout=1, **photparams)
    
    #get the photometric data from the apphot output    
    if outkeys is not None:
	if ap_list == 'default':
		for i in np.arange(100)[::-1]:
			ap_list = image+'.mag.'+str(i)
			if os.path.exists(ap_list):
				break
	photret = ascii.read(ap_list, format='daophot')
	seletedphot = photret[outkeys]
	for colname in seletedphot.keys():
		seletedphot[colname]=seletedphot[colname].filled(fill_value=99.99)
    else:	
	seletedphot = None
	
    return seletedphot


if __name__ == "__main__":	
	import optparse
	parser = optparse.OptionParser()

	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_x = None
	parser.add_option('-x', dest = 'x', type = float, default= def_x, help='x position of target; default: %s'%def_x )

	def_y = None
	parser.add_option('-y', dest = 'y', type = float, default= def_x, help='y position of target; default: %s'%def_y )

	def_coor = 'default'
	parser.add_option('--coor', dest='coor', type='string', default=def_coor, help='input coordinate list(s); default: %s.coo.?'%def_inputimage)

	def_apphotoutput = 'default'
	parser.add_option('--apphot_output', dest='apphot_output', type='string', default=def_apphotoutput, help='iraf/apphot internal output file; default: %s.mag.?'%def_inputimage)

	def_centroid = 'centroid'
	parser.add_option('--centroid', dest='centroid', type='string', default=def_centroid, help='centroid algorithm, valid inputs are centroid, none, gauss and ofilter; default: %s'%def_centroid)

	def_app = 8
	parser.add_option('-a','--aperture_radius', dest='app', type=float, default=def_app, help='aperture radius(radii) in scale unit; default: %s'%def_app)


	def_skyin = 15
	parser.add_option('--skyin', dest='skyin', type = float, default=def_skyin, help='inner sky annulus radius; default: %s'%def_skyin)


	def_wsky = 15
	parser.add_option('--sky_width', dest='sky_width', type = float, default=def_wsky, help='the width of sky annulus; default: %s'%def_wsky)


	def_fwhmpsf = 5
	parser.add_option('--fwhmpsf', dest = 'fwhmpsf', type = float, default= def_fwhmpsf, help='FWHM of the PSF in scale unit; default: %s'%def_fwhmpsf )


	def_sigma = 10
	parser.add_option('--sigma', dest = 'sigma', type = float , default= def_sigma, help='standard of deviation of background in counts; default: %s'%def_sigma )
	
	def_datamin = 'INDEF'
	parser.add_option('--datamin', dest = 'datamin', type = "string", default= def_datamin, help='minimum good value; default: %s'%def_datamin )

	def_datamax = 'INDEF'
	parser.add_option('--datamax', dest = 'datamax', type = "string", default= def_datamax, help='maximum good value; default: %s'%def_datamax )


	def_readnoise = 5
	parser.add_option('--readnoise', dest = 'readnoise', type = float, default= def_readnoise, help='CCD readnoise in electrons; default: %s'%def_readnoise )

	def_epadu = 1
	parser.add_option('--epadu', dest = 'epadu', type = float, default= def_epadu, help='gain in electrons per count; default: %s'%def_epadu )

	def_okeys = 'XCENTER,YCENTER,MSKY,STDEV,NSKY,AREA,SUM,MAG,MERR'
	parser.add_option('--colnames', dest='record_names', type='string', default=def_okeys, help='the names of the record for which to save to photometry table; default:%s'%def_okeys)

	def_output = ''
	parser.add_option('-o','--output', dest='output', type='string', default=def_output, help='the filename of the output file where the selected record of the photometry will be saved; default: %s'%def_output)

	def_daophot_format= False
	parser.add_option('-d', '--daophot', dest='daophot', action="store_true", default=def_daophot_format, help="whether the input of the source list are in default daophot format")

	options, remainder = parser.parse_args()

	input_image = options.input_image
	centroid_algorithm = options.centroid
	coo_list = options.coor
	coo_list_format=options.daophot
	if coo_list_format:
		coo_list_format = 'daophot'
	else:
		coo_list_format = 'simptxt'

	x = options.x
	y = options.y

	ap_list = options.apphot_output
	fwhmpsf= options.fwhmpsf
	sigma= options.sigma

	app = options.app
	skyin = options.skyin
	wsky = options.sky_width

	outkeys = options.record_names
	outkeys = outkeys.split(',')
	if len(outkeys)==0:
		outkeys = None
	outfile = options.output


	datamin= options.datamin
	if datamin != 'INDEF':
		datamin = float(datamin)

	datamax= options.datamax
	if datamax != 'INDEF':
		datamax = float(datamax)

	readnoise= options.readnoise
	epadu= options.epadu

	emission = 'yes'
	sky_sigma_down = 3
	sky_sigma_up = 3
	zeropt= 25

	if not os.path.exists(input_image):
		raise ValueError("input image %s doesn't exist"%input_image)

	single_point = False
	if x is not None and y is not None:
		coo_list = 'temp.coo'
		fid = open(coo_list,'w')
		fid.write(str(x)+' ')
		fid.write(str(y))
		fid.close()
		single_point = True
	elif coo_list == 'default':
		print "the input coordinate list: %s.coo.?"%input_image
	else:
		print "the input coordinate list: %s"%coo_list


	rettable = apphot_phot_iraf(input_image, centroid_algorithm=centroid_algorithm, coo_list=coo_list, ap_list=ap_list, coo_list_format=coo_list_format, fwhmpsf=fwhmpsf, sigma=sigma, app=app, skyin=skyin, wsky=wsky, sky_sigma_down=sky_sigma_down, sky_sigma_up=sky_sigma_up, readnoise=readnoise, epadu=epadu, datamin =datamin, datamax = datamax, single_point = single_point, outkeys=outkeys)

	print rettable

	if outfile != '':
		rettable.write(outfile, format='ascii.commented_header')
