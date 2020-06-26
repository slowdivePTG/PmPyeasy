#!/usr/bin/python

#pipeline for photometry by Ping Chen

#modified by Ping on 2016-06-19: make aperture and PSF photometry compossible
#philosophy: the default is aperture photometry; PSF photometry needs specific flag

#modified by Ping on 2016-07-04: remove things about triming images

#SN position (x,y) comes from fits image style, which treat the first pixel as (1,1)
#Many data process (sfind, dophot, iraf apphot etc) treat the first poiont as (0,0)

#modification on 2016-10-11: (for aperture photometry,) aperture size used to be fixed to a value for all images;
#now adjust aperture size as well as background according to the image fwhm: aperture size 2*fwhm, sky_in 3*fwhm, sky_out 5*fwhm

#modification on 2016-10-28: aperture size 1.5*fwhm
#modify matplotlib related part under crontab job restriction(GUI not supported)
#self.standards type changed from ndarray to astropy.table

#modification on 2018-02-26: add functions to integrate photometry for WFCAM images
#add new image information element airmass to the photometry record table

#todo: estimate image background and improve get_fwhm_fistar
#todo: use astropy.stats.sigma_clipped_stats

#from __future__ import print_function

import numpy as np
import os
import shutil
import wget
from astropy.table import Table,Column, vstack, hstack
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from scipy.optimize import curve_fit
from collections import OrderedDict
import ConfigParser

try:
	import pyds9
except:
	print "pyds9 import failure..."

try:
	import matplotlib.pylab as plt
except:
	print "no display mode"

from extern import cosmics
from utils.ccdproc_cosmicray_removal_lacosmic import remove_CR_ccdproc_cosmicray_lacosmic
from utils.common import __get_fits_info as get_fits_info
from utils.common import sigma_clipping,get_best, radec_format_transformation

from utils.prepare_work import setting_dir
from utils.get_image_info_difftelescope import get_obstime,get_exptime,get_filter,get_bitpix, get_airmass_telescope
from utils.timeout_handle import Command
from utils.catalog_search import query_local_APASS,query_VO_SCS,query_General_MAST,panstarrs_query_Vizier,apass_query_Vizier,twomass_query_Vizier
from utils.photometry_collection import apphot_iraf

from utils.mag_flux_convertion import Vega_AB_mag_convertion, mag2flux, flux2mag
from utils.get_airmass_obs import get_airmass_given_time_target_observatory
from utils.Iowa_fits_header_simplify_standardize import simplify_and_standardize_fits_header

from mpltk.event_handler import mouse_pick

from iraftk.iraf_apphot_daofind import daofind_iraf
from iraftk.iraf_imstat import imstat_iraf
from iraftk.iraf_imarith import imarith_iraf
from iraftk.iraf_imcopy import imcopy_iraf
from iraftk.iraf_imcombine import imcombine_iraf
from iraftk.iraf_xytran import geomap_iraf, geoxytran_iraf

from pyds9tk.select_or_delete_from_given_list_ds9_display import get_xy_on_image_from_ds9, input_xys_through_ds9_get_mags
from pyds9tk.select_source_on_ds9 import get_lines_by_selecting_xys_on_image_from_ds9
from pyds9tk.ds9_display_sources import create_ds9_region_file, pylab_talk_to_ds9


import warnings
warnings.filterwarnings('ignore')

from photbasesetup import basesetup

import pickle

def load_phot_target(picklefile):
	'''
	Load pickled photometry object
	'''
	try:
		sn = pickle.load(open(picklefile,'r'))
	except:
		sn = None
		print "failed to load picklefile %s"%picklefile
	return sn

class photometry():
	'''
	photometry class
	'''
	def __init__(self,sn,telescope,pipe=True, auto_init = True, result_dir=None,data_repository_dir=None,photometry_method = 'apphot'):
		'''
		inputs:
			sn: the target name
			telescope: the telescope name
			pipe: whether in pipeline mode
			auto_init:
			result_dir: the master workplace folder for the given target
			data_repository_dir: the input data repository directory
			photometry_method:
		'''
		self.base = basesetup()
		self.pipe = pipe
		self.auto_init = auto_init
		self.current_sn = sn
		self.current_telescope = telescope
		self.photometry_method = photometry_method
		self.readme = "Please keep the notes below regarding the reduction: "

		self.ds9D = None #store the current instance of ds9
		#set up data repository and photometry directory system
		if not pipe:
			self.repository_dir_current_sn =  os.path.abspath(data_repository_dir)
                        self.current_sn_dir = os.path.abspath(result_dir)
		else:
			self.repository_dir_current_telescope = os.path.join(self.base.data_repository_master_dir,self.current_telescope)
			self.repository_dir_current_sn = os.path.join(self.repository_dir_current_telescope,self.current_sn)
			self.current_sn_master_dir = os.path.join(self.base.workplace_master_dir,self.current_telescope)
			self.current_sn_dir = os.path.join(self.current_sn_master_dir, self.current_sn)
		if not os.path.exists(self.repository_dir_current_sn):
			raise IOError("directory %s not found"%self.repository_dir_current_sn)
		if not os.path.exists(self.current_sn_dir):
			os.mkdir(self.current_sn_dir)

		self.__init_settings__()
		self.__init_dophot_pm()
		self.dophot_select_stars_method = 'eq' #
		self.dophot_select_star_value = 1
		self.__init_renew_trigger()
		self.__init_fitsh_grmatch_options()

		self.imagepreprocess_verbose = False
		self.cal_plot = False
		self.cal_offset_method='median'#'median','mean' or 'funcfit';the method adopted to get the offset between instrumental mags and std mags
		self.cal_offset_funcfit_type = 'constant' #'constant' or 'o1poly'

		self.cal_check_stds = False

		if self.photometry_method == 'apphot':
			self.result_table_file = os.path.join(self.result_dir,self.current_sn+'_photometry_info_AP.txt')
		elif self.photometry_method == 'psfphot':
			self.result_table_file = os.path.join(self.result_dir,self.current_sn+'_photometry_info_PSF.txt')
		else:
			raise IOError("Invalid input for photometry_method...")

		self.notesfile = os.path.join(self.result_dir, 'README.txt')

		self.new = True
		self.result_table = None

		if os.path.isfile(self.result_table_file):
			self.result_table =Table.read(self.result_table_file,format='ascii.fixed_width')
			self.new = False

		if not self.new:
			self.reln_image_to_datadir()

		self.ln_image_to_datadir(self.repository_dir_current_sn,self.raw_image_dir)

		self.images = self.get_images_dict(self.raw_image_dir)
		self.apphots  = {}
		self.psfphots = {}
		for img in self.images.keys():
			img_apphot  = os.path.join( self.aperture_photometry_dir, img.split('.')[0]+self.apphot_ret_file_suffix)
			img_psfphot = os.path.join( self.psf_photometry_dir, img.split('.')[0]+self.psfphot_ret_file_suffix)
			self.apphots[img] = img_apphot
			self.psfphots[img]= img_psfphot

		self.stdobs_match_files	= {}
		self.stdobs_transcoef_files = {}

		self.__init_photometry_info()
		self.__init_iraf_apphot_config()


		if self.new:
			self.__prepare_parameter_files()
			for key in self.photometry_info.keys():
				realimg_abs = os.path.realpath(self.images[key])
				self.photometry_info[key]['realimg'] = realimg_abs.split('/')[-1]

			if self.auto_init:
				for key in self.photometry_info.keys():
					self.__renew_flt_info(key)
				for key in self.photometry_info.keys():
					self.__renew_obstime_info(key)
				for key in self.photometry_info.keys():
					self.__renew_exptime_info(key)
				for key in self.photometry_info.keys():
					self.__renew_bitpix_info(key)

		else:
			self.__load_old_results(self.result_table,self.photometry_info)

			for key in self.photometry_info.keys():
				if key not in self.result_table['name']:
					realimg_abs = os.path.realpath(self.images[key])
					self.photometry_info[key]['realimg'] = realimg_abs.split('/')[-1]

			for key in self.photometry_info.keys():
				if key not in self.result_table['name']:
					self.__renew_flt_info(key)

			for key in self.photometry_info.keys():
				if key not in self.result_table['name']:
					self.__renew_obstime_info(key)

			for key in self.photometry_info.keys():
				if key not in self.result_table['name']:
					self.__renew_exptime_info(key)

			for key in self.photometry_info.keys():
				if key not in self.result_table['name']:
					self.__renew_bitpix_info(key)


	def __init_settings__(self):

		dirdict = setting_dir(self.current_sn_dir)
		self.raw_image_dir     = dirdict['data_dir']
		self.warehouse_dir     = dirdict['warehouse_dir']
		self.modified_image_dir= dirdict['modified_image_dir']
		self.stars_dir         = dirdict['stars_dir']
		self.template_dir      = dirdict['template_dir']
		self.result_dir        = dirdict['result_dir']
		self.log_dir           = dirdict['log_dir']
		self.std_ref_dir       = dirdict['std_ref_dir']
		self.subtraction_dir   = dirdict['subtraction_dir']
		self.parafile_dir      = dirdict['para_dir']
		self.psf_photometry_dir= dirdict['psfphot_dir']
		self.ds9_region_dir    = dirdict['ds9reg_dir']
		self.aperture_photometry_dir = dirdict['apphot_dir']

		self.parafile_template_dir = os.path.join(self.base.conffile_dir, 'para_templates')


		self.starfile_suffix = '.star'                        #This is the simplied version of source list with the format:
		self.iraf_daofind_output_suffix = '_iraf_daofind.out' #This is the output filename suffix for iraf/daofind, for example 000_iraf_daofind.out
		self.regfile_suffix = '.reg'
		self.apphot_ret_file_suffix = '.apphot'   	   #This is apphot photometry output modified/simplified version: x,y,mag,magerr
		self.phot_iraf_apphot_suffix = '_iraf_apphot.out'  #This is original apphot photometry output, refer to iraf for the detailed format
		self.apphot_iraf_parafile = 'apphot_iraf_param.par'
		self.psfphot_ret_file_suffix = '.psfphot'


		self.templates = {}
		self.templates_after_astrometry = {}
                self.template_factors = [0.7,0.3,-0.5]		#the weight in choosing the best image

		try:
			if self.current_sn in self.base.first_targets_info['Name']:
				RA_str  = self.base.first_targets_info[self.base.first_targets_info['Name']==self.current_sn]['RA'].data[0]
				Dec_str = self.base.first_targets_info[self.base.first_targets_info['Name']==self.current_sn]['Dec'].data[0]
			else:
				RA_str  = self.base.second_targets_info[self.base.second_targets_info['Name']==self.current_sn]['RA'].data[0]
				Dec_str = self.base.second_targets_info[self.base.second_targets_info['Name']==self.current_sn]['Dec'].data[0]

			RA_deg,Dec_deg = radec_format_transformation(RA_str,Dec_str)
			print "%s (RA,DEC)=(%s,%s): (%s,%s)"%(self.current_sn, RA_str, Dec_str, RA_deg, Dec_deg)

			self.sn_ra_world_deg = RA_deg
			self.sn_dec_world_deg = Dec_deg
			self.sn_ra_str = RA_str
			self.sn_dec_str = Dec_str
		except Exception as e:
			print e
			self.sn_ra_world_deg = None
			self.sn_dec_world_deg = None
			self.sn_ra_str = None
			self.sn_dec_str = None

		self.current_ccd_readnoise = None
		self.current_ccd_epadu = None

		self.pixscale = None
		self.std_region_radius = 0.33	#degree

		#local APASS data directories
		self.local_apass_dir = os.path.join(self.base.stdcat_dir,'APASS')

		self.std_ref_stars_file = {}
		self.std_ref_stars_file['apass'] = 'standard_reference_star_apass.txt'
		self.std_ref_stars_file['2mass'] = 'standard_reference_star_2mass.txt'
		self.std_ref_stars_file['panstarrs'] = 'standard_reference_star_panstarrs.txt'

		self.flt_std_ref_stars_file = {}
		self.flt_std_ref_stars_file_big = {}
		self.flt_std_ref_and_obs_measurement_stars_file = {}


		self.standards = None	#all information for standard stars

		self.panstarrs_stds_gp = None #only (x,y, mag_flt, magerr_flt)
		self.panstarrs_stds_rp = None
		self.panstarrs_stds_ip = None
		self.panstarrs_stds_zp = None
		self.panstarrs_stds_yp = None

		self.apass_stds_B = None
		self.apass_stds_V = None
		self.apass_stds_rp= None
		self.apass_stds_ip= None

		self.twomass_stds_J  = None
		self.twomass_stds_H  = None
		self.twomass_stds_K  = None

		self.reference_flt = None	# if template images don't have wcs info and we don't want API astrometry for all template then only one image with this 'reference_flt' will be given wcs and this wcs info will be used for other templates

		self.apass_catalog_method = 2
		self.apass_catalog_method_notes = "multiple methods to retrive APASS catalog provided: [1. local data (the whole APASS distribution) extraction; 2, apass_query_Vizier]"

		#The followings are related with stds from PanSTARRS catalog
		self.panstarrs_catalog_method = 2
		self.panstarrs_catalog_method_notes = "multiple methods to retrive PS1 catalog provided: [1, query_General_MAST (http://gsss.stsci.edu); 2, panstarrs_query_Vizier (based on astroquery.Vizier)]"

		self.twomass_catalog_method = 2
		self.twomass_catalog_method_notes = "multiple methods to retrive 2MASS catalog provided: [1, query_VO_SCS (https://irsa.ipac.caltech.edu); 2, twomass_query_Vizier (based on astroquery.Vizier)]"

		self.panstarrs_mag_photometry_method = 'PSF' # PanSTARRS catalog provides magnitudes from two photometry kinds: PSF or AP, only works for method 1
		self.panstarrs_min_nx = None		     # colname 'ng', 'nr', 'ni', 'nz', 'ny'
		self.panstarrs_min_nstackDetections = None

		self.selfphot_mag_faint_cut = None	#the magnitudes from private photometry results used as for calibration
		self.selfphot_mag_saturate_cut = None
		self.selfphot_magerr_up_cut = None #upper limit for magnitude uncertainty
		self.selfphot_match_method_for_RI_stds = 'surrounding_search' #if cal std from selfphot, when prepare std for R,I band then r,i std are needed, then matching

		self.panstarrs_mag_faint_cut = None
		self.panstarrs_mag_saturate_cut = None

		self.apass_mag_faint_cut = None
		self.apass_mag_saturate_cut = None
		self.apass_nobs_min = None
		self.apass_mobs_min = None
		self.apass_remove_single_measurement = True

		self.twomass_mag_faint_cut = None
		self.twomass_mag_saturate_cut = None

		self.physical_image_offset_x = 0	#difference between images coordinate and physical coordinate
		self.physical_image_offset_y = 0

		self.stdcal_xytran_stdnum = 3	#how many stars to select when xytran method used to link std and obs

		self.trim_before_grmatch = False
		self.trim_left = None		# For example, self.trim_left = n1, self.trim_right=n2, self.trim_bottom = n3, self.trim_up =n4
		self.trim_right = None		# the sources within [n1:n2] in x direction and [n3:n4] in y direction survive
		self.trim_bottom = None 	# self.trim_right and self.trim_up can be negative value
		self.trim_up = None


		self.filter_after_fistar = False #select identidied sources within specfic region
		self.filter_after_daofind = False #select identidied sources within specfic region
		self.edge_xl = None		#
		self.edge_xu = None
		self.edge_yl = None
		self.edge_yu = None


		self.magnitude_cut_before_grmatch = False
		self.reference_mag_min = None	#the faint end of reference magnitudes for calibration
		self.reference_mag_max = None	#the bright end of reference magnitudes for calibration
		self.reference_mag_min_note = 'the magnitude cut in the faint end'
		self.reference_mag_max_note = 'the magnitude cut in the bright end'

		self.reference_magerr_cut = None	#the threshold for filtering magnitudes with large uncertainties
		self.input_magerr_cut = None

		self.stdref_tpl_match_criteria = 10	# stdref and tpl match in 'surrouding_search' method
		self.tpl_obs_match_criteria = 10

		self.psfphot_sn_match_criteria = 2	# find corresponding target in psfphot results for SN
		#self.apphot_xy_centroid	= True		#

		self.maxtime_on_single_fwhm = 60        #maximum time allowed for obtaining fwhm from single image
		self.maxtime_on_astrometry = 300	#maxmum time for obtaining astrometry solution

		self.tpl_obs_before_match_display = False
		self.tpl_obs_match_result_display     = False	#display the match result of stars on template and input image
		self.stdref_tpl_match_result_display  = False   #display the match result of stars on template and standard calibration database


		#if os.path.isfile(os.path.join(self.std_ref_dir,self.std_ref_stars_file)):
		#	self.standards = np.loadtxt(os.path.join(self.std_ref_dir,self.std_ref_stars_file))

		self.drop_status = {}
		self.drop_status['0'] = 'safe'
		self.drop_status['1'] = 'bad image quality: background/bad flat'
		self.drop_status['2'] = 'bad image quality: sparse stars/low transparency/star finding error'
		self.drop_status['3'] = 'fwhm'
		self.drop_status['4'] = 'target not found on image'
                self.drop_status['5'] = 'grmatch failure'
                self.drop_status['6'] = 'standard calibration failure'
		self.drop_status['7'] = 'reference match failure'
		self.drop_status['8'] = 'photometry result anomaly/photometry result records less than 10'
		self.drop_status['9'] = 'more than one sources detected at given position'
		self.drop_status['10'] = 'target affccted by cosmic ray or other artifacts'
		self.drop_status['11'] = 'uncharted new phenomenon'
		self.drop_status['12'] = 'get filter information error'
		self.drop_status['13'] = 'image subtraction failure'


		self.dophot_version = 'fortran' #the version of dophot used for PSF photometry: 'fortran' or 'C'
		self.dophot_verbose = False
		self.dophot_image_prepare = 'imcopy' #or softlink

		#light curves
		self.lcs_raw = {}

		self.host_mags_file = 'host_mags.txt'
		self.host_mags_dict = {}
		self.lcs_hostfree = {}


		#insulate the pmfile to prevent from deleting by another process when multiple psf photometry running
		self.insulate_pmfile = False


	def __init_fitsh_grmatch_options(self):
		'''
		See docs/fitsh_grmatch.help for details on the options
		'''
		self.fitsh_grmatch_type = 'point'  #'point'--point matching, 'coord'--coordinate matching, or 'id'--identifier matching

		self.fitsh_grmatch_pointmatch_pars = {}
		self.fitsh_grmatch_pointmatch_pars['--col-ref'] = '1,2' #The index of the first column is always 1
		self.fitsh_grmatch_pointmatch_pars['--col-inp'] = '1,2'
		self.fitsh_grmatch_pointmatch_pars['--order'] = 1 #If the order is A, >= (A+1)*(A+2)/2 valid points are needed to fit the transformation
		self.fitsh_grmatch_pointmatch_pars['--max-distance'] = 1
		self.fitsh_grmatch_pointmatch_pars['--triangulation'] = 'auto,mixed,maxnumber=200'
		self.fitsh_grmatch_pointmatch_pars['--col-ref-ordering'] = -3 #negative sign indicates ascending, small values first
		self.fitsh_grmatch_pointmatch_pars['--col-inp-ordering'] = -3
		self.fitsh_grmatch_pointmatch_pars['--fit'] = None  #iterations=<N>,firstrejection=<F>,sigma=<S>
		self.fitsh_grmatch_pointmatch_pars['--weight'] = None

		self.fitsh_grmatch_coordmatch_pars = {}

		self.fitsh_grmatch_idmatch_pars = {}



	def __init_renew_trigger(self):
		'''
		controller on whether renew some properties
		'''
		self.renew_template_xys = False			#
		self.renew_target_xys = False
		self.renew_aperture_photometry = False		#renew the instmag of self.photometry_info['xxx']['instmag']
		self.renew_aperture_photometry_retfile = False	#renew the aperture photometry result file
		self.renew_psf_photometry = False
		self.renew_psf_photometry_retfile = False
		self.renew_relative_mag = False
		self.renew_stdcal_mag = False
		self.renew_std_ref_match = False
		self.renew_stdcal_xytran = False
		self.renew_standards = False

		self.renew_apphot_parfile = True #overwrite existing iraf-apphot parameter file?

	def __init_iraf_apphot_config(self):
		'''
		parameters init for iraf.apphot
		'''
		self.apphot_iraf_options = {}
		self.apphot_iraf_calgorithm = 'centroid' #the default centering algorithm for apphot, other options 'gauss', 'ofilter'
		self.apphot_iraf_datamin = None
		self.apphot_iraf_datamax = None

		self.apphot_iraf_options['fwhmpsf'] = None
		self.apphot_iraf_options['app'] = 8
		self.apphot_iraf_options['skyin'] = 10
		self.apphot_iraf_options['skywidth'] = 20
		self.apphot_iraf_options['sky_sigma_down'] = 3
		self.apphot_iraf_options['sky_sigma_up'] = 3
		self.apphot_iraf_options['def_zeropt'] = 25

		self.apphot_iraf_autoscale_aperture = True
		self.apphot_iraf_app_nfwhm = 2
		self.apphot_iraf_skyin_nfwhm = 3
		self.apphot_iraf_skywidth_nfwhm = 3


	def __init_dophot_pm(self):
		'''
		initiate the pm data for dophot_C and dophot_fortran
		'''
		self.dophot_C_pm = OrderedDict()

		self.dophot_C_pm['FWHM']    = None    #Approx FWHM of objects (pixels) along major axis.
		self.dophot_C_pm['SKY']     = None
		self.dophot_C_pm['EPERDN']  = None
		self.dophot_C_pm['RDNOISE'] = None
		self.dophot_C_pm['AXIS_RATIO']  = 1.0 #The initial guess of the ratio of the minor axis divided by the major axis for star objects.
		self.dophot_C_pm['TILT']        = 0.0 #The initial guess of the position angle of the major axis with respect to the positive x-axis for a typical unblended stellar object on an image

		self.dophot_C_pm['APBOX_X']     = 16.0#The sizes of the sides of the aperture photometry box in the x direction.
		self.dophot_C_pm['APBOX_Y']     = 16.0#The sizes of the sides of the aperture photometry box in the y direction.
		self.dophot_C_pm['MASKBOX_X']   = 5   # Size of mask box size in x.
		self.dophot_C_pm['MASKBOX_Y']   = 5   # Size of mask box size in y.
		self.dophot_C_pm['NFITBOX_X']   = 12.0#The sizes of the sides of the fitting box used for finding objects in the x direction.
		self.dophot_C_pm['NFITBOX_Y']   = 12.0#The sizes of the sides of the fitting box used for finding objects in the y direction.

		self.dophot_C_pm['IBOTTOM'] = -50        #Lowest allowed data value in data numbers.
		self.dophot_C_pm['ITOP'] = 65535         #Maximum allowed data value in data numbers.
		self.dophot_C_pm['THRESHMIN'] = 100.0    #The value above sky of the lowest threshold that DoPHOT uses to search for objects, in DN
		self.dophot_C_pm['THRESHMAX'] = 40000.0  #The maximum **possible** first threshold that DoPHOT uses to search for objects, in DN.
		self.dophot_C_pm['THRESHDEC'] = 1.0      #Threshold decrement in powers-of-2. For THRESHMAX=40000,THRESHMIN=100, the thresholds (relative to the sky value) would be
							 #25600,12800,6400,3200,1600,800,400,200,100
		self.dophot_C_pm['DOFINALFIT'] = 'YES'   #Do a last fitting iteration after all objects found  fully fits faintest objects of type 3, but slows DoPHOT
		self.dophot_C_pm['THRESHEMP'] = 0.0      #Fit empirical PSF at and below this value.
		self.dophot_C_pm['MAX_SOUGHT'] = 32768   #Quit after this number of improved stars.
		self.dophot_C_pm['MAX_PERF']  = 2000     #Only average up to this number of stars to get shape parameters.
		self.dophot_C_pm['RANGE_MAG'] = 30.      #Sets threshhold some magnitudes fainter than brightest fit.

		self.dophot_C_pm['AUTOSCALE'] = 'NO'     #Auto-scaling of sizes by FWHM.
		self.dophot_C_pm['AUTOTHRESH'] = 'NO'    #Auto-scaling of thresholds.
		self.dophot_C_pm['SCALEFITBOX'] = 3.0    #Size of fit box in units of FWHM.
		self.dophot_C_pm['FITBOXMIN'] = 5.0      #Smallest allowed fit box size.
		self.dophot_C_pm['SCALEAPBOX'] = 6.0     #Size of aperture phot box in units of FWHM.
		self.dophot_C_pm['APBOXMIN'] = 7.0       #Smallest allowed aperture phot box size.
		self.dophot_C_pm['SCALEMASKBOX'] = 1.5   #Size of mask box in units of FWHM.
		self.dophot_C_pm['AMASKBOXMIN'] = 5.0    #Smallest allowed mask box size.
		self.dophot_C_pm['SIGMAIBOTTOM'] = 10.0  #Level of IBOTTOM below sky in units of noise.
		self.dophot_C_pm['SIGMATHRESHMIN'] = 2.0 #Level of THRESHMIN above sky in units of noise.
		self.dophot_C_pm['FIXPOS'] = 'NO'        #Fix star positions?

		self.dophot_C_pm['PARAMS_DEFAULT'] = 'conf/paramdefault_dophot_C'
		self.dophot_C_pm['PARAMS_OUT']    = None            #Output parameters file name.
		self.dophot_C_pm['IMAGE_IN']      = None            #Input image name.
		self.dophot_C_pm['IMAGE_OUT']     = None            #Output image name.
		self.dophot_C_pm['EMP_SUBRAS_OUT']= None            #Empirical PSF subraster (most recent).
		self.dophot_C_pm['OBJECTS_IN']    = None            #Input object list file name.
		self.dophot_C_pm['OBJECTS_OUT']   = None            #Output object list file name.
		self.dophot_C_pm['SHADOWFILE_IN'] = None            #Input shadow file name.
		self.dophot_C_pm['SHADOWFILE_OUT']= None            #Output shadow file name.
		self.dophot_C_pm['ERRORS_OUT']    = None            #Errors of fit to be output if shadow file is requested
		self.dophot_C_pm['LOGFILE'] = None                #Log file name.  'TERM' for screen but failed in practice
		self.dophot_C_pm['LOGVERBOSITY'] = 1

		self.dophot_C_pm['ICRIT']      = 10        #Obliterate if # of pixels > ITOP exceeds this.
		self.dophot_C_pm['CENTINTMAX'] = 40000.0   #Obliterate if central intensity exceeds this.
		self.dophot_C_pm['CTPERSAT']   = 6.0e4     #Assumed intensity for saturated pixels.
		self.dophot_C_pm['NBADLEFT']   = 0         #Ignore pixels closer to the left edge than this.
		self.dophot_C_pm['NBADRIGHT']  = 0         #Ignore pixels closer to the right edge than this.
		self.dophot_C_pm['NBADTOP']    = 0         #Ignore pixels closer to the top edge than this.
		self.dophot_C_pm['NBADBOT']    = 0         #Ignore pixels closer to the bottom edge than this.

		self.dophot_C_pm['PSFTYPE']     = 'PGAUSS'   #PSF type: (PGAUSS, GAUSS, EXTPGAUSS)
		self.dophot_C_pm['SKYTYPE']     = 'PLANE'    #SKY type: (PLANE, HUBBLE, MEDIAN)
		self.dophot_C_pm['OBJTYPE_OUT'] = 'COMPLETE' #INTERNAL, COMPLETE, INCOMPLETE
		self.dophot_C_pm['OBJTYPE_IN']  = 'COMPLETE' #INTERNAL, COMPLETE, INCOMPLETE

		self.dophot_C_pm['EMP_STAR_X'] = 0         #X position of empirical template.
		self.dophot_C_pm['EMP_STAR_Y'] = 0         #Y position of empirical template.
		self.dophot_C_pm['EMP_STAR_Z'] = 0         #Central intensity of empirical template.

		self.dophot_fortran_pm = OrderedDict()
        	self.dophot_fortran_pm['SKY']     =  None	#Approximate mean sky value in data numbers.
        	self.dophot_fortran_pm['FWHM']    =  None	#Approx FWHM of objects (pixels) along major axis.
        	self.dophot_fortran_pm['EPERDN']  =  None	#Electrons per data number.
        	self.dophot_fortran_pm['RDNOISE'] =  None	#Readout noise in electrons.
        	self.dophot_fortran_pm['TOP'] = 65535		#Maximum allowed data value in data numbers.
		self.dophot_fortran_pm['AUTOTHRESH'] = 'NO'		#psf photometry control parameters
		self.dophot_fortran_pm['THRESHMIN'] = 100.0    #The value above sky of the lowest threshold that DoPHOT uses to search for objects, in DN
		self.dophot_fortran_pm['THRESHMAX'] = 40000.0  #The maximum **possible** first threshold that DoPHOT uses to search for objects, in DN.
		self.dophot_fortran_pm['PARAMS_DEFAULT'] = 'conf/paramdefault_dophot_fortran'
		self.dophot_fortran_pm['PSFTYPE']        = 'PGAUSS' #for C version, there is a new model EXTPGAUSS
		self.dophot_fortran_pm['OBJTYPE_OUT']    = 'COMPLETE' #INTERNAL, COMPLETE, INCOMPLETE
		self.dophot_fortran_pm['OBJTYPE_IN']     = 'COMPLETE' #INTERNAL, COMPLETE, INCOMPLETE


	def __init_photometry_info(self):
		'''
		initiate the photometry info dictionary
		'''
		self.photometry_info = OrderedDict()
		self.photometry_info_keys = ['name','realimg','flt','obstime','exptime','bitpix','flt2int','fixbad','rmCR',\
					     'fwhm','bkg','airmass','nstar','template','x','y','instmag','instmagerr','relmag','relmagerr', 'magzpt','magzpterr', 'calmag','calmagerr','drop']
		self.photometry_info_dtypes = ['S10','S60','S10','f8','f8','i4','i4','i4','i4','f8','f8','f8','i4','i4','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','i4']

		self.photometry_info_init_values = {}
		self.photometry_info_init_values['realimg'] = '?'
		self.photometry_info_init_values['flt'] = '?'
		self.photometry_info_init_values['obstime'] = 0
		self.photometry_info_init_values['exptime'] = 0
		self.photometry_info_init_values['bitpix'] = 0
                self.photometry_info_init_values['flt2int'] = 0
                self.photometry_info_init_values['fixbad'] = 0
                self.photometry_info_init_values['rmCR'] = 0	#remove cosmic ray
		self.photometry_info_init_values['fwhm'] = 0
                self.photometry_info_init_values['bkg'] = 0
		self.photometry_info_init_values['airmass'] = 0
                self.photometry_info_init_values['nstar'] = 0
                self.photometry_info_init_values['template'] = 0
                self.photometry_info_init_values['x'] = 0.0
                self.photometry_info_init_values['y'] = 0.0
                self.photometry_info_init_values['instmag'] = 99.99
                self.photometry_info_init_values['instmagerr'] = 99.99
                self.photometry_info_init_values['relmag'] = 99.99
                self.photometry_info_init_values['relmagerr'] = 99.99
                self.photometry_info_init_values['magzpt'] = 99.99
                self.photometry_info_init_values['magzpterr'] = 99.99
                self.photometry_info_init_values['calmag'] = 99.99
                self.photometry_info_init_values['calmagerr'] = 99.99
                self.photometry_info_init_values['drop'] = 0

		for img in self.images:
			self.photometry_info[img] = OrderedDict()
                	self.photometry_info[img]['name'] = img
			self.photometry_info[img]['realimg']   = self.photometry_info_init_values['realimg']
			self.photometry_info[img]['flt']       = self.photometry_info_init_values['flt']
			self.photometry_info[img]['obstime']   = self.photometry_info_init_values['obstime']
			self.photometry_info[img]['exptime']   = self.photometry_info_init_values['exptime']
			self.photometry_info[img]['bitpix']    = self.photometry_info_init_values['bitpix']
			self.photometry_info[img]['flt2int']   = self.photometry_info_init_values['flt2int']
			self.photometry_info[img]['fixbad']    = self.photometry_info_init_values['fixbad']
			self.photometry_info[img]['rmCR']      = self.photometry_info_init_values['rmCR']
			self.photometry_info[img]['fwhm']      = self.photometry_info_init_values['fwhm']
			self.photometry_info[img]['bkg']       = self.photometry_info_init_values['bkg']
			self.photometry_info[img]['airmass']       = self.photometry_info_init_values['airmass']
			self.photometry_info[img]['nstar']     = self.photometry_info_init_values['nstar']
			self.photometry_info[img]['template']  = self.photometry_info_init_values['template']
			self.photometry_info[img]['x']         = self.photometry_info_init_values['x']
			self.photometry_info[img]['y']         = self.photometry_info_init_values['y']
			self.photometry_info[img]['instmag']   = self.photometry_info_init_values['instmag']
			self.photometry_info[img]['instmagerr']= self.photometry_info_init_values['instmagerr']
			self.photometry_info[img]['relmag']    = self.photometry_info_init_values['relmag']
			self.photometry_info[img]['relmagerr'] = self.photometry_info_init_values['relmagerr']
			self.photometry_info[img]['magzpt'] = self.photometry_info_init_values['magzpt']
			self.photometry_info[img]['magzpterr'] = self.photometry_info_init_values['magzpterr']
			self.photometry_info[img]['calmag']    = self.photometry_info_init_values['calmag']
			self.photometry_info[img]['calmagerr'] = self.photometry_info_init_values['calmagerr']
			self.photometry_info[img]['drop']      = self.photometry_info_init_values['drop']

		self.anomaly_signs = {}
		self.anomaly_signs['bkg'] =  -9999.99
		self.anomaly_signs['fwhm'] = -99.99
		self.anomaly_signs['nstar'] = -9999


	def __prepare_parameter_files(self):
		'''
		prepare a copy of parameter files for the working target
		'''
		if not os.path.exists(self.parafile_dir):
			os.mkdir(self.parafile_dir)

		for parafile in os.listdir(self.parafile_template_dir):
			parafile_abs_from = os.path.join(self.parafile_template_dir,parafile)
			parafile_abs_to = os.path.join(self.parafile_dir,parafile)
			if not os.path.isfile(parafile_abs_to):
			#self.__delete_file_if_exist(parafile_abs_to)
				shutil.copy(parafile_abs_from,parafile_abs_to)

	def check_notes(self ):
		'''
		there is a notes file under result directory, try to look into that for reduction notes
		'''
		if not os.path.exists(self.notesfile):
			print "Notes file not exists..."
		else:
			notes = open(self.notesfile).readlines()
			print notes

	def add_note(self, note):
		'''
		add note to self.readme
		'''
		self.readme = self.readme + note+'\n'

	def save_notes(self, outfile=None):
		'''
		save notes in self.readme to the outfile; default outfile: self.notesfile, reduction README file
		'''
		if outfile is None:
			outfile = self.notesfile
		fid = open(outfile, 'awt')
		fid.write(self.readme)
		fid.close()

	def get_single_star_SNR(self,flt,photometry_function = 'iraf'):
		'''
		simple statistics on SNR of the photometry targets in given band
		'''
		obstimes = []
		snrs = []

		for img in self.images.keys():
			if self.photometry_info[img]['drop'] != 0:
				continue
			if self.photometry_info[img]['flt'] != flt:
				continue

			obstime = self.photometry_info[img]['obstime']
			snr = self.__single_star_SNR(img,photometry_function=photometry_function)

			obstimes.append(obstime)
			snrs.append(snr)

		SNR_table = Table()
		obstime_col = Column(data=obstimes,name='obstime')
		snr_col = Column(data=snrs,name='snr')
		SNR_table.add_columns([obstime_col,snr_col])

		return SNR_table




	def __single_star_SNR(self,image_key, photometry_function='iraf', which_dir = 'raw_image', fixed_app_size = False):
		'''
		estimate the signal to noise ratio of photometry target
		The SNR calculation is overly simplified
		Please refer to /Users/chenping/Documents/instrument_data_analysis/SNR_CCD_optical.pdf for improvement

		INPUTS:
			image_key:
			photometry_function:
			which_dir:
			fixed_app_size:
		'''

		telname = self.current_telescope
		GAIN = self.base.telescopes_info[telname]['gain']
		RDNS = self.base.telescopes_info[telname]['readnoise']

		#if gain is not found in registered telescope information, fill 1.0 for this
		if GAIN == '':
			GAIN = 1.0
		else:
			GAIN = float(GAIN)

		if RDNS == '':
			RDNS = 0
		else:
			RDNS = float(RDNS)

		if photometry_function == 'iraf':
                        x = self.photometry_info[image_key]['x']
                        y = self.photometry_info[image_key]['y']
                        xy = np.array([x,y]).reshape((1,2))
			image = self.__get_internal_image(image_key, which_dir=which_dir)

			self.get_apphot_iraf_parameters(image_key)
                        options = self.apphot_iraf_options.copy()

                       	photret_singlept = self.__aperture_photometry_apphot_iraf_single_image(image,xy,options)
                        #xpos, ypos, flux, bkg,bkg_stdev, sky_counts, target_area,mag,mag_err

			photoret_singlept = photret_singlept.ravel()
                        print photret_singlept

                        if len(photret_singlept) == 0:
                                SNR = -99.99
                                return SNR

			flux = photret_singlept[0,2]
			bkg = photret_singlept[0,3]
			target_area = photret_singlept[0,6]

			bkg_total = bkg*target_area
			signal = flux*GAIN
			noise = np.sqrt(bkg_total*GAIN + flux*GAIN)

                        print xy
			print "Gain:",GAIN
			print "flux:",flux
			print "bkg_total:", bkg_total
			print "signal:", signal
			print "noise:", noise

			SNR = signal/noise

		return SNR



	def get_photometry_uncertainty_flt(self,flt, phot_method = 'apphot', mode='relmag', updatetable=1):
		'''
		get photometry unertainties for source of different magnitudes

		This is used to estimate systematic uncertainties by studying the scatter of magnitudes from all images for each star

		How? Register all stars in the template image and then find all matches in other images for each star. Then study the scatter of stars with different magnitudes

		'''
		self.__find_template_imagekey(updatetable=updatetable)
		tpl_image = self.templates[flt]
		tpl_imgkey = tpl_image.split('.')[0]
		if updatetable:
               		self.__dict2table()
                photret_table_all = self.result_table_newsaved

		photret_flt = self.__select_rows_from_table(photret_table_all,'flt',flt)
		photret_flt_nondrop = self.__select_rows_from_table(photret_flt,'drop',0)

		if phot_method == 'apphot':
			photret_dir = self.aperture_photometry_dir
			photfile_suffix = self.apphot_ret_file_suffix
		elif phot_method == 'psfphot':
			photfile_suffix = self.psfphot_ret_file_suffix
			photret_dir = self.psf_photometry_dir
		else:
			raise IOError("Invalid input for phot_method...")

		tpl_apphot = os.path.join(photret_dir, tpl_imgkey+photfile_suffix)


		xys_dict = {}
		mags_dict = {}

		offset_tpl = self.photometry_info[tpl_image][mode] - self.photometry_info[tpl_image]['instmag']
		tpl_data = np.loadtxt(tpl_apphot)

		print tpl_data

		for i,phot_single in enumerate(tpl_data):
			mags_dict[i] = [phot_single[2]+offset_tpl]
			xys_dict[i] = [phot_single[0],phot_single[1]]

		for image in photret_flt_nondrop['name']:
			if image == tpl_image:
				continue

			if mode != 'calmag' and mode !='relmag':
				raise IOError("Invlalid input for mode")

			offset_mag = self.photometry_info[image][mode] - self.photometry_info[image]['instmag']

			imgkey = image.split('.')[0]
			this_tpl_match_file = os.path.join(photret_dir, imgkey+'_tpl.match')
			this_phot_data = np.loadtxt(this_tpl_match_file)

			tpl_xys_match = this_phot_data[:,[0,1]]

			for dict_key in xys_dict.keys():
				xy = xys_dict[dict_key]
				exist_flag,xy_match,match_index = self.__find_corresponding(tpl_xys_match,xy,0.5)
				if not exist_flag:
					continue

				mag = this_phot_data[match_index,6] + offset_mag
				mags_dict[dict_key].append(mag)

		mags = []
		stds  = []
		for dict_key in mags_dict.keys():
			mags_this_star = mags_dict[dict_key]


			if len(mags_this_star)<4:
				continue
			mags.append(np.mean(mags_this_star))
			stds.append(np.std(mags_this_star))

		mags_array = np.array(mags).reshape(len(mags),1)
		stds_array = np.array(stds).reshape(len(stds),1)
		out = np.hstack((mags_array,stds_array))

		outsort = out[np.argsort(out[:,0]),:]
		#print outsort

		out_file_name = os.path.join(self.result_dir,'mags_stds_'+flt+'.txt')
		self.__delete_file_if_exist(out_file_name)

		np.savetxt(out_file_name,outsort,fmt="%5.2f %5.2f")

		plt.plot(mags,stds,'o')
		plt.xlabel('mag')
		plt.ylabel('mag scatter')
		plt.show()



	def get_light_curves_ref_stars(self, flt, mag_min_ref=None, mag_max_ref=None, phot_method='apphot', lctype='calmag', only_caldata=True, verbose=0):
		'''
		extract light curves for all reference stars which meet the given criteria


		INPUTS:
			flt:
			mag_min_ref: faintest end
			mag_max_ref: brightest end
			phot_method: 'apphot' or 'psfphot'
			lctype: 'relmag' or 'calmag'
		'''

		if flt not in self.templates.keys():
			self.__find_template_imagekey(updatetable=updatetable)

		tpl_image = self.templates[flt]
		if verbose:
			print "the template image: %s %s"%(tpl_image, self.images[tpl_image])
		tpl_imgkey = tpl_image.split('.')[0]

		if phot_method == 'apphot':
			photret_dir = self.aperture_photometry_dir
			photfile_suffix = self.apphot_ret_file_suffix
		elif phot_method == 'psfphot':
			photfile_suffix = self.psfphot_ret_file_suffix
			photret_dir = self.psf_photometry_dir
		else:
			raise IOError("Invalid input for phot_method...")

		tpl_photfile = os.path.join(photret_dir, tpl_imgkey+photfile_suffix)
		if verbose:
			print "The photometry file of tempalte image: %s"%tpl_photfile

		tpl_photdata = np.loadtxt(tpl_photfile)
		if verbose:
			print "photometry data of tempalte image:\n"
			print tpl_photdata

		lcs_ref_dir = os.path.join(self.warehouse_dir, 'lcs_ref')
		if not os.path.exists(lcs_ref_dir):
			os.mkdir(lcs_ref_dir)

		lcs_ref_spm_dir = os.path.join(lcs_ref_dir, phot_method)
		if not os.path.exists(lcs_ref_spm_dir):
			os.mkdir(lcs_ref_spm_dir)

		for photline in tpl_photdata:
			x,y,mag,magerr=photline

			if mag_min_ref is not None:
				if mag > mag_min_ref:
					continue

			if mag_max_ref is not None:
				if mag < mag_max_ref:
					continue

			reflcfile =os.path.join( lcs_ref_spm_dir,  'ref_'+flt+'_' + str(np.round(x,1)).replace('.','p') + '_' + str(np.round(y,1)).replace('.', 'p') + '.txt')

			lcdata = self.get_light_curve_for_control_star(flt, x, y, phot_method=phot_method, lctype=lctype, only_caldata=only_caldata)


			calstar_index = np.where(lcdata['mag']!=99.99)[0]
			if len(calstar_index)>1:
				lcdata.write(reflcfile, format='ascii.commented_header')
			else:
				print "the star at image coordinate (%s,%s) only exist on template image"%(x,y)


	def get_light_curve_for_control_star(self,flt,tpl_x,tpl_y, phot_method='apphot', lctype='calmag', only_caldata=True,  verbose=0, updatetable=1):
		'''
		extract light curve for a given star with position on reference image (tpl_x, tpl_y)

		INPUTS:
			flt:
			tpl_x:
			tpl_y:
			phot_method: 'apphot' or 'psfphot'
			lctype: 'calmag' or 'relmag'
			only_caldata: only extract control which have been used in the relative calibration
		'''
		self.__find_template_imagekey(updatetable=updatetable)
		tpl_image = self.templates[flt]
		tpl_imgkey = tpl_image.split('.')[0]
		if updatetable:
                	self.__dict2table()
                photret_table_all = self.result_table_newsaved

		photret_flt = self.__select_rows_from_table(photret_table_all,'flt',flt)
		photret_flt_nondrop = self.__select_rows_from_table(photret_flt,'drop',0)



		if phot_method == 'apphot':
			photret_dir = self.aperture_photometry_dir
			photfile_suffix = self.apphot_ret_file_suffix
		elif phot_method == 'psfphot':
			photfile_suffix = self.psfphot_ret_file_suffix
			photret_dir = self.psf_photometry_dir
		else:
			raise IOError("Invalid input for phot_method...")

		tpl_apphot = os.path.join(photret_dir, tpl_imgkey+photfile_suffix)


		tpl_data = np.loadtxt(tpl_apphot)
		if verbose:
			print tpl_data

		mags = []
		magerrs = []
		obstimes = []
		for i,image in enumerate(photret_flt_nondrop['name']):

			obstimes.append(photret_flt_nondrop['obstime'][i])

			offset_mag = self.photometry_info[image][lctype] - self.photometry_info[image]['instmag']

			imgkey = image.split('.')[0]
			if imgkey == tpl_imgkey:
				if phot_method == 'apphot':
					tpl_phot_suffix = '.apphot'
				else:
					tpl_phot_suffix = '.psfphot'
				this_tpl_match_file = os.path.join(photret_dir,imgkey+tpl_phot_suffix)

			else:
				if only_caldata:
					this_tpl_match_file = os.path.join(photret_dir, imgkey+'_tpl_caldata.match')
				else:
					this_tpl_match_file = os.path.join(photret_dir, imgkey+'_tpl.match')

			this_phot_data = np.loadtxt(this_tpl_match_file)

			tpl_xys_match = this_phot_data[:,[0,1]]

			xy = [tpl_x,tpl_y]
			exist_flag,xy_match,match_index = self.__find_corresponding(tpl_xys_match,xy,1.5)
			if not exist_flag:
				mags.append(99.99)
				magerrs.append(99.99)
				continue
			if imgkey == tpl_imgkey:
				mag = this_phot_data[match_index,2] + offset_mag
				magerr_inst = this_phot_data[match_index,3]
			else:
				mag = this_phot_data[match_index,6] + offset_mag
				magerr_inst = this_phot_data[match_index,7]

			magerr_relcal = self.photometry_info[image]['relmagerr']
			magerr = np.sqrt(magerr_inst**2 + magerr_relcal**2)


			mags.append(mag[0])
			magerrs.append(magerr[0])

		ref_star_data  = [obstimes, mags, magerrs]
		control_star_table = Table(ref_star_data, names = ['obstime', 'mag', 'magerr'])

		#control_star_table = Table()
		#obstime_col = Column(data=obstimes,name='obstime')
		#mag_col     = Column(data=mags,name='mag')
		#magerr_col  = Column(data=magerrs,name='magerr')
		#control_star_table.add_columns([obstime_col,mag_col,magerr_col])

		return control_star_table


	def check_battle_damage(self,flts='all', updatetable=1):
		'''
		Check images which have been droped with non-zero drop flag
		'''
		if updatetable:
			self.__dict2table()

		if flts == 'all':
			self.__find_template_imagekey(updatetable=updatetable)
			flts_tocheck = self.templates.keys()

			if len(flts_tocheck) == 0:
				flts_tocheck = np.unique(self.result_table_newsaved['flt'])

		elif isinstance(flts,str):
			flts_tocheck = [flts]
		elif isinstance(flts,list):
			flts_tocheck = flts
		else:
			raise IOError('Invalid input for flts...')


		droped_images = []
		for flt in flts_tocheck:

			print "Check %s band image:"%flt

			info_flt = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)
			for img in info_flt['name']:
				drop_status =  self.photometry_info[img]['drop']
				if drop_status != 0:
					print img,drop_status
					droped_images.append(img)
		return droped_images


	def block_specific_image(self,colname,criteria,mode='eq', updatetable=1):
		'''
		block images which you don't want to process

		other mode rather than 'eq' is still under construction
		'''
		if updatetable:
			self.__dict2table()
		colnames = self.result_table_newsaved.colnames

		if colname not in colnames:
			raise IOError("%s not exist in table"%colname)


		for i,value in enumerate(self.result_table_newsaved[colname]):
			imgkey = self.result_table_newsaved['name'][i]
			print imgkey
			if mode == 'eq':
				if value == criteria:
					self.photometry_info[imgkey]['drop'] = 10



	def bring_dead_back_to_life(self,drop_flag='all', updatetable=1):
		'''
		renew drop flag for images with drop signal as 'drop_flag'

		Inputs:
			drop_flag: default value is 'all' which bring all dead to life
				   You can rescue specific images by specify drop_flags,
				   for example drop_flag = '3'
				   Please see self.drop_status for the details on drop flags

		'''
		drop_flags_allowed = self.drop_status.keys()
		if updatetable:
			self.__dict2table()
		imginfo_table = self.result_table_newsaved

		if drop_flag == 'all':
			save_images = imginfo_table['name']
		elif str(drop_flag) in drop_flags_allowed:
			info_table_specific_drop = self.__select_rows_from_table(imginfo_table,'drop',drop_flag)
			N_save = len(info_table_specific_drop)
			save_images = info_table_specific_drop['name']
		else:
			raise IOError('Sorry, invalid input for drop_flag')

		if drop_flag != 'all' and N_save == 1:
			self.photometry_info[save_images]['drop'] = 0
		else:
			for img in save_images:
				self.photometry_info[img]['drop'] = 0



	def init_phot_info_specific_col(self,colname,flts='all', updatetable=1):
		'''
		initiate specific info of photometry infomation

		Inputs:
			colname: check self.photometry_info_keys for valid inputs
			flts:	default 'all',
				non-default input should be a list or np.ndarray containing the
				ilters, eg flts= ['B'] or flts = ['B','V']
		'''
		if updatetable:
			self.__dict2table()
		imginfo_table = self.result_table_newsaved


		if flts == 'all':
			flts = np.unique(imginfo_table['flt'])

		if isinstance(flts,str):
			flts = [flts]

		for flt in flts:
			init_images = self.__select_rows_from_table(imginfo_table,'flt',flt)
			for img in init_images['name']:
				self.photometry_info[img][colname] = self.photometry_info_init_values[colname]


	def select_photometry_measurement_given_region(self, image_key,xl=None,xh=None,yl=None,yh=None, photometry_method='apphot', outfile=None):
		'''
		select photometry points within [xl:xh,yl:yh]

		The selection works on the phootometry output file (self.apphot_ret_file_suffix for aperture photometry and self.psfphot_ret_file_suffix for psf photometry)
		and the selected result will be written into 'outfile' or overwritten into original photometry output file if 'outfile' is not provided

		'''
		imgkey = image_key.split('.')[0]
		if photometry_method == 'apphot':
			filesuffix  = self.apphot_ret_file_suffix
			photret_dir = self.aperture_photometry_dir
		elif photometry_method == 'psfphot':
			filesuffix = self.psfphot_ret_file_suffix
			photret_dir = self.psf_photometry_dir
		else:
			raise ValueError("not support")

		infile = os.path.join(photret_dir, imgkey+filesuffix)
		data = np.loadtxt(infile)

		if xl is not None:
			data = data[data[:,0]>xl,:]
		if xh is not None:
			data = data[data[:,0]<xh, :]

		if yl is not None:
			data = data[data[:,1]>yl,:]
		if yh is not None:
			data = data[data[:,1]<yh, :]

		if outfile is None:
			outfile = infile
			print "orginal file %s will be overwrited"%infile

		np.savetxt(outfile, data, fmt='%6.2f %6.2f %5.2f %5.2f')



	def smarts_ccd_remove_photometry_in_bad_region(self,xlow=105, ylow=30):
		'''
		SMARTS 1.3m ANDICAM CCD has some region which are not used
		'''
		photmethod = self.photometry_method
		for img in self.images.keys():
			if self.photometry_info[img]['drop'] == 0:
				self.select_photometry_measurement_given_region(img, xl=xlow, yl=ylow, photometry_method=photmethod)




	def delete_photometry_measurement_given_region(self, imgkey, x_center, y_center, radius, infile = None, xcol_infile= 0, ycol_infile = 1,  photometry_method='apphot', outfile=None):

		'''
		delete photometry points within circle at given center with given radius

		'''
		if infile is None:

			if photometry_method == 'apphot':
				filesuffix  = self.apphot_ret_file_suffix
				photret_dir = self.aperture_photometry_dir
			elif photometry_method == 'psfphot':
				filesuffix = self.psfphot_ret_file_suffix
				photret_dir  = self.psf_photometry_dir
			else:
				raise ValueError("not support")

			infile = os.path.join(photret_dir, imgkey+filesuffix)


		data = np.loadtxt(infile)
		x_image = data[:,xcol_infile]
		y_image = data[:,ycol_infile]

		dists = np.sqrt((x_image - x_center)**2 + (y_image - y_center)**2)
		mask = dists > radius
		data = data[mask,:]
		if outfile is None:
			outfile = infile
			print "orginal file %s will be overwrited"%infile

		M,N = data.shape
		np.savetxt(outfile, data, fmt='%6.2f '*(N-1)+'%6.2f')




	def check_images_picked_from_lc_flt(self,flt,lctype = 'calmag', image_dir = 'raw_image', updatetable=1):
		'''
		display the image corresponding to the point (obstime,mag) in light curve

		lctype:	'calmag', 'relmag', 'instmag'
		'''
		if updatetable:
			self.__dict2table()
		photret_table_all = self.result_table_newsaved

		photret_table_nondrop = self.__select_rows_from_table(photret_table_all,'drop',0)
		flt_photret = self.__select_rows_from_table(photret_table_nondrop, 'flt', flt)
		obstimes = flt_photret['obstime']

		if lctype not in ['calmag','relmag','instmag']:
			raise IOError("Invalid input for lctype, calmag, relmag or instmag are allowed")

		mags = flt_photret[lctype]
		magerrs = flt_photret[lctype+'err']
		obstimes_selected,mags_selected,index_selected=mouse_pick(obstimes,mags,err=magerrs) #mouse_pick can pick up multiple points

		for obstime_wanted,mag_wanted in zip(obstimes_selected, mags_selected):

			table_selected_from_obstime = self.__select_rows_from_table(flt_photret,'obstime',obstime_wanted,mode='lgt',criteria=0.001)
			table_selected_from_mag = self.__select_rows_from_table(flt_photret, lctype, mag_wanted,mode='lgt',criteria=0.01)

			for name in table_selected_from_obstime['name']:
				if name in table_selected_from_mag['name']:
					final_selected_imgkey = name

			d = self.display_image_with_ds9(final_selected_imgkey,which_dir=image_dir)
			x = self.photometry_info[final_selected_imgkey]['x']
	        	y = self.photometry_info[final_selected_imgkey]['y']
        		xy_reg = [x,y]
	        	self.__load_source_regions_to_ds9(d, xy_reg, radius=25, color='red', width=2)
			drop_signal_default = 0
			drop_signal = raw_input("Enter the drop signal for image %s (default: %s):" %(final_selected_imgkey,drop_signal_default)) or drop_signal_default

			if not str(drop_signal).isdigit():
				raise IOError("Invalid input for drop_signal")

			drop_sinal = int(drop_signal)
			if drop_signal:
				self.photometry_info[final_selected_imgkey]['drop'] = drop_signal




	def lc_plot(self,mode,flts = 'all',bining =False, bin_width = 0.5,remove_inbin_outliers=False, plot_droped = False,plot_bined_only = False, updatetable=1):
		'''
		Inputs:
			mode: instmag; relmag; calmag
			flts: 'all' or 'V' or ['B','V']
			bining: bin the data of not if multiple measurements obtained within the one day
			bin_width: the widt for the bin
		'''
		if updatetable:
			self.__dict2table()
		self.__find_template_imagekey(updatetable=updatetable)


		if mode == "instmag":
			magkey = 'instmag'
			magerrkey = 'instmagerr'
		elif mode == "relmag":
			magkey = 'relmag'
			magerrkey = 'relmagerr'
		elif mode == "calmag":
			magkey = 'calmag'
			magerrkey = 'calmagerr'
		else:
			raise KeyError("Acceptable mode inputs are 'instmag','relmag' and 'calmag'")

		if flts == 'all':
			flts_plot = self.templates.keys()
		else:
			if isinstance(flts,str):
				flts  = [flts]
			flts_plot = flts

		for flt in flts_plot:
			print flt
			table_flt =self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)
			if not plot_droped:
				table_flt = self.__select_rows_from_table(table_flt,'drop',0)

			table_flt.sort('obstime')

			obstimes = table_flt['obstime']
			mags = table_flt[magkey]
			magerrs = table_flt[magerrkey]

			if not plot_bined_only:
				plt.errorbar(obstimes,mags,yerr=magerrs,fmt='o')

			if bining:
				obstimes_bined = []
				mags_bined = []
				magerrs_bined = []

				bined_already_index = []
				for i,obstime in enumerate(obstimes):
					if i in bined_already_index:
						continue

					distances = obstimes - obstime
					indexs = np.where(np.abs(distances)< bin_width)[0]
					print indexs
					bined_already_index = np.append(bined_already_index,indexs)

					obstime_bined = np.mean(obstimes[indexs])
					mags_tobin = mags[indexs]
					magerrs_tobin = magerrs[indexs]

					#remove outlier measurements
					if remove_inbin_outliers and len(mags_tobin)>1:
						mags_tobin_nonoutlier,index_keep,stdev = sigma_clipping(mags_tobin,sig=2)
						magerrs_tobin_nonoutlier = magerrs_tobin[index_keep]
					else:
						mags_tobin_nonoutlier = mags_tobin
						magerrs_tobin_nonoutlier = magerrs_tobin

					print mags_tobin_nonoutlier

					weights = 1./magerrs_tobin_nonoutlier**2
					mag_bined = np.sum(weights*mags_tobin_nonoutlier) / np.sum(weights)
					magerr_bined = 1./ np.sqrt(np.sum(weights))

					obstimes_bined.append(obstime_bined)
					mags_bined.append(mag_bined)
					magerrs_bined.append(magerr_bined)
				plt.errorbar(obstimes_bined,mags_bined, yerr=magerrs_bined,fmt='s')

		plt.gca().invert_yaxis()
		plt.legend(flts_plot)
		plt.show()


	def shift_obstime_from_exposure_start_to_exposure_middle(self):
		'''
		when accurate obstimes are required, we need to make sure the reported obstime is at the middle of the exposure
		'''

		for img in self.images.keys():
			self.photometry_info[img]['obstime'] = self.photometry_info[img]['obstime'] + self.photometry_info[img]['exptime']/2.0/3600/24


	def save_lc_data(self,flts='all', binned =False, hostfree=False, bin_width=0.5,remove_inbin_outlier = False, updatetable=1):
		'''
		INPUTS:
			binned: binning the data first before save?
			hostfree: subtract host flux before save?
		'''
		self.__find_template_imagekey(updatetable=updatetable)

		if flts == 'all':
			flts_save = self.templates.keys()
		elif isinstance(flts,str):
			flts_save = [flts]
		elif isinstance(flts,list):
			flts_save = flts
		else:
			raise IOError('Invalid input for flts...')

		for flt in flts_save:
			self.save_lc_data_flt(flt,binned = binned, hostfree=hostfree, bin_width=bin_width,remove_inbin_outlier = remove_inbin_outlier)


	def __get_lc_flt(self, flt):
		'''
		get light curve (t,mag,magerr) from self.photometry_info
		'''
		self.__dict2table()
		tpl_table = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)
		tpl_nondrop_table = self.__select_rows_from_table(tpl_table,'drop',0)

		lc_flt  = Table()
		lc_flt.add_columns([tpl_nondrop_table['obstime'],tpl_nondrop_table['calmag'],tpl_nondrop_table['calmagerr']])
		lc_flt.sort('obstime')

		self.lcs_raw[flt] = lc_flt

		return lc_flt

	def __binning_data(self, obstimes, mags, magerrs, bin_width=0.5, remove_inbin_outlier=False):
		'''
		Data binning
		'''
                obstimes_bined = []
                mags_bined = []
                magerrs_bined = []

                bined_already_index = []
                for i,obstime in enumerate(obstimes):
                        if i in bined_already_index:
                                continue

                        distances = obstimes - obstime
                        #indexs = find(np.abs(distances)< bin_width)
                        indexs = np.where(np.abs(distances)< bin_width)[0]
                        bined_already_index = np.append(bined_already_index,indexs)

                        obstime_bined = np.mean(obstimes[indexs])
                        mags_tobin = mags[indexs]
                        magerrs_tobin = magerrs[indexs]

                        #remove outlier measurements
                        if remove_inbin_outlier and len(mags_tobin)>1:
                                mags_tobin_nonoutlier,index_keep,stdev = sigma_clipping(mags_tobin,sig=2)
                                magerrs_tobin_nonoutlier = magerrs_tobin[index_keep]
                        else:
                                mags_tobin_nonoutlier = mags_tobin
                                magerrs_tobin_nonoutlier = magerrs_tobin
                        weights = 1./magerrs_tobin_nonoutlier**2
                        mag_bined = np.sum(weights*mags_tobin_nonoutlier) / np.sum(weights)
                        magerr_bined = 1./ np.sqrt(np.sum(weights))

                        obstimes_bined.append(obstime_bined)
			mags_bined.append(mag_bined)
                        magerrs_bined.append(magerr_bined)

		N = len(obstimes_bined)
		obstimes_bined_array = np.array(obstimes_bined).reshape(N,1)
		mags_bined_array    = np.array(mags_bined).reshape(N,1)
		magerrs_bined_array = np.array(magerrs_bined).reshape(N,1)
		lc_flt = np.hstack((obstimes_bined_array,mags_bined_array,magerrs_bined_array))

		return lc_flt

	def save_lc_data_flt(self,flt,binned=False, hostfree=False,bin_width=0.5,remove_inbin_outlier=False):
		'''
		save light curve
		'''
		if hostfree:
			self.__get_host_mags_dict()
			mag_host = self.host_mags_dict[flt][0]
			magerr_host = self.host_mags_dict[flt][1]
			lc_flt = self.subtract_host_magnitude_flt(flt, mag_host, magerr_host = magerr_host, A_flt=0, display_result = False)
			flt = flt + '-HFree'
		else:
			lc_flt = self.__get_lc_flt(flt)

		obstimes = lc_flt['obstime']
                mags = lc_flt['calmag']
                magerrs = lc_flt['calmagerr']

		if binned:
			lc_flt = self.__binning_data(obstimes, mags,magerrs,bin_width=bin_width, remove_inbin_outlier=remove_inbin_outlier)
			flt = flt+'-ave'

		sn_name = self.current_sn
		telescope_code =self.base.telescopes_info[self.current_telescope]['code']
		if self.photometry_method == 'apphot':
			lc_filename = sn_name+'_'+telescope_code+'-'+flt+'-AP.txt'
		elif self.photometry_method == 'psfphot':
			lc_filename = sn_name+'_'+telescope_code+'-'+flt+'-PSF.txt'
		else:
			raise IOError("Invalid input for photometry_method...")
		outfile = os.path.join(self.result_dir,lc_filename)
		np.savetxt(outfile,lc_flt,fmt="%12.4f %8.3f %8.3f")


		if self.base.save_data_local_dir != '' and os.path.exists(self.base.save_data_local_dir):
			lc_dir_archive = os.path.join(self.base.save_data_local_dir,sn_name)
			if not os.path.exists(lc_dir_archive):
				os.mkdir(lc_dir_archive)

			outfile = os.path.join(lc_dir_archive, lc_filename)
			np.savetxt(outfile,lc_flt,fmt="%12.4f %8.3f %8.3f")

	def __get_host_mags_dict(self):
		'''
		Read in the magnitudes of host galaxy from the externally prepared and well stored file
		'''
		host_mags_file  = os.path.join( self.std_ref_dir , self.host_mags_file )

		if not os.path.exists(host_mags_file):
			raise IOError("host magnitudes not available: %s"%host_mags_file)

		host_mags = np.loadtxt(host_mags_file, dtype= {'names':('flt','mag','magerr'), 'formats':('S8','f4','f4')})

		for host_flt_mag in host_mags:
			flt = host_flt_mag[0]
			mag  = host_flt_mag[1]
			magerr = host_flt_mag[2]
			self.host_mags_dict[flt] = [mag, magerr]



	def subtract_host_magnitude(self, flts='all', display_result =False ):
		'''
		For nucleus targets, if archive host magnitude exist, we can subtract the host flux to get rough estimation on transient light curve

		'''
		self.__get_host_mags_dict()

		if flts == 'all':
			self.__find_template_imagekey()
			flts = self.templates.keys()

		for flt in flts:
			host_mag = self.host_mags_dict[flt][0]
			host_magerr = self.host_mags_dict[flt][1]

			lc_hostfree = self.subtract_host_magnitude_flt(flt, host_mag, magerr_host = host_magerr, display_result = display_result)

			self.lcs_hostfree[flt]  = lc_hostfree



	def subtract_host_magnitude_flt(self, flt, mag_host, magerr_host = 0, A_flt=0, display_result = True):
		'''
		subtract host flux

		INPUTS:
			flt:
			mag_host:
			magerr_host:
			A_flt: extinction in filter 'flt'
		'''
		obstime_new = []
		mag_new = []
		magerr_new = []

		flts_Vega =  ['U','B','V']
		lc_flt = self.__get_lc_flt(flt)

		for lcp in lc_flt:
			obstime = lcp['obstime']
			mag     = lcp['calmag']
			magerr  = lcp['calmagerr']

			print obstime, mag, magerr

			if mag > mag_host:
				continue

			if flt in flts_Vega:
				mag_tot_AB  = Vega_AB_mag_convertion(mag,flt,mode='Bessell',direction='Vega2AB')
	              		mag_host_AB = Vega_AB_mag_convertion(mag_host,flt,mode='Bessell',direction='Vega2AB')
			else:
				mag_tot_AB  = mag
				mag_host_AB = mag_host


			if flt == 'gp':
				flt_temp = 'g'
			elif flt == 'rp':
				flt_temp = 'r'
			elif flt == 'ip':
				flt_temp = 'i'
			elif flt == 'zp':
				flt_temp = 'z'
			elif flt == 'yp':
				flt_temp = 'y'
			else:
				flt_temp = flt

			lamb_flt,flux_tot  = mag2flux(mag_tot_AB - A_flt, flt_temp)
                        lamb_flt,flux_host = mag2flux(mag_host_AB - A_flt,flt_temp)
                        flux_transient = flux_tot - flux_host

			magerr_tot = magerr
                        fluxerr_transient   = np.sqrt(magerr_tot**2 + magerr_host**2)*1.086*flux_transient
			fluxupper_transient = flux_transient + fluxerr_transient

			mag_transient    = flux2mag(flux_transient,flt_temp,mode='wavelength',I_unit='cgs',wavelength_units='Angstroms')
			magerr_transient = mag_transient - flux2mag(fluxupper_transient,flt_temp,mode='wavelength',I_unit='cgs',wavelength_units='Angstroms')
			obstime_new.append(obstime)
			mag_new.append(mag_transient)
			magerr_new.append(magerr_transient)

		if display_result:
			plt.errorbar(obstime_new, mag_new, yerr=magerr_new, fmt='o')
			plt.gca().invert_yaxis()
			plt.xlabel('JD')
			plt.ylabel('mag')
			plt.title('%s band transient light curve'%flt)
                	plt.show()

		lc_hostfree_data = [obstime_new, mag_new, magerr_new]
		lc_hostfree_table = Table(data=lc_hostfree_data, names=['obstime', 'calmag', 'calmagerr'])

		return lc_hostfree_table


	def Iowa_prepare_template_image_wcs(self):
		'''
		The Iowa images have the non-standard header issue which will cause the problems when processing by the astropy.wcs function
		To get rid of the problem, only nessary keywords are kept when preparing Iowa template images
		'''
		for flt in self.templates.keys():
			tplimg = self.templates[flt]
			tempimg = self.Iowa_prepare_template_image_wcs_single(tplimg)


	def Iowa_prepare_template_image_wcs_single(self, tplimg, outimg=None):
		'''
		solve the non-standard issue of Iowa image
		'''
		inimg = os.path.join(self.raw_image_dir, tplimg)
		if outimg is None:
			outimg = os.path.join(self.template_dir, 'cal_%s'%tplimg)
		simplify_and_standardize_fits_header(inimg, outimg, out_header_only=False, overwrite=True)
		return outimg



	def wfcam_photometry_calibration_casu_zpt(self):
		'''
		For photometry on WFCAM images, the CASU data reduction pipeline gives the photometry zero point which is derived from calibration with 2MASS, details here http://www3.interscience.wiley.com/cgi-bin/fulltext?ID=122210794&PLACEBO=IE.pdf&mode=pdf

		'''

		for img in self.photometry_info.keys():
			if self.photometry_info[img]['airmass'] == 0:
				self.__renew_airmass_info(img)

			if self.photometry_info[img]['magzpt'] == 99.99:
				self.__get_magzpt_wfcam(img)

			instmag = self.photometry_info[img]['instmag']
			instmagerr = self.photometry_info[img]['instmagerr']
			magzpt = self.photometry_info[img]['magzpt']
			magzpterr = self.photometry_info[img]['magzpterr']

			exptime = self.photometry_info[img]['exptime']
			photmethod = self.photometry_method

			airmass = self.photometry_info[img]['airmass']

			mag, magerr = self.__wfcam_stdcal_casu_magzpt(instmag, instmagerr, magzpt, magzpterr, exptime, photmethod, airmass)

			self.photometry_info[img]['calmag'] = mag
			self.photometry_info[img]['calmagerr'] = magerr


	def __get_magzpt_wfcam(self, imgkey):
		'''
		Get the WFCAM photometry zeropoint from the image header which has been obtained by CASU data reduction pipeline
		'''
		input_image = os.path.join(self.modified_image_dir, imgkey)
		fitsinfo = get_fits_info(input_image, ['MAGZPT','MAGZRR'])
		magzpt = fitsinfo['MAGZPT']
		magzpterr = fitsinfo['MAGZRR']

		self.photometry_info[imgkey]['magzpt'] = magzpt
		self.photometry_info[imgkey]['magzpterr'] = magzpterr


	def __wfcam_stdcal_casu_magzpt(self, instmag, instmagerr, zpt, zpterr, exptime, photmethod, airmass):
		'''
		Apply the CASU photometry zeropoint to instrumental magnitude
		'''
		if photmethod == 'apphot':
			mag = zpt + instmag - 25 + 2.5*np.log10(exptime) - 0.05*(airmass-1)
		elif photmethod == 'psfphot':
			mag = zpt + instmag + 2.5*np.log10(exptime) - 0.05*(airmass-1)
		else:
			raise ValueError("photometry method %s not supported..."%photmethod)

		magerr = np.sqrt(instmagerr**2 + zpterr**2)

		return mag, magerr


	def explore_color_term_effect_in_calibrtion(self, flt, flt2):
		'''
		study the color term coefficients
		'''
		ps  = np.array([])
		eps = np.array([])

		for img in self.images.keys():
			if self.photometry_info[img]['flt'] == flt and self.photometry_info[img]['drop']==0:
				print img
				figoutfile = os.path.join(self.std_ref_dir, flt+'_'+img.split('.')[0]+'_cal.pdf')
				popt,perr = self.photcal_with_colorterm_single(img, flt, flt2, figoutfile=figoutfile)
				ps = np.append(ps,popt)
				eps = np.append(eps, perr)

		return ps, eps

	def photcal_with_colorterm_single(self, image_key, flt, flt2, sigclipfirst=True, sig1=6, sig2=3, figoutfile=None):
		'''
		This is created for the purpose of exploring the potential color term effect in calibration against to given standard stars

		INPUTS:
			image_key:
			flt:
			flt2: another flt to set the color
			sigclipfirst: do sigma clipping first before the calibration process
			sig1: sigma parameter for self.__sigma_clipping
			sig2: sigma parameter for self.__funcfit_remove_outlier
		'''
		imgkey = image_key.split('.')[0]

		refimg_abs = self.images[image_key]
		stds = self.__get_standard_reference_star_APASS(flt,wcs=True, single_use = True, refimg = refimg_abs)
		ref_list_file_raw = os.path.join(self.std_ref_dir, 'std_ref_whole_info_' + imgkey + '.txt')

		if flt in ['B','V']:
			magcol = flt+'mag'
		else:
			magcol = flt+'_mag'

		if flt2 in ['B','V']:
			magcol2 = flt2+'mag'
		else:
			magcol2 = flt2+'_mag'

		magerrcol = 'e_' + magcol
		magerrcol2 = 'e_' + magcol2

		wantcols = ['x','y',magcol, magerrcol, magcol2, magerrcol2]
		refdata = Table.read(ref_list_file_raw, format='ascii.csv')
		refdata = refdata[(refdata[magcol2]!=0)*(refdata[magcol]!=0)*(refdata[magcol]<50)*(refdata[magcol2]<50)]
		ref_list_file = os.path.join(self.std_ref_dir, 'std_ref_'+flt+flt2+'_'+imgkey + '.txt')
		self.__select_columns_from_table_to_ascii(refdata, wantcols, ref_list_file)

		match_output_file = os.path.join(self.std_ref_dir, imgkey +'_std_'+flt+flt2+'.match')
		data = self.__stdcal_link_std_obs_surrounding_search_single_image(image_key, flt, ref_list_file, stds_shift_x=0, stds_shift_y=0)

		mags_ref    = data[:,2]
		magserr_ref  = data[:,3]
		mags_input  = data[:,8]
		magserr_input= data[:,9]
		color = data[:,2] - data[:,4]
		colorerr =  np.sqrt(data[:,3]**2+data[:,5]**2)

		mags_offset = mags_ref - mags_input
		magserr_offset = np.sqrt(magserr_ref**2+magserr_input**2)



		magoffset0 = np.median(mags_offset)
		if self.cal_offset_funcfit_type == 'o1poly':
			def func(x,a,b):
				return a*x+b
			def yerr(x,aerr,berr):
				return np.sqrt((x*aerr)**2+berr**2)
			p0 = np.array([0,magoffset0])
		elif self.cal_offset_funcfit_type == 'constant':
			def func(x,c):
				return x*0+c
			def yerr(x,cerr):
				return x*0+cerr
			p0 = np.array([magoffset0])
		else:
			raise ValueError('%s for self.cal_offset_funcfit_type not supported yet'%self.cal_offset_funcfit_type)

		if sigclipfirst:
			survive_mask,discard_mask = self.__sigma_clipping(mags_offset,variability_method='mad', sig=sig1)
			color_rmd = color[discard_mask]
			mags_offset_rmd = mags_offset[discard_mask]
			magserr_offset_rmd = magserr_offset[discard_mask]
			color_svv = color[survive_mask]
			mags_offset_svv = mags_offset[survive_mask]
			magserr_offset_svv = magserr_offset[survive_mask]

		fitdata, popt, perr = self.__funcfit_remove_outlier(color_svv, mags_offset_svv, magserr_offset_svv,func, p0, nsig=sig2, rm_mode = 'all')
		offset = func(fitdata[:,0],*popt)

		cut = int(np.floor((len(fitdata)*0.688)))
		offset_dev = np.abs(fitdata[:,1]-offset)
		indice = np.argsort(offset_dev)[cut]
		offset_err = offset_dev[indice]

		if self.cal_plot:
			from matplotlib import gridspec
			fig = plt.figure(figsize=(9, 6))
			gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
			ax0 = plt.subplot(gs[0])
			plt.xlabel('X')
			plt.ylabel('magnitude offset (reference - input)')
			ax1 = plt.subplot(gs[1],sharey=ax0)

			ax0.errorbar(color_svv, mags_offset_svv, yerr=magserr_offset_svv, fmt='gx',label='raw')
			ax0.errorbar(fitdata[:,0],fitdata[:,1],yerr=fitdata[:,2],fmt='ro',label='fitdata')
			if sigclipfirst:
				ax0.errorbar(color_rmd, mags_offset_rmd, yerr=magserr_offset_rmd, fmt='bo',label='very raw', alpha=0.5)

			ax0.plot(fitdata[:,0], offset,'k')
			ax0.plot(fitdata[:,0], offset+offset_err,'b',alpha=0.4)
			ax0.plot(fitdata[:,0], offset-offset_err,'b',alpha=0.4)

			hist,edges = np.histogram(fitdata[:,1],bins=15)
			bin_centers = edges[:-1] + (edges[1:]-edges[:-1])/2.
			ax1.step(hist,bin_centers)

			if figoutfile is not None:
				plt.title(image_key)
				plt.savefig(figoutfile)
				plt.close()
			else:
				plt.show()

		if self.cal_check_stds:
			pylab_talk_to_ds9(color, mags_offset, data, refimg_abs, xcol=0,ycol=1, yerrplot=magserr_offset)


		return popt,perr


	def standard_calibration_direct(self, flts='all', std_catalog='apass', auto_shift_xy=True, verbose=0):
		'''
		directly calibrate the science image to the standard catalogue without the relative calibration within the same band

		See self.stdcal_single_image for detail
		'''
		if flts == 'all':
			if len(self.templates)==0:
				self.__find_template_imagekey(updatetable=1)
			fltstocal = self.templates.keys()
		elif isinstance(flts,list) or isinstance(flts,np.ndarray):
			fltstocal = flts
		else:
			raise IOError('Invalid input for flts')

		for flt in fltstocal:
			self.__dict2table()
			fltdata = self.result_table_newsaved[self.result_table_newsaved['flt']==flt]
			for img in fltdata['name']:
				if self.photometry_info[img]['drop'] == 0 and (self.photometry_info[img]['calmag'] == 99.99 or self.renew_stdcal_mag):
					self.stdcal_single_image(img, flt, std_catalog=std_catalog, link_stdref_obs_method='surrounding_search', auto_shift_xy=auto_shift_xy, verbose=verbose)



	def standard_calibration(self,flts='all', std_catalog = None,  tpl_has_wcs=False,all_tpl_astrometry=True,renew_tpl_astrometry=False,link_stdref_obs_method=None, std_obs_matched = False, stds_shift_x = 0, stds_shift_y = 0, updatetable=1):
		'''


		Inputs:
			flts: 'all', all filters whose template image found
			      list or np.ndarray, for example flts=['B','V','rp']

			std_catalog: 'apass', '2mass' or 'panstarrs'

			tpl_has_wcs: template images already have wcs information?
					If True, the wcs information will be applied to get the image positions
					of standard stars and use grmatch to get more precise match between
					stars on template image and stars from standard database;
				    	If False, astrometry will be needed to put standard stars on the template image;
					whether all templates need astrometry solution
					depend on 'all_tpl_astrometry'

			all_tpl_astrometry: get astrometry for all templates?
					If True, get astrometry for each template;
					If False, astrometry will be obtained for the best template, and other
					templates use that solution to get a rough image positions and
					 apply grmatch further.

			renew_tpl_astrometry:


			link_stdref_obs_method: the following options available
						'surrounding_search': this require the template image has correct wcs solution; and stars at the same position within given criteria are treated as iddentical stars
						'grmatch': match stdref stars and obs stars with 'grmatch' task; template image has valid wcs solution but not correct/precise
						'world_image_grmatch': the template image doesn't have wcs solution and difficult to obtain the wcs; then perform some simple world coordinate to image coordinate transformation based on knowledge of image size and CCD pixel scale and the [RA, Dec] of one pixel on the image, then perform grmatch the same as option 'grmatch'

			std_obs_matched:


		Please see self.standard_calibration_flt for details
		'''
		self.__find_template_imagekey(updatetable=updatetable)
		flts_obs = self.templates.keys()

		if flts == 'all':
			for flt in flts_obs:
				print "working on %s band"%flt
				self.standard_calibration_flt(flt, std_catalog = std_catalog, tpl_has_wcs=tpl_has_wcs,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry,std_obs_matched = std_obs_matched,link_stdref_obs_method=link_stdref_obs_method, stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y)
		elif isinstance(flts,list) or isinstance(flts,np.ndarray):
			for flt in flts:
				print "working on %s band"%flt
				self.standard_calibration_flt(flt, std_catalog = std_catalog, tpl_has_wcs=tpl_has_wcs,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry, std_obs_matched = std_obs_matched,link_stdref_obs_method = link_stdref_obs_method, stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y)
		else:
			raise IOError('Invalid input for flts')


	def standard_calibration_flt_simple(self, flt, m_zpt, em_zpt):
		'''
		calmag = instmag + m_zpt
		'''
		self.__dict2table()
		table_flt = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)
		for image_key in table_flt['name']:
			if self.photometry_info[image_key]['drop'] > 0:
				continue
			if (self.photometry_info[image_key]['calmag'] != 99.99) and (not self.renew_stdcal_mag):
				continue

			self.photometry_info[image_key]['calmag'] = self.photometry_info[image_key]['instmag'] + m_zpt
			calmagerr = np.sqrt(em_zpt**2 + self.photometry_info[image_key]['instmagerr']**2)
			self.photometry_info[image_key]['calmagerr'] = calmagerr



	def standard_calibration_flt(self,flt, std_catalog = None,  tpl_has_wcs=False,all_tpl_astrometry=True,renew_tpl_astrometry=False,link_stdref_obs_method = None,std_obs_matched = False, stds_shift_x =0, stds_shift_y =0, updatetable=1):
		'''
		standard calibration for image with filter 'flt'
		The current philosophy here is to perform calibration for the reference image and apply the calibration result to all images

		check the self.photometry_info[template_imgkey]['calmag'] first. If the value is not 99.99, then skip the calibration for reference image

		Inputs:
			flt:

			tpl_has_wcs: same as self.standard_calibration
			all_tpl_astrometry: same as self.standard_calibration
			renew_tpl_astrometry: same as self.standard_calibration


		Involved functions:
			self.__stdcal_template_prepare_tplimg_flt
			self.__stdcal_template_prepare_stdfile_flt
			self.__stdcal_template_flt

		'''
		apass_flts = ['B','V','R','I','gp','rp','ip']
		twomass_flts = ['J', 'H', 'K']
		self.__find_template_imagekey(updatetable=updatetable)

		flt_tplkey = self.templates[flt]
		if self.photometry_info[flt_tplkey]['calmag'] == 99.99:
			if not std_obs_matched:
				if link_stdref_obs_method != 'world_image_grmatch':
					self.__stdcal_template_prepare_tplimg_flt(flt,tpl_has_wcs=tpl_has_wcs,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry)
				if std_catalog is not None:
					which_stdcatalog = std_catalog
				else:
					if flt in apass_flts:
						which_stdcatalog = 'apass'
					elif flt in twomass_flts:
						which_stdcatalog = '2mass'
					else:
						raise ValueError('%s band photometry calibration is not supported'%flt)

				self.__stdcal_template_prepare_stdfile_flt(flt,link_stdref_obs_method = link_stdref_obs_method, which_catalog=which_stdcatalog)
				if link_stdref_obs_method is None:
					if not tpl_has_wcs and all_tpl_astrometry:
						link_stdref_obs_method = 'surrounding_search'
					else:
						link_stdref_obs_method = 'grmatch'

			self.__stdcal_template_flt(flt, link_stdref_obs_method = link_stdref_obs_method, std_obs_matched = std_obs_matched, stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y )
		elif self.renew_stdcal_mag:
			if not std_obs_matched:
				if link_stdref_obs_method != 'world_image_grmatch':
					self.__stdcal_template_prepare_tplimg_flt(flt,tpl_has_wcs=tpl_has_wcs,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry)
				if std_catalog is not None:
					which_stdcatalog = std_catalog
				else:
					if flt in apass_flts:
						which_stdcatalog = 'apass'
					elif flt in twomass_flts:
						which_stdcatalog = '2mass'
					else:
						raise ValueError('%s band photometry calibration is not supported'%flt)

				self.__stdcal_template_prepare_stdfile_flt(flt,link_stdref_obs_method = link_stdref_obs_method, which_catalog=which_stdcatalog)
				if link_stdref_obs_method is None:
					if not tpl_has_wcs and all_tpl_astrometry:
						link_stdref_obs_method = 'surrounding_search'
					else:
						link_stdref_obs_method = 'grmatch'

			self.__stdcal_template_flt(flt, link_stdref_obs_method = link_stdref_obs_method, std_obs_matched = std_obs_matched , stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y )

		else:
			print "old standard calibration result for images in filter %s will be used"%flt

		if updatetable:
			self.__dict2table()
		table_flt = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)

		for image_key in table_flt['name']:
			tpl_imgkey = self.templates[flt]
			if self.photometry_info[image_key]['template'] == 1:
				continue
			if self.photometry_info[image_key]['drop'] > 0:
				continue
			if (self.photometry_info[image_key]['calmag'] != 99.99) and (not self.renew_stdcal_mag):
				continue

			offset = self.photometry_info[tpl_imgkey]['calmag'] - self.photometry_info[tpl_imgkey]['relmag']
			self.photometry_info[image_key]['calmag'] = self.photometry_info[image_key]['relmag'] + offset

			instmagerr = self.photometry_info[image_key]['instmagerr']
			relmagerr = self.photometry_info[image_key]['relmagerr']
			stdmagerr = np.sqrt(self.photometry_info[tpl_imgkey]['calmagerr']**2 - self.photometry_info[tpl_imgkey]['instmagerr']**2)

			#Improvement needed!!!
			if self.photometry_info[tpl_imgkey]['drop'] >0 or stdmagerr==99.99:
				stdmagerr = 2*relmagerr

			calmagerr = np.sqrt(instmagerr**2 + relmagerr**2 + stdmagerr**2)
			self.photometry_info[image_key]['calmagerr'] = calmagerr

	def __get_astrometry_solution_grtrans(self, ref_catalog_input, ref_catalog_output, input_source_list, racol, deccol, magcol, raref=None, decref=None, verbose=0):
		'''
		get astrometry with programs grmatch and grtrans from fitsh
		INPUTS:
			ref_catalog_input: the input reference star catalog with ra,dec columns
			ref_catalog_output: the output reference star catalog with image x,y columns appended
			input_source_list: formated as x,y,mag,magerr
			racol: RA column index in input_catalog
			deccol: Dec column index in input_catalog
			magcol: magnitude column index in input_catalog
			raref: the reference RA for sky projection to image plane; around the center of FOV
			decref: the reference Dec for sky projection to image plane; around the center of FOV
		'''
		if raref is None:
			raref = self.sn_ra_world_deg
		if decref is None:
			decref = self.sn_dec_world_deg
		if self.pixscale is None:
			self.pixscale = 1.0
		scale = 1.0/(self.pixscale/3600.0/180*np.pi)
		grtrans = os.path.join(self.base.fitsh_dir, 'grtrans')
		temp = np.loadtxt(ref_catalog_input)
		M,N = temp.shape
		xrefcol = N+1
		yrefcol = N+2
		grtranscmd1 = "%s --input %s --wcs tan,scale=%s,ra=%s,dec=%s --col-radec %s,%s --col-out %s,%s --output %s"%(ref_catalog,scale,raref,decref,racol, deccol,xrefcol,yrefcol,ref_catalog_output)
		try:
			if verbose:
				print grtranscmd1
			os.system(grtranscmd1)
		except:
			raise OSError("grtranscmd1 failure, pls check")

		self.fitsh_grmatch_pointmatch_pars['--col-ref'] = '%s,%s'%(xrefcol,yrefcol) #The index of the first column is always 1
		#self.fitsh_grmatch_pointmatch_pars['--col-inp'] = '1,2'
		#self.fitsh_grmatch_pointmatch_pars['--order'] = 1 #If the order is A, >= (A+1)*(A+2)/2 valid points are needed to fit the transformation
		#self.fitsh_grmatch_pointmatch_pars['--max-distance'] = 1
		#self.fitsh_grmatch_pointmatch_pars['--triangulation'] = 'auto,mixed,maxnumber=200'
		self.fitsh_grmatch_pointmatch_pars['--col-ref-ordering'] = -magcol #negative sign indicates ascending, small values first
		#self.fitsh_grmatch_pointmatch_pars['--col-inp-ordering'] = -3

		self.__fitsh_grmatch(ref_catalog_output, input_source_list, match_output, trans_output, mode ='file')



	def __get_astrometry_solution_API_client(self,img,wcs_out, time_out_max = 300, apikey ='aubtqramhfmriufp' ):
		'''
		Get astrometry for image 'img' with nova.astrometry.net API and the wcs result will be saved in 'wcs_out'

		Please see 'astrometry_API/client.py' for details
		'''

		api_dir = os.path.join(self.base.base_dir, 'extern/nova_astrometry_api')
		API_client = os.path.join(api_dir, "client.py")
		if not os.path.exists(API_client):
			raise IOError("API client for astrometry not available")

		img_key = os.path.basename(img)
		print "working on astrometry for %s"%img

		snra = self.sn_ra_world_deg
		sndec = self.sn_dec_world_deg
		imgfr = 0.006  #image field radius 21.6 armin

		command_line = "python %s --server http://nova.astrometry.net/api/ --apikey %s --wcs %s --upload %s --ra %s --dec %s --radius %s -p >/dev/null"%(API_client, apikey,  wcs_out, img, snra, sndec, imgfr)
		command = Command(command_line)
		failure_signal = command.run(timeout=time_out_max)
		#failure_signal  non-zero value means failure

		if not failure_signal:
			success = True
		else:
			success = False

		return success




	def __stdcal_template_prepare_tplimg(self,flts='all',tpl_has_wcs=False,all_tpl_astrometry=True,renew_tpl_astrometry=False, updatetable=1):
		'''
		Prepare template images
		'''
		self.__find_template_imagekey(updatetable=updatetable)
		flts_obs = self.templates.keys()

		if flts == 'all':
			for flt in flts_obs:
				self.__stdcal_template_prepare_tplimg_flt(flt,tpl_has_wcs=tpl_has_wcs,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry)
		elif isinstance(flts,list) or isinstance(flts,np.ndarray):
			for flt in flts:
				self.__stdcal_template_prepare_tplimg_flt(flt,tpl_has_wcs=tpl_has_astrometry,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry)
		else:
			raise IOError('Invalid input for flts')




	def __stdcal_template_prepare_tplimg_flt(self,flt,tpl_has_wcs=False,all_tpl_astrometry=True,renew_tpl_astrometry=False):
		'''
		prepare reference image which will be used to link the calibration reference stars
		from both standard datebase and new measurements.
		Inputs:
			refer to 'standard_calibration' for explanation
		'''
		maxtime_on_astrometry = self.maxtime_on_astrometry
		tpl_key = self.templates[flt]

        	if self.current_telescope == 'WFCAM':
			tpl_img_orig = os.path.join(self.modified_image_dir, tpl_key)
		else:
	        	tpl_img_orig = os.path.join(self.raw_image_dir,tpl_key)
	        tpl_img_for_stdcal = self.templates_after_astrometry[flt]

		finish_code = 1
                if tpl_has_wcs:
                	if not os.path.exists(tpl_img_for_stdcal) or renew_tpl_astrometry:
				shutil.copy(tpl_img_orig, tpl_img_for_stdcal)

                elif all_tpl_astrometry:
                	if not os.path.exists(tpl_img_for_stdcal) or renew_tpl_astrometry:
	        		success = self.__get_astrometry_solution_API_client(tpl_img_orig,tpl_img_for_stdcal, time_out_max = maxtime_on_astrometry)
				if not success:
					finish_code = 0
		else:
			if self.reference_flt is not None:
				ref_flt = self.reference_flt
				tpl_ref_flt_key  = self.templates[ref_flt]
				tpl_img_ref_flt_orig = os.path.join(self.raw_image_dir,tpl_ref_flt_key)
				tpl_img_for_stdcal_ref_flt = self.templates_after_astrometry[ref_flt]

				if not os.path.exists(tpl_img_for_stdcal_ref_flt) or renew_tpl_astrometry:
					success = self.__get_astrometry_solution_API_client(tpl_img_ref_flt_orig,tpl_img_for_stdcal_ref_flt, time_out_max = maxtime_on_astrometry)
					if not success:
						finish_code = 0
				if flt != ref_flt:
					shutil.copy(tpl_img_for_stdcal_ref_flt,tpl_img_for_stdcal)
			else:
				raise KeyError("Specify self.reference_flt please!")

		return finish_code




	def __stdcal_template_prepare_stdfile_flt(self,flt,link_stdref_obs_method=None, which_catalog = 'apass'):
		'''
		Get the reference stars for magnitude calibration

		INPUTs:
			flt:
			link_stdref_obs_method: see self.standard_calibration
			which_catalog: 'apass' or '2mass'

		'''

		if link_stdref_obs_method == "world_image_grmatch":
			wcs = False
		else:
			wcs = True

		if which_catalog == 'apass':
			self.__get_standard_reference_star_APASS(flt, wcs=wcs)
		elif which_catalog == '2mass':
			self.__get_standard_reference_star_2MASS(flt, wcs=wcs)
		elif which_catalog == 'panstarrs':
			magtype = self.panstarrs_mag_photometry_method
			self.__get_standard_reference_star_PanSTARRS(flt, wcs=wcs, photmethod = magtype)
		else:
			raise ValueError('catalog %s is not supported...'%which_catalog)



	def __stdcal_template_link_std_obs_flt(self,flt, method = 'grmatch', stds_shift_x=0, stds_shift_y = 0):
		'''
		get the one-one mapping between standards and new measurements

		Inputs:
			methods: available options currently available:
				grmatch:
				selection_on_ds9:
				external:
				surrounding_search:
		'''

		if method == 'grmatch':
			self.__stdcal_template_link_std_obs_grmatch_flt(flt)
		elif method == 'world_image_grmatch':
			self.__stdcal_template_link_std_obs_grmatch_flt(flt)
		elif method == 'surrounding_search':
			self.__stdcal_template_link_std_obs_surrounding_search_flt(flt, stds_shift_x=stds_shift_x, stds_shift_y= stds_shift_y)
		elif method == 'selection_on_ds9':
			self.__stdcal_template_link_std_obs_selection_on_ds9_flt(flt)
		elif method == 'external':
			self.__stdcal_template_link_std_obs_selection_external_flt(flt)
		else:
			raise IOError("Invalid input for method, available options are grmatch, selection_on_ds9 and external")


	def __stdcal_link_std_obs_single_image(self,imgkey,flt,method = 'grmatch', single_use = False, stds_shift_x=0, stds_shift_y =0, external_refimg=None):
		'''
		match sources on given image with standards

		INPUTS:
			image:
			flt:
			method: 'grmatch', 'surrounding_search', 'xytran', 'selection_on_ds9', 'external'
			single_use:
			stds_shift_x: shift standards in x direction by this value if method is 'surrounding_search'
			stds_shift_y: shift standards in y direction by this value if method is 'surrounding_search'
		'''
		imgkey_s = imgkey.split('.')[0]
		if not single_use:
      		  		ref_list_file = os.path.join(self.std_ref_dir, self.flt_std_ref_stars_file[flt])	#x y mag magerr
        	else:
				ref_list_file = os.path.join(self.std_ref_dir, 'std_ref_' + imgkey_s + '.txt')

		if method == 'grmatch':
			match_ret = self.__stdcal_template_link_std_obs_grmatch_single_image(imgkey,flt, single_use = single_use)
		elif method == 'surrounding_search':
			match_ret = self.__stdcal_link_std_obs_surrounding_search_single_image(imgkey,flt, ref_list_file, stds_shift_x=stds_shift_x, stds_shift_y= stds_shift_y)
		elif method == 'xytran':
			match_ret = self.__stdcal_template_link_std_obs_xytrans_single_image(imgkey, flt, single_use = single_use, verbose = 1, external_refimg=external_refimg)
		elif method == 'selection_on_ds9':
			match_ret = self.__stdcal_template_link_std_obs_selection_on_ds9_single_image(imgkey,flt)
		elif method == 'external':
			match_ret = self.__stdcal_template_link_std_obs_selection_external_single_image(imgkey,flt)
		else:
			raise IOError("Invalid input for method, available options are grmatch, surrouding_search, xytran, selection_on_ds9 and external")

		return match_ret

	def __stdcal_template_link_std_obs_surrounding_search_flt(self,flt, stds_shift_x=0, stds_shift_y = 0):
		'''
		Get one-one mapping between standard reference stars and sources detected in template image
		'''
		tpl_imgkey = self.templates[flt]
      		ref_list_file = os.path.join(self.std_ref_dir, self.flt_std_ref_stars_file[flt])	#x y mag magerr
		match_ret = self.__stdcal_link_std_obs_surrounding_search_single_image(tpl_imgkey,flt, ref_list_file)

		return match_ret


	def __stdcal_template_link_std_obs_xytrans_single_image(self, imgkey,flt, single_use = False, verbose = False, external_refimg=None):
		'''
		First manually select 3 pair of sources from std catalog and new measurements; then use the selected pairs to do transformation for the whole catalog; finally do 		  the normal surrounding search
		'''
		imgkey_s = imgkey.split('.')[0]

		if not single_use:
      		  	ref_list_file = os.path.join(self.std_ref_dir, self.flt_std_ref_stars_file[flt])	#x y mag magerr
			stdxymags_regionfile =  os.path.join(self.std_ref_dir, flt + '.reg')
        	else:
			ref_list_file = os.path.join(self.std_ref_dir, 'std_ref_' + imgkey_s + '.txt')
			stdxymags_regionfile = os.path.join(self.std_ref_dir, imgkey_s + '.reg')

        	if self.photometry_method == 'apphot':
        		input_list_file = os.path.join(self.aperture_photometry_dir, imgkey_s + self.apphot_ret_file_suffix)
        		#column names in ref_list and input_list: xpos, ypos, mag,mag_err
			inxymags_regionfile = os.path.join(self.aperture_photometry_dir, imgkey_s + '.reg')
        	elif self.photometry_method == 'psfphot':
        		input_list_file = os.path.join(self.psf_photometry_dir, imgkey_s + self.psfphot_ret_file_suffix)
        		inxymags_regionfile = os.path.join(self.psf_photometry_dir, imgkey_s + '.reg')
        	else:
        		raise IOError("Invalid input for photometry_method")

		inxymags = np.loadtxt(input_list_file)
		stdxymags = np.loadtxt(ref_list_file)
		fitsimg = self.images[imgkey]


		match_listfile = os.path.join(self.std_ref_dir, imgkey_s + '_tpl_stdref_xytrans.input')
        	match_database = os.path.join(self.std_ref_dir, imgkey_s + '_tpl_stdref_xytrans.coef')

		if os.path.exists(match_listfile) and os.path.exists(match_database) and (not self.renew_stdcal_xytran):
			print "%s and %s exist; if you want to renew, modify self.renew_stdcal_xytran"%(match_listfile,match_database)
		else:
			if not os.path.exists(stdxymags_regionfile):
				create_ds9_region_file(stdxymags, stdxymags_regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='green', width=1, load_text=False,  text_uoffset=25, text_loffset=25, textcol1=2,textcol2=3)
			if external_refimg is not None:
				dispimg = external_refimg
			else:
				dispimg = fitsimg

			select_num = self.stdcal_xytran_stdnum

			stdlist = self.__input_xys_through_ds9_get_mags(dispimg, stdxymags, regionfile=stdxymags_regionfile, wantnum=select_num, ds9name=None, verbose=1)#output table colnames: num, x, y, mag, magerr
			if not os.path.exists(inxymags_regionfile):
				create_ds9_region_file(inxymags, inxymags_regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='green', width=1, load_text=False,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)
			inlist = self.__input_xys_through_ds9_get_mags(fitsimg, inxymags, regionfile=inxymags_regionfile, wantnum=select_num, ds9name=None, verbose=1)

			stdlist.write(os.path.join(self.std_ref_dir, 'std_ref_pickups_for_%s.txt'%imgkey_s), format='ascii.fixed_width')
			inlist.write(os.path.join(self.std_ref_dir, 'obs_pickups_on_%s.txt'%imgkey_s), format='ascii.fixed_width')

			reference_list = Table()
			reference_list.add_columns([stdlist['x'], stdlist['y']])
			reference_list.rename_column('x','xstd')
			reference_list.rename_column('y','ystd')
			reference_list.add_columns([inlist['x'], inlist['y']])
			reference_list.write(match_listfile, format='ascii.no_header')


		geomap_iraf(match_listfile, match_database)

		tg_xys = inxymags[:,0:2]
		ref_xys_beforetrans = stdxymags[:,0:2]
		ref_xys_beforetrans_file = os.path.join(self.std_ref_dir, 'std_ref_'+imgkey_s + '_xy_beforetran.txt')
		np.savetxt(ref_xys_beforetrans_file, ref_xys_beforetrans)
		ref_xys_aftertran_file   = os.path.join(self.std_ref_dir, 'std_ref_'+imgkey_s + '_xy_aftertran.txt')
		geoxytran_iraf(ref_xys_beforetrans_file, ref_xys_aftertran_file, match_database, match_listfile)

		ref_xys = np.loadtxt(ref_xys_aftertran_file)
		criteria = self.stdref_tpl_match_criteria
		target_stars,indexmat,ref_stars,indexref = self.__crude_match(tg_xys,ref_xys,criteria)

		if verbose:
			print indexmat

		match_ret = np.hstack((stdxymags[indexref], inxymags[indexmat]))

        	std_ref_and_obs_measurement_stars_file = imgkey_s +'_std_'+flt+'.match'
        	match_output_file  = os.path.join(self.std_ref_dir, std_ref_and_obs_measurement_stars_file)
		np.savetxt(match_output_file,match_ret)

		match_result_plot = self.stdref_tpl_match_result_display

		if match_result_plot:
			tg_xys_matched  = tg_xys[indexmat]
			ref_xys_matched = ref_xys[indexref]
			plt.plot(tg_xys_matched[:,0]-ref_xys_matched[:,0],tg_xys_matched[:,1]-ref_xys_matched[:,1],'o')
			plt.show()

		return match_ret



	def __stdcal_link_std_obs_surrounding_search_single_image(self, imgkey, flt, ref_list_file, stds_shift_x=0, stds_shift_y=0, match_output_file=None, verbose = False):
		'''
		match two catalog according to the (x,y) coordinates.
		The two catalogs (std and obs) are obtained for image with imgkey.
		The first two columns of both catalogs are image coordinate x and y.

		INPUTS:
			imgkey:
			flt:
			single_use: if True, 'the ref_list_file' will be slightly different
			stdwhole:if True, the std file contain whole package of catalog information
			stds_shift_x:
			stds_shift_y:
		'''

		imgkey_s = imgkey.split('.')[0]
		if match_output_file is None:
        		match_output_file  = os.path.join(self.std_ref_dir, imgkey_s +'_std_'+flt+'.match')

        	if self.photometry_method == 'apphot':
        		input_list_file = os.path.join(self.aperture_photometry_dir, imgkey_s + self.apphot_ret_file_suffix)
        		#column names in ref_list and input_list: xpos, ypos, mag,mag_err
        	elif self.photometry_method == 'psfphot':
        		input_list_file = os.path.join(self.psf_photometry_dir, imgkey_s + self.psfphot_ret_file_suffix)
        	else:
        		raise IOError("Invalid input for photometry_method")
		if verbose:
			print input_list_file

		input_list = np.loadtxt(input_list_file)
		tg_xys = input_list[:,0:2]
		ref_list = np.loadtxt(ref_list_file)
		ref_xys = ref_list[:,0:2]
		ref_xys[:,0] = ref_xys[:,0] + stds_shift_x
		ref_xys[:,1] = ref_xys[:,1] + stds_shift_y

		criteria = self.stdref_tpl_match_criteria
		target_stars,indexmat,ref_stars,indexref = self.__crude_match(tg_xys,ref_xys,criteria)

		if verbose:
			print indexmat

		match_ret = np.hstack((ref_list[indexref],input_list[indexmat]))
		np.savetxt(match_output_file,match_ret)

		match_result_plot = self.stdref_tpl_match_result_display
		if match_result_plot:
			tg_xys_matched  = tg_xys[indexmat]
			ref_xys_matched = ref_xys[indexref]
			fig = plt.figure(figsize=(6,6))
			plt.plot(tg_xys_matched[:,0]-ref_xys_matched[:,0],tg_xys_matched[:,1]-ref_xys_matched[:,1],'o')
			plt.show()

		return match_ret

	def __stdcal_template_link_std_obs_selection_on_ds9_flt(self,flt):
		tpl_imgkey = self.templates[flt]

		self.__stdcal_template_link_std_obs_selection_on_ds9_single_image(tpl_imgkey,flt)

		raise IOError("under construction...")




	def get_template_reference_mags_flt(self, flt,  which_dir='raw_image', renew_refmags=False):
		'''
		get the template mags table. The output format is: id, x, y, mag, magerr


		INPUTS:
			flt:
			which_dir:
			renew_refmags:
		'''

		if flt not in self.templates.keys():
			raise ValueError("no template image avaiable for filter %s"%flt)

		image_key = self.templates[flt]
		tpl_imgkey_s = image_key.split('.')[0]

		if self.photometry_method is None:
			self.photometry_method = raw_input("please input the photometry method used(apphot/psfphot):")

		if self.photometry_method   == 'apphot':
			ref_list_filename      = os.path.join(self.aperture_photometry_dir, tpl_imgkey_s + self.apphot_ret_file_suffix)
		elif self.photometry_method == 'psfphot':
			ref_list_filename      = os.path.join(self.psf_photometry_dir, tpl_imgkey_s + self.psfphot_ret_file_suffix)
		else:
			raise IOError("Invalid input for photometry_method")


		if not os.path.exists(ref_list_filename):
			raise IOError("PSF photometry result file %s not exists"%ref_list_filename)

		psfphot_ret = np.loadtxt(ref_list_filename)
		xys = psfphot_ret[:,0:2]
		mags = psfphot_ret[:,2]
		magerrs = psfphot_ret[:,3]

		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		regionfile = os.path.join(self.template_dir, "%s_%s_xymag.reg"%(flt, tpl_imgkey_s))
		if (not os.path.exists(regionfile)) or renew_refmags:
			create_ds9_region_file(psfphot_ret, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='red', width=1, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)


		refmag_table = self.__input_xys_through_ds9_get_mags(input_image, psfphot_ret, regionfile=regionfile)

		refmag_filename = os.path.join(self.template_dir, "%s_%s_idxymag.txt"%(flt, tpl_imgkey_s))
		refmag_table.write(refmag_filename, format='ascii.fixed_width')

		return refmag_table



	def __input_xys_through_ds9_get_mags(self, fitsimg, xysmags, regionfile=None, wantnum=None, newds9=True, verbose=1):
		'''
		pick up stars from ds9 display and extract mags from file xysmagfile

		INPUTS:
			fitsimg:
			xysmags: format x,y,mag,magerr
		'''
		if verbose:
			print "Now start clicking on ds9 to pick up sources"
			print "You can modify self.stdref_tpl_match_criteria to define the searching cone radius"


		xys = xysmags[:,0:2]
		print xys

		i = 0

		num_flags = []
		xs = []
		ys = []
		mags = []
		magerrs = []

		while True:
			if wantnum is not None:
				if i>(wantnum-1):
					break

			num_default = i
			if wantnum is None:
				num_flag = raw_input("please input the num flag for the star (default:%s )"%num_default) or num_default
			else:
				num_flag = num_default

			num_int = int(num_flag)
			i = num_int + 1

			if num_int < 0:
				break

			xy = self.__get_xy_on_image_from_ds9(input_image=fitsimg, regionfile=regionfile, newds9=newds9)
			print "xy from mouse pick up:", xy

			x = xy[0] - 1
			y = xy[1] - 1

			xy_target = [x,y]

			criteria  = self.stdref_tpl_match_criteria
			yesfound, xy_match,index = self.__find_corresponding(xys, xy_target, criteria)

			if (not yesfound):
				print "no object found within criteeia..."
				mag = 99.99
				magerr = 99.99
			elif len(index)>1:
                        	print 'more than one objects fall into the criteria region of reference star. This object is droped!'
				mag = 99.99
				magerr = 99.99
			else:
				print "the corresponding object found at (%s)"%xy_match
				mag = xysmags[index,2][0]
				magerr = xysmags[index,3][0]

			num_flags.append(num_int)
			xs.append(x)
			ys.append(y)

			mags.append(mag)
			magerrs.append(magerr)

		num_flags_col = Column(name='num',data=num_flags)
		xs_col = Column(name='x', data = xs)
		ys_col = Column(name='y', data = ys)
		mags_col = Column(name = 'mag', data = mags)
		magerrs_col = Column(name = 'magerr',data = magerrs)

		out_table = Table()
		out_table.add_columns([num_flags_col, xs_col, ys_col, mags_col, magerrs_col])

		print out_table

		return out_table



	def __stdcal_template_link_std_obs_external_flt(self,flt):
		raise IOError("under construction...")


	def __stdcal_template_link_std_obs_grmatch_flt(self,flt):

		tpl_imgkey = self.templates[flt]
		match_ret = self.__stdcal_template_link_std_obs_grmatch_single_image(tpl_imgkey,flt)

		return match_ret

	def __stdcal_template_link_std_obs_grmatch_single_image(self, imgkey,flt,single_use= False):
		'''
		match standard stars with sources detected on given image with grmatch

		INPUTS:
			imgkey:
			flt:
		'''
		imgkey_s = imgkey.split('.')[0]

		if not single_use:
			ref_list = os.path.join(self.std_ref_dir, self.flt_std_ref_stars_file[flt])	#x y mag magerr
		else:
			ref_list = os.path.join(self.std_ref_dir, 'std_ref_' + imgkey_s + '.txt')

	        std_ref_and_obs_measurement_stars_file  = imgkey_s +'_std_'+flt+'.match'
		match_output  = os.path.join(self.std_ref_dir, std_ref_and_obs_measurement_stars_file)
		trans_fitting = os.path.join(self.std_ref_dir, imgkey_s + '_tpl_stdref.coef')


		if self.photometry_method == 'apphot':
			input_list = os.path.join(self.aperture_photometry_dir, imgkey_s + self.apphot_ret_file_suffix)
			#column names in ref_list and input_list: xpos, ypos, mag,mag_err

		elif self.photometry_method == 'psfphot':
			input_list = os.path.join(self.psf_photometry_dir, imgkey_s + self.psfphot_ret_file_suffix)

		else:
			raise IOError("Invalid input for photometry_method")

		if (not os.path.exists(match_output)) or self.renew_std_ref_match:
			self.__delete_file_if_exist(match_output)
			self.__delete_file_if_exist(trans_fitting)
			matched_table = self.__fitsh_grmatch(ref_list, input_list, match_output, trans_fitting)

			match_result_plot = self.stdref_tpl_match_result_display
			if match_result_plot:
				if not os.path.exists(match_output):
					raise IOError("grmatch failure...")
				trans_output = os.path.join(self.std_ref_dir, imgkey_s + '_tpl_stdref.trans')
				self.__fitsh_grtrans(match_output, trans_output, trans_fitting)
				trans_result = np.loadtxt(trans_output)
				ref_xys = trans_result[:,[0,1]]
				input_xys  = trans_result[:,[4,5]]
				plt.plot(ref_xys[:,0]-input_xys[:,0],ref_xys[:,1]-input_xys[:,1],'o')
				plt.show()
		else:
			print "std and obs match result exists %s"%match_output
			matched_table = np.loadtxt(match_output)

		return matched_table


	def __stdcal_template_flt(self,flt,link_stdref_obs_method='grmatch', std_obs_matched = False, stds_shift_x = 0, stds_shift_y = 0):
		'''
		Find one-one corresponding between standard reference stars and measurements in template image
		perform relative calibration to get the mag offset between standards and new measurements.

		Main functions involved:
			self.__stdcal_template_link_std_obs_flt
			self.__relative_calibration_single

		Inputs:
			flt:
			link_stdref_obs_method: See self.__stdcal_template_link_std_obs_flt for allowed options
			std_obs_matched: if True, then one mapping file containing std mags and new measurements is given
					on which the calibration will be based

		'''

		self.__get_stdref_filenames()

		tpl_imgkey = self.templates[flt]
		tpl_key = tpl_imgkey.split('.')[0]

		if self.photometry_info[tpl_imgkey]['calmag'] != 99.99:
			if not self.renew_stdcal_mag:
				return

		if std_obs_matched:
			self.get_std_obs_match_filenames()
			std_obs_match = self.stdobs_match_files[tpl_imgkey]
			if not os.path.exists(std_obs_match):
				raise IOError("the std_obs_match file %s not exist"%std_obs_match)
			matched_ret = np.loadtxt(std_obs_match)
			M,N = matched_ret.shape
			if N<5:
				mags_ref_matched   = matched_ret[:,[0,1]]
				mags_input_matched = matched_ret[:,[2,3]]
			else:
				mags_ref_matched   = matched_ret[:,[2,3]]
				mags_input_matched = matched_ret[:,[6,7]]

		else:
			print link_stdref_obs_method
			self.__stdcal_template_link_std_obs_flt(flt,method=link_stdref_obs_method, stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y )
			std_obs_match = os.path.join(self.std_ref_dir,self.flt_std_ref_and_obs_measurement_stars_file[flt])
			matched_ret   = np.loadtxt(std_obs_match)

			if self.reference_mag_min is not None:
				mag_faint_cut = self.reference_mag_min
				matched_ret = matched_ret[matched_ret[:,2]<mag_faint_cut,:]

			if self.reference_mag_max is not None:
				mag_bright_cut  = self.reference_mag_max
				matched_ret = matched_ret[matched_ret[:,2]>mag_bright_cut,:]

			if self.reference_magerr_cut is not None:
				magerr_cut_ref = self.reference_magerr_cut
				matched_ret = matched_ret[matched_ret[:,3]<magerr_cut_ref,:]

			if self.input_magerr_cut is not None:
				magerr_cut_input = self.input_magerr_cut
				matched_ret = matched_ret[matched_ret[:,7]<magerr_cut_input,:]

			ref_matched   = matched_ret[:,[0,1,2,3]]
			input_matched = matched_ret[:,[4,5,6,7]]
			mags_ref_matched = ref_matched[:,[2,3]]
			mags_input_matched = input_matched[:,[2,3]]

		tpl_instmag = self.photometry_info[tpl_imgkey]['instmag']
		tpl_instmagerr = self.photometry_info[tpl_imgkey]['instmagerr']
		cal_mag,cal_magerr, temp_future_use = self.__relative_calibration_single(mags_ref_matched,mags_input_matched,tpl_instmag)

		self.photometry_info[tpl_imgkey]['calmag'] = cal_mag
		self.photometry_info[tpl_imgkey]['calmagerr'] =  np.sqrt( cal_magerr**2 + tpl_instmagerr**2 )




	def stdcal_single_image(self, imgkey,flt, instmag =None, instmagerr=None, link_stdref_obs_method='grmatch', std_catalog = 'apass', panstarrs_photmethod = 'PSF',  std_world2img_refimage=None,  stds_shift_x=0, stds_shift_y=0, std_obs_matched = False, verbose = 0, update_stdmag=True, auto_shift_xy = False):
		'''
		This function directly calibration given image to standard stars

		Find one-one corresponding between standard reference stars and measurements in given image
		perform relative calibration to get the mag offset between standards and new measurements.


		Main functions involved:
			self.__stdcal_template_link_std_obs_flt
			self.__relative_calibration_single

		Inputs:
			flt:
			link_stdref_obs_method: See self.__stdcal_template_link_std_obs_flt for allowed options
			std_obs_matched: if True, then one mapping file containing std mags and new measurements is given
					on which the calibration will be based

			std_catalog: 'apass', 'panstarrs', '2mass', or from private photometry result for example 'selfphot,LCOGT,001,psfphot' or 'selfphot,LCOGT,001,010,psfphot'
				    where LCOGT is telescope code, 001 is image index, and PSF is photometry method for which apphot is another option
			panstarrs_photmethod: 'PSF' or 'AP'

			stds_shift_x:
			stds_shift_y:

			std_world2img_refimage: the image containing WCS information acoording which transform the standard stars from world coordinate to image coordinate


		'''

		print imgkey
		key_s = imgkey.split('.')[0]

		if self.photometry_info[imgkey]['calmag'] != 99.99:
			if not self.renew_stdcal_mag:
				if verbose>0:
					print "calibration result exists for %s"%imgkey
				return
			else:
				if verbose > 0:
					print "renew calibration result for %s"%imgkey

		if std_obs_matched:
			self.get_std_obs_match_filenames()
			std_obs_match = self.stdobs_match_files[imgkey]
			matched_ret = np.loadtxt(std_obs_match)
			M,N = matched_ret.shape
			if N<5:  #only magnitude pairs
				mags_ref_matched   = matched_ret[:,[0,1]]
				mags_input_matched = matched_ret[:,[2,3]]
			else:   # x, y, mag, magerr, x, y, mag,magerr
				mags_ref_matched   = matched_ret[:,[2,3]]
				mags_input_matched = matched_ret[:,[6,7]]

		else:
			#prepare the standard stars
			print link_stdref_obs_method

			if std_world2img_refimage is None:
				if self.current_telescope == 'Iowa':
					outimg = os.path.join(self.template_dir, imgkey)
					refimg_abs = self.Iowa_prepare_template_image_wcs_single(imgkey, outimg=outimg)
				else:
					refimg_abs = self.images[imgkey]
				external_refimg = None
			else:
				refimg_abs = std_world2img_refimage
				external_refimg = std_world2img_refimage

			if std_catalog == 'apass':
				stds = self.__get_standard_reference_star_APASS(flt,wcs=True, single_use = True, refimg = refimg_abs)
			elif std_catalog == 'panstarrs':
				stds = self.__get_standard_reference_star_PanSTARRS(flt, wcs = True, photmethod = panstarrs_photmethod, reference_image = refimg_abs, single_use=True, savedatakey=imgkey, external_refimg = external_refimg, save_std_flt_file = True, save_std_flt_reg = True, r2R_flts = 'gr', i2I_flts = 'ri')
			elif std_catalog == '2mass':
				stds = self.__get_standard_reference_star_2MASS(flt,wcs=True)
			elif std_catalog.split(',')[0] == 'selfphot':
				std_catalog_segs = std_catalog.split(',')
				if len(std_catalog_segs)<4:
					raise ValueError('invalid input for std_catalog')
				elif len(std_catalog_segs) == 4:
					stdcatcode, telcode, imgindex, photmethod = std_catalog_segs
					imgindex2 = None
				elif len(std_catalog_segs) == 5:
					stdcatcode, telcode, imgindex, imgindex2, photmethod = std_catalog_segs
					imgindex2 = imgindex2+'.fits'
				else:
					raise ValueError('invalid input for std_catalog')

				refimgkey = imgindex+'.fits'
				refxymagsfile, refradecmagfile, bridge_image = self.__build_bridge_data_for_stdcal_from_selfphot(refimgkey, flt, telcode, photmethod, secondary_imgkey=imgindex2, saverefimg=None, xymagfile=None, radecmagfile=None, verbose=1)

				if std_world2img_refimage is None:
					refimg_abs = bridge_image
					external_refimg = bridge_image

				stds = self.__get_standard_reference_star_selfphot(refradecmagfile, refimg_abs, single_use=1, savedatakey=imgkey)
			else:
				raise ValueError('std catalog %s not supported'%std_catalog)

			if auto_shift_xy:
				inputimage = self.images[imgkey]
				target_xy_wcs = self.__get_xy_on_image_from_wcsinfo(inputimage)
				stds_shift_x =  self.photometry_info[imgkey]['x']- target_xy_wcs[0]
				stds_shift_y =  self.photometry_info[imgkey]['y']- target_xy_wcs[1]
				if verbose:
					print "(x, y) transformation from (RA, Dec) = (%s, %s)"%(target_xy_wcs[0], target_xy_wcs[1])
					print "(x, y) target coordinate on the image = (%s, %s)"%(self.photometry_info[imgkey]['x'], self.photometry_info[imgkey]['y'])
					print "%s stds shift x: %s; shift y:%s"%(imgkey, stds_shift_x, stds_shift_y)

			matched_ret = self.__stdcal_link_std_obs_single_image(imgkey,flt,method=link_stdref_obs_method, single_use = True, stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y, external_refimg=external_refimg)

			if self.reference_mag_min is not None:
				mag_faint_cut = self.reference_mag_min
				matched_ret = matched_ret[matched_ret[:,2]<mag_faint_cut,:]
			if self.reference_mag_max is not None:
				mag_bright_cut  = self.reference_mag_max
				matched_ret = matched_ret[matched_ret[:,2]>mag_bright_cut,:]
			if self.reference_magerr_cut is not None:
				magerr_cut_ref = self.reference_magerr_cut
				matched_ret = matched_ret[matched_ret[:,3]<magerr_cut_ref,:]
			if self.input_magerr_cut is not None:
				magerr_cut_input = self.input_magerr_cut
				matched_ret = matched_ret[matched_ret[:,7]<magerr_cut_input,:]

			ref_matched   = matched_ret[:,[0,1,2,3]]
			input_matched = matched_ret[:,[4,5,6,7]]
			mags_ref_matched = ref_matched[:,[2,3]]
			mags_input_matched = input_matched[:,[2,3]]

		if instmag is None:
			instmag = self.photometry_info[imgkey]['instmag']
		if instmagerr is None:
			instmagerr = self.photometry_info[imgkey]['instmagerr']

		offset_method = self.cal_offset_method
		cal_mag,cal_magerr, temp_future_use = self.__relative_calibration_single(mags_ref_matched,mags_input_matched, instmag, offset_method=offset_method)

		if update_stdmag:
			self.photometry_info[imgkey]['calmag'] = cal_mag
			self.photometry_info[imgkey]['calmagerr'] = np.sqrt( cal_magerr**2 + instmagerr**2 )

		return cal_mag, cal_magerr




	def __get_image_size(self, img):
		'''
		get image size from header 'NAXIS1' and 'NAXIS2'
		'''
		hdu = fits.open(img)
		header = hdu['PRIMARY'].header

		NX = header['NAXIS1']
		NY = header['NAXIS2']

		return NX,NY

	def __get_tpl_image_size(self,flt, tplimg=None):
		'''
		get the image size
		'''
		if tplimg is not None:
			NX,NY = self.__get_fits_image_size(tplimg)
		else:
			cal_template = self.templates_after_astrometry[flt]
			image  = fits.open(cal_template)
			header = image['PRIMARY'].header
			data = image['PRIMARY'].data

			#LCOGT use new CCD with different NAXIS1 and NAXIS2
			NX = None
			NY = None

			self.__find_template_imagekey()
        		ori_template = os.path.join(self.raw_image_dir,self.templates[flt])
			hdu_ori = fits.open(ori_template)
			header_ori = hdu_ori['PRIMARY'].header

			try:
				NX = header_ori['NAXIS1']
				NY = header_ori['NAXIS2']
			except:
				print "image size not found in the original tempalte image"

			print NX,NY


			if NX is None or NY is None:
				try:
					NX = header['NAXIS1']
					NY = header['NAXIS2']
				except:
					print "image size not found in the astrometry solution"

			if NX is None or NY is None:
				try:
					NX = int(self.base.telescopes_info[self.current_telescope]['nx'])
	        	                NY = int(self.base.telescopes_info[self.current_telescope]['ny'])
				except:
					raise ValueError("NAXIS1 and NAXIS2 not available in image header and please specify it though telescope infomation file")

			print NX,NY

		return NX,NY


	def __get_fits_image_size(self, image):

		hdu = fits.open(image)
		hdr = hdu['PRIMARY'].header
		NX = hdr['NAXIS1']
		NY = hdr['NAXIS2']
		return NX,NY



	def get_standard_reference_star_PanSTARRS(self, reference_image,  flt, wcs = True, photmethod = 'PSF'):
		'''
		INPUTS:
			reference_image:
			flt:
			wcs:
			photmethod:
		'''
		std_table = self.__get_standard_reference_star_PanSTARRS(flt, wcs = wcs, photmethod = photmethod, reference_image = reference_image, save_std_flt_file = False, save_std_flt_reg = False)

		if flt == 'gp':
			self.panstarrs_stds_gp = std_table
		elif flt == 'rp':
			self.panstarrs_stds_rp = std_table
		elif flt == 'ip':
			self.panstarrs_stds_ip = std_table
		elif flt == 'zp':
			self.panstarrs_stds_zp = std_table
		elif flt == 'yp':
			self.panstarrs_stds_yp = std_table
		else:
			raise ValueError("flt %s not supported"%flt)


	def __get_standard_reference_star_PanSTARRS(self, flt, wcs = True, photmethod = 'PSF', reference_image = None, single_use = False, imgkey=None, savedatakey=None,external_refimg=None, save_std_flt_file = True, save_std_flt_reg = True, r2R_flts = 'gr', i2I_flts = 'ri', wcsinfo_hdu='PRIMARY'):
		'''
		prepare PanSTARRS std stars and the expected product is a file containing img_x, img_y, mag_flt, magerr_flt where img_x and img_y are image
		coordinate related to the provided reference image and mag_flt is magnitudes in filter of interest. The reference image is provided in
		different ways under different circumstances, which are listed below.

		Reference image:
			(1) single_use --> False: the reference image is related to self.templates[flt]
			(2) single_use --> True:  'reference_image' is provided from internal data of current telescope+sn, 						  for example '001.fits'
			(3) 'external_refimg' is not None then will be used at higher priorities than (1), (2)

		INPUTS:
			flt:
			wcs: whether the input reference image has WCS info in header
			photmethod: 'PSF' or 'AP'; 'PSF' not supported for self.panstarrs_catalog_method=2
			r2R_flts: gr or ri
			i2I_flts: ri or iz

		Transformation:
		       R = r - 0.1837*(g - r) - 0.0971;  sigma = 0.0106
               	       R = r - 0.2936*(r - i) - 0.1439;  sigma = 0.0072
                       I = r - 1.2444*(r - i) - 0.3820;  sigma = 0.0078
                       I = i - 0.3780*(i - z)  -0.3974;  sigma = 0.0063
		'''
		#check the input filter
		flts_acceptable = ['gp','rp','ip','zp','yp', 'R', 'I']
		if flt not in flts_acceptable:
			raise IOError,'For now only %s acceptable !' % flts_acceptable

		#get the related panstarrs catalog data
		std_ref_stars_file_abs = os.path.join(self.std_ref_dir,self.std_ref_stars_file['panstarrs'])
		refstars = Table.read(std_ref_stars_file_abs, format='ascii.csv')
		refstars = refstars.filled(fill_value=99.99)

		#Below are for catalog obtained from MAST
		if self.panstarrs_catalog_method == 1:
			if photmethod == 'PSF':
				gmags = refstars['gMeanPSFMag']
				gmagerrs = refstars['gMeanPSFMagErr']
				rmags = refstars['rMeanPSFMag']
				rmagerrs = refstars['rMeanPSFMagErr']
				imags = refstars['iMeanPSFMag']
				imagerrs = refstars['iMeanPSFMagErr']
				zmags = refstars['zMeanPSFMag']
				zmagerrs = refstars['zMeanPSFMagErr']
			elif photmethod == 'AP':
				gmags = refstars['gMeanApMag']
				gmagerrs = refstars['gMeanApMagErr']
				rmags = refstars['rMeanApMag']
				rmagerrs = refstars['rMeanApMagErr']
				imags = refstars['iMeanApMag']
				imagerrs = refstars['iMeanApMagErr']
				zmags = refstars['zMeanApMag']
				zmagerrs = refstars['zMeanApMagErr']
			else:
				raise ValueError("invalid input for photmethod...")
		elif self.panstarrs_catalog_method == 2:
			gmags = refstars['gmag']
			gmagerrs = refstars['e_gmag']
			rmags = refstars['rmag']
			rmagerrs = refstars['e_rmag']
			imags = refstars['imag']
			imagerrs = refstars['e_imag']
			zmags = refstars['zmag']
			zmagerrs = refstars['e_zmag']
		else:
			raise ValueError('self.panstarrs_catalog_method valid option: 1,2 ')

		#prepare R/I mag from gri transformation
		if r2R_flts == 'gr':
			Rmags = rmags - 0.1837*(gmags - rmags) - 0.0971
			Rmagerrs = np.sqrt( rmagerrs**2 + (0.1837*gmagerrs)**2 + (0.1837*rmagerrs)**2 + 0.0106**2 )
		elif r2R_flts == 'ri':
			Rmags = rmags - 0.2936*(rmags - imags) - 0.1439
			Rmagerrs = np.sqrt( rmagerrs**2 + (0.2936*gmagerrs)**2 + (0.2936*rmagerrs)**2 + 0.0072**2 )

		#Rmag_col = Column(data=Rmags, name='Rmag')
		#Rmagerr_col = Column(data=Rmagerrs, name='Rmagerr')
		Rmag_col = Column(data=Rmags, name='Rmag')
		if self.panstarrs_catalog_method == 1:
			Rmagerr_col = Column(data=Rmagerrs, name='Rmagerr')
		else:
			Rmagerr_col = Column(data=Rmagerrs, name='e_Rmag')


		if i2I_flts == 'ri':
			Imags = rmags - 1.2444*(rmags - imags) - 0.3820		#I = r - 1.2444*(r - i) - 0.3820;  sigma = 0.0078
			Imagerrs = np.sqrt( rmagerrs**2 + (1.2444*rmagerrs)**2 + (1.2444*imagerrs)**2 + 0.0078**2 )
		elif i2I_flts == 'iz':
			Imags = imags - 0.3780*(imags - zmags) - 0.3974		#I = i - 0.3780*(i - z)  -0.3974;  sigma = 0.0063
			Imagerrs = np.sqrt( imagerrs**2 + (0.3780*imagerrs)**2 + (0.3780*zmagerrs)**2 + 0.0063**2 )

		#Imag_col = Column(data=Imags, name='Imag')
		#Imagerr_col = Column(data=Imagerrs, name='Imagerr')
		Imag_col = Column(data=Imags, name='Imag')

		if self.panstarrs_catalog_method == 1:
			Imagerr_col = Column(data=Imagerrs, name='Imagerr')
		else:
			Imagerr_col = Column(data=Imagerrs, name='e_Imag')


		refstars.add_columns([Rmag_col, Rmagerr_col, Imag_col, Imagerr_col])


		if single_use:
			if savedatakey is not None:
				savedata_key_s = savedatakey.split('.')[0]
			else:
				savedata_key = reference_image.split('/')[-1]
				savedata_key_s = savedatakey.split('.')[0]

		self.__get_stdref_filenames() #template images
		std_ref_big_file = self.flt_std_ref_stars_file_big[flt]
		std_ref_big_absfile = os.path.join(self.std_ref_dir, std_ref_big_file)
		self.__delete_file_if_exist(std_ref_big_absfile)

		#check the get the wcs info of the reference image
		if external_refimg is not None:
			if not os.path.isfile(external_refimg):
				raise IOError('external reference image %s does not exist'%external_refimg)
			image = fits.open(external_refimg)
			header = image[wcsinfo_hdu].header
		elif single_use:
			if reference_image is None:
				raise ValueError('the reference image refimg is required...')
			elif not os.path.isfile(reference_image):
				raise IOError('reference image %s does not exist'%reference_image)
			else:
				image = fits.open(reference_image)
				header = image[wcsinfo_hdu].header
		else:

			cal_template = self.templates_after_astrometry[flt]
			image  = fits.open(cal_template)
			header = image[wcsinfo_hdu].header

		if self.panstarrs_catalog_method == 1:
			if photmethod == 'PSF':
				colnames_dict = {'ra':'raMean', 'dec':'decMean','gp':'gMeanPSFMag', 'gperr':'gMeanPSFMagErr', 'rp':'rMeanPSFMag', 'rperr':'rMeanPSFMagErr', 'ip':'iMeanPSFMag','iperr':'iMeanPSFMagErr','zp':'zMeanPSFMag', 'zperr':'zMeanPSFMagErr', 'yp':'yMeanPSFMag', 'yperr':'yMeanPSFMagErr', 'R':'Rmag', 'Rerr':'Rmagerr', 'I':'Imag', 'Ierr':'Imagerr'}
			elif photmethod == 'AP':
				colnames_dict = {'ra':'raMean','dec':'decMean', 'gp':'gMeanApMag', 'gperr':'gMeanApMagErr', 'rp':'rMeanApMag', 'rperr':'rMeanApMagErr', 'ip':'iMeanApMag', 'iperr':'iMeanApMagErr', 'zp':'zMeanApMag', 'zperr':'zMeanApMagErr', 'yp':'yMeanApMag', 'yperr':'yMeanApMagErr', 'R':'Rmag', 'Rerr':'Rmagerr', 'I':'Imag', 'Ierr':'Imagerr'}
			else:
				raise ValueError('photmethod has to be PSF or AP')
		elif self.panstarrs_catalog_method == 2:
			colnames_dict = {'ra':'RAJ2000', 'dec':'DEJ2000', 'gp':'gmag', 'rp':'rmag', 'ip':'imag','zp':'zmag','yp':'ymag','gperr':'e_gmag', 'rperr':'e_rmag', 'iperr':'e_imag','zperr':'e_zmag','yperr':'e_ymag','R':'Rmag',		  'Rerr':'e_Rmag', 'I':'Imag', 'Ierr':'e_Imag'}
		else:
			raise ValueError('self.panstarrs_catalog_method valid option: 1, 2')

		if wcs:
			#stds filtering to get stars with good magnitudes
			#the current values can be modified to have real filtering effect
			mag_up = 99
			magerr_up = 99
			mag_low = -100
			magerr_low = -100
			refstars = refstars[refstars[colnames_dict[flt]]>mag_low]
			refstars = refstars[refstars[colnames_dict[flt+'err']]>magerr_low]

			refstars = refstars[refstars[colnames_dict[flt]]<mag_up]
			refstars = refstars[refstars[colnames_dict[flt+'err']]<magerr_up]
			print refstars

			if self.panstarrs_catalog_method ==1:
				if self.panstarrs_min_nx is not None:
					colname = 'n'+flt[0]
					refstars = refstars[refstars[colname]>self.panstarrs_min_nx]
				if self.panstarrs_min_nstackDetections is not None:
					refstars = refstars[refstars['nStackDetections'] > self.panstarrs_min_nstackDetections]

			w = WCS(header)
			ra  = refstars[colnames_dict['ra']].data
			dec = refstars[colnames_dict['dec']].data
			world = np.concatenate((ra.reshape(len(ra),1),dec.reshape(len(dec),1)),axis=1)
			pix = w.wcs_world2pix(world,1) # Pixel coordinates of (RA, DEC)

			#fits image treat first pixel as (1,1)
			#to make the star pixel position consistent with that of the result from sfind
			pix = pix-1
			xy_table = Table(pix,names=['img_x','img_y'])
			refstars_new = hstack((xy_table,refstars))

			refstars_new.write(std_ref_big_absfile,format='ascii.csv')
			print refstars_new

			cols_wanted = [colnames_dict[flt],colnames_dict[flt+'err']]
			mags_table = refstars[cols_wanted]
			std_table = hstack((xy_table, mags_table))
			print std_table

			#select stds in science image field
			if external_refimg is not None:
				NX = 9999
				NY = 9999
				lowcut = -9999
			if reference_image is None:
				NX,NY = self.__get_tpl_image_size(flt)
				lowcut = 0
			else:
				NX,NY = self.__get_fits_image_size(reference_image)
				lowcut = 0
			print NX,NY

			mask1 = np.logical_and(xy_table['img_x']>lowcut, xy_table['img_x']<NX)
			mask2 = np.logical_and(xy_table['img_y']>lowcut, xy_table['img_y']<NY)
			mask = mask1*mask2
			std_table = std_table[mask]
			print std_table

		else:
			raise IOError('under construction...')


                if self.panstarrs_mag_faint_cut is not None:
                        faint_mag_cut = self.panstarrs_mag_faint_cut
			std_table = std_table[std_table[colnames_dict[flt]]<faint_mag_cut]

                if self.panstarrs_mag_saturate_cut is not None:
                        saturate_mag_cut = self.panstarrs_mag_saturate_cut
			std_table = std_table[std_table[colnames_dict[flt]]>saturate_mag_cut]


		if save_std_flt_file:
			if not single_use:
				std_ref_flt = self.flt_std_ref_stars_file[flt]
				std_ref_filename_abs = os.path.join(self.std_ref_dir,std_ref_flt)
			else:
				std_ref_filename_abs = os.path.join(self.std_ref_dir, 'std_ref_'+savedata_key_s+'.txt' )

			self.__delete_file_if_exist(std_ref_filename_abs)
			std_table.write(std_ref_filename_abs, format='ascii.fast_no_header')

		if save_std_flt_reg:
		#save the sources as a ds9  reg file
			if not single_use:
				source_regfile =  os.path.join(self.std_ref_dir, flt + '.reg')
			else:
				source_regfile = os.path.join(self.std_ref_dir, savedata_key_s + '.reg')
			create_ds9_region_file(pix, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		return std_table




	def __get_standard_reference_star_2MASS(self,flt,wcs=True, eliminate =True, reference_image=None):
		'''
		self.standards have multiple colnames; This function is used to extract image position (x,y) and single band magnitude

		INPUTS:
			eliminate: eliminate stars which are outside the field of view of the image

		OUTPUTS:
			output file contains single band standard stars, with columns of 'x','y','mag','magerr'
		'''
		self.__get_stdref_filenames()

		#check the unput filter
		flts_acceptable = ['J','H','K']
		if flt not in flts_acceptable:
			raise IOError,'For now only %s acceptable !' % flts_acceptable


		std_ref_stars_file_abs = os.path.join(self.std_ref_dir,self.std_ref_stars_file['2mass'])
		refstars = Table.read(std_ref_stars_file_abs, format='ascii.csv')

		std_ref_big_file = self.flt_std_ref_stars_file_big[flt]
		std_ref_big_absfile = os.path.join(self.std_ref_dir, std_ref_big_file)
		self.__delete_file_if_exist(std_ref_big_absfile)

		if reference_image is None:
			cal_template = self.templates_after_astrometry[flt]
		else:
			if os.path.exists(reference_image):
				cal_template = reference_image
			else:
				raise ValueError("reference image provided %s not exists"%reference_image)

		image  = fits.open(cal_template)
		header = image['PRIMARY'].header

		if wcs:
			NX,NY = self.__get_tpl_image_size(flt, tplimg=reference_image)

			w = WCS(header)
			ra  = refstars['ra'].data.data	#refstars['ra'] is masked table column
			dec = refstars['dec'].data.data
			world = np.concatenate((ra.reshape(len(ra),1),dec.reshape(len(dec),1)),axis=1)

			pix = w.wcs_world2pix(world,1) # Pixel coordinates of (RA, DEC)

			#fits image treat first pixel as (1,1)
			#to make the star pixel position consistent with that of the result from sfind
			pix = pix-1
			xy_table = Table(pix,names=['img_x','img_y'])
			refstars_new = hstack((xy_table,refstars))

			refstars_new.write(std_ref_big_absfile,format='ascii.csv')

			colnames_dict = {'J':'j_m', 'Jerr':'j_msig', 'H':'h_m', 'Herr':'h_msig', 'K':'k_m', 'Kerr':'k_msig'}
			cols_wanted = [colnames_dict[flt],colnames_dict[flt+'err']]
			mags_table = refstars[cols_wanted]

			#empty values exist in the table, fill the empty with 99.99
			mags_table = mags_table.filled(fill_value=99.99)

			std_table = hstack((xy_table, mags_table))

			if eliminate:
				mask1 = np.logical_and(xy_table['img_x']>0, xy_table['img_x']<NX)
				mask2 = np.logical_and(xy_table['img_y']>0, xy_table['img_y']<NY)
				mask = mask1*mask2
				std_table = std_table[mask]

		else:
			raise IOError('under construction...')

		std_ref_flt = self.flt_std_ref_stars_file[flt]
		std_ref_flt_abs = os.path.join(self.std_ref_dir,std_ref_flt)
		self.__delete_file_if_exist(std_ref_flt_abs)
		std_table.write(std_ref_flt_abs, format='ascii.fast_no_header')

		#save the sources as a ds9  reg file
		source_regfile =  os.path.join(self.std_ref_dir, flt + '.reg')
		create_ds9_region_file(pix, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		return std_table


	def __get_standard_reference_star_selfphot(self, radecmagfile, refimg, single_use=1, wcs=1, external_refimg =None,  verbose=1, savedatakey=None):
		'''
		INPUTS:
			radecmagfile:
			refimg:
			single_use:
			wcs: whether the refimg has wcs info in header;
			verbose:
		'''

		if single_use:
			if savedatakey is not None:
				savedata_key_s = savedatakey.split('.')[0]
			else:
				savedata_key = refimg.split('/')[-1]
				savedata_key_s = savedata_key.split('.')[0]
			if verbose:
				print "image %s will be used to relate the std catalog and observed sources"%refimg


		refstars = np.loadtxt(radecmagfile, skiprows=0)
		if verbose:
			print "the primitive std stars list are below:"
			print refstars

		#you can exclude standard stars with too big uncertainties
		if self.selfphot_magerr_up_cut is not None:
			refstars = refstars[refstars[:,3]<self.selfphot_magerr_up_cut]
			if verbose:
				print "stds with uncertainties >%s removed"%self.selfphot_magerr_up_cut

		if single_use:
			if refimg is None:
				raise ValueError('the reference image refimg is required...')
			elif not os.path.exists(refimg):
				raise IOError('reference image %s not exist'%refimg)
			else:
				image = fits.open(refimg)
				header = image['PRIMARY'].header
		else: #dead code block
			std_ref_big_file = self.flt_std_ref_stars_file_big[flt]
			std_ref_big_absfile = os.path.join(self.std_ref_dir, std_ref_big_file)
			self.__delete_file_if_exist(std_ref_big_absfile)
			cal_template = self.templates_after_astrometry[flt]
			image  = fits.open(cal_template)
			header = image['PRIMARY'].header

		if wcs:
			NX,NY = self.__get_image_size(refimg)
			w = WCS(header)
			ra  = refstars[:,0]
			dec = refstars[:,1]
			world = np.concatenate((ra.reshape(len(ra),1),dec.reshape(len(dec),1)),axis=1)
			pix = w.wcs_world2pix(world,1) # Pixel coordinates of (RA, DEC)

			#fits image treat first pixel as (1,1)
			#to make the star pixel position consistent with that of the result from sfind
			pix = pix-1
			refstars_new = np.concatenate((pix,refstars),axis=1)

			mask  = np.logical_and(refstars_new[:,4]<90,refstars_new[:,4]<90)
			std = refstars_new[mask,:]
			std = std[:,[0,1,4,5]]

			if verbose:
				print "updated std stars list are below:"
				print std

		if self.selfphot_mag_faint_cut is not None:
			faint_mag_cut = self.selfphot_mag_faint_cut
			std  = std[std[:,2]<faint_mag_cut]
			if verbose:
				print "brighter than %s sources survived"%faint_mag_cut

		if self.selfphot_mag_saturate_cut is not None:
			saturate_mag_cut = self.selfphot_mag_saturate_cut
			std = std[std[:,2]>saturate_mag_cut]
			if verbose:
				print "fainter than %s sources survived"%saturate_mag_cut

		if wcs:
			#only the standard stars within the input image region are needed
			mask1 = np.logical_and(std[:,0]>0,std[:,0]<NX)
			mask2 = np.logical_and(std[:,1]>0,std[:,1]<NY)
			mask = mask1*mask2
			std = std[mask,:]


		if not single_use:
			std_ref_flt = self.flt_std_ref_stars_file[flt]
			std_ref_filename_abs = os.path.join(self.std_ref_dir,std_ref_flt)
		else:
			std_ref_filename_abs = os.path.join(self.std_ref_dir, 'std_ref_'+savedata_key_s+'.txt' )

		self.__delete_file_if_exist(std_ref_filename_abs)
		fid = open(std_ref_filename_abs,'w')
		np.savetxt(fid,std,fmt="%8.4f %8.4f %8.3f %8.3f")
		fid.close()

		#save the sources as a ds9  reg file
		if not single_use:
			source_regfile = os.path.join(self.std_ref_dir, flt + '.reg')
		else:
			source_regfile = os.path.join(self.std_ref_dir, savedata_key_s + '.reg')

		create_ds9_region_file(std, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)
		return std



	def __build_bridge_data_for_stdcal_from_selfphot(self, imgkey, flt, telcode, photmethod, secondary_imgkey=None, saverefimg=None, xymagfile=None, radecmagfile=None, verbose=1):
		'''
		prepare standard stars for calibration from photometry result for images of other telescope

		reply on the consistent file structure

		For R, I band, we need conversion from r,i band for which the secondary image from another filter is required

		INPUTS:
			imgkey: for example, '001.fits'
			telcode: for example, 'LCOGT', 'LT'
			photmethod: 'apphot', 'psfphot'
			saverefimg: the filename of the image to be saved
			xymagfile: the std catalog with format x,y,mag,magerr
			radecmagfile: the std catalog with format ra, dec, mag, magerr
		'''
		srcimgdir = self.raw_image_dir.replace(self.current_telescope, telcode)
		srcretdir = self.result_dir.replace(self.current_telescope, telcode)

		if photmethod == 'apphot':
			suffix1 = '.txt'	# photometry results table file
			photretdir = self.aperture_photometry_dir.replace(self.current_telescope, telcode)
		elif photmethod == 'psfphot':
			suffix1 = '_PSF.txt'
			photretdir = self.psf_photometry_dir.replace(self.current_telescope, telcode)
		else:
			raise ValueError('invalid input for photmethod')

		photrettable = os.path.join(srcretdir, '%s_photometry_info%s'%(self.current_sn, suffix1))
		photfile = os.path.join(photretdir, imgkey.replace('fits', photmethod))
		photret = Table.read(photrettable, format='ascii.fixed_width')

		imgphot = photret[photret['name']==imgkey]
		zptmag = imgphot['magzpt'][0]
		zptmagerr = imgphot['magzpterr'][0]
		flt1 = imgphot['flt'][0]

		if zptmag == 99.99 or zptmagerr == 99.99:
			instmag = imgphot['instmag'][0]
			instmagerr = imgphot['instmagerr'][0]
			calmag = imgphot['calmag'][0]
			calmagerr = imgphot['calmagerr'][0]
			zptmag = calmag - instmag
			zptmagerr = np.sqrt(calmagerr**2-instmagerr**2)

		xymags1 = np.loadtxt(photfile)
		xymags1[:,2] = xymags1[:,2] + zptmag
		xymags1[:,3] = np.sqrt(xymags1[:,3]**2 + zptmagerr**2)

		if secondary_imgkey is not None:
			photfile2 = os.path.join(photretdir, secondary_imgkey.replace('fits', photmethod))
			imgphot = photret[photret['name']==secondary_imgkey]
			zptmag = imgphot['magzpt'][0]
			zptmagerr = imgphot['magzpterr'][0]
			flt2 = imgphot['flt'][0]

			if zptmag == 99.99 or zptmagerr == 99.99:
				instmag = imgphot['instmag'][0]
				instmagerr = imgphot['instmagerr'][0]
				calmag = imgphot['calmag'][0]
				calmagerr = imgphot['calmagerr'][0]

				zptmag = calmag - instmag
				zptmagerr = np.sqrt(calmagerr**2-instmagerr**2)

			xymags2 = np.loadtxt(photfile2)
			xymags2[:,2] = xymags2[:,2] + zptmag
			xymags2[:,3] = np.sqrt(xymags2[:,3]**2 + zptmagerr**2)


			if self.selfphot_match_method_for_RI_stds == 'grmatch':
				raise ValueError("grmatch method for matching %s and %s under development"%(imgkey,secondary_imgkey))
				#grmatch_option = self.grmatch_option_obsimg2refimg
				#matched_table = self.__fitsh_grmatch(photfile, photfile2, match_output_filename, trans_fitting_filename, option=grmatch_option)
			elif self.selfphot_match_method_for_RI_stds == 'surrounding_search':
				tg_xys  = xymags2[:,0:2]
	        	        ref_xys = xymags1[:,0:2]
	        	        criteria = self.tpl_obs_match_criteria
        		        target_stars,indexmat,ref_stars,indexref = self.__crude_match(tg_xys,ref_xys,criteria)
				if len(indexmat) != 0:
			                match_ret = np.hstack((xymags1[indexref],xymags2[indexmat]))
				else:
					print "No match result obtained !!! You can:"
					print "Try to modify self.tpl_obs_match_criteria and see if improve"
					print "Try to modify self.selfphot_match_method_for_RI_stds to grmatch and run agian"
			else:
				raise ValueError('invalid input for self.selfphot_match_method_for_RI_stds')

			if flt == 'R':
				if flt1 == 'rp' and flt2 == 'ip':
					xymags = match_ret[:,[0,1,2,3]]
					rmag = match_ret[:,2]
					rmagerr = match_ret[:,3]
					imag = match_ret[:,6]
					imagerr = match_ret[:,7]
					Rmag = rmag - 0.2936*(rmag - imag) - 0.1439
					sigma = 0.0072
					Rmagerr = np.sqrt(sigma**2 + np.sqrt((0.2936*imagerr)**2+(0.7064*rmagerr)**2))  #uncertainty come from two parts, ri measurements and transformation systematic error
					xymags[:,2] = Rmag
					xymags[:,3] = Rmagerr
				else:
					raise ValueError("when prepare for R band std, r band should be main filter for transformation data and i band be secondary")

			elif flt == 'I':
				if flt1 == 'ip' and flt2 == 'rp':

					if verbose:
						print "SDSS-i band stds:"
						print match_ret[:, [0,1,2,3]]
						print "SDSS-r band stds:"
						print match_ret[:,[4,5,6,7]]

					xymags = match_ret[:,[0,1,2,3]]
					imag = match_ret[:,2]
					imagerr = match_ret[:,3]
					rmag = match_ret[:,6]
					rmagerr = match_ret[:,7]
					Imag = rmag - 1.2444*(rmag - imag) - 0.3820
					sigma = 0.0078
					Imagerr = np.sqrt(sigma**2 + np.sqrt((1.2444*imagerr)**2+(0.2444*rmagerr)**2))  #uncertainty come from two parts, ri measurements and transformation systematic error
					xymags[:,2] = Imag
					xymags[:,3] = Imagerr
				else:
					raise ValueError("when prepare for I band std, i band should be main filter for transformation data and r band be secondary")

			else:
				raise ValueError("std transformation only works for R and I")

		else:
			if flt != flt1:
				raise ValueError("The obtained filter doesn't match the desired one")
			xymags = xymags1



		if not xymagfile:
			savefile = os.path.join(self.std_ref_dir, 'xymag_'+telcode+imgkey.replace('fits', 'txt'))
		else:
			savefile = os.path.join(self.std_ref_dir, xmagfile)
		self.__delete_file_if_exist(savefile)
		fid = open(savefile,'w')
		np.savetxt(fid, xymags, fmt="%8.4f %8.4f %8.3f %8.3f")
		fid.close()

		if not saverefimg:
			saveimg =  os.path.join(self.std_ref_dir, telcode+imgkey)
		else:
			saveimg = os.path.join(self.std_ref_dir, saverefimg)
		self.__delete_file_if_exist(saveimg)
		if verbose:
			print "Image %s from telescope %s will be saved as %s"%(imgkey, telcode, saveimg)
		shutil.copy(os.path.join(srcimgdir, imgkey), saveimg)


		radecmags = xymags.copy()
		radecmags[:,0:2] = self.__image2world(saveimg, xymags[:,0:2], hduindex=None)
		if not radecmagfile:
			savefile2 = os.path.join(self.std_ref_dir, 'radecmag_'+telcode+imgkey.replace('fits', 'txt'))
		else:
			savefile2 = os.path.join(self.std_ref_dir, radecmagfile)
		self.__delete_file_if_exist(savefile2)
		fid2 = open(savefile2,'w')
		np.savetxt(fid2, radecmags, fmt="%9.6f %9.6f %8.3f %8.3f")
		fid2.close()

		return savefile, savefile2, saveimg


	def __get_standard_reference_star_APASS(self,flt,wcs=True, single_use = False, refimg = None, center_ra = None, center_dec = None, center_X = None, center_Y = None, verbose=0):
		'''
		Currently the calibrated star are from APASS DR9 which have Johson BV band and SDSS gri band photometry
		And we extend the use to Cousin R I band by transforming the SDSS r i band to I R according to
		'https://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php#Lupton2005'

		    R = r - 0.1837*(g - r) - 0.0971;  sigma = 0.0106
		    R = r - 0.2936*(r - i) - 0.1439;  sigma = 0.0072

		    I = r - 1.2444*(r - i) - 0.3820;  sigma = 0.0078
		    I = i - 0.3780*(i - z)  -0.3974;  sigma = 0.0063        (no z band in APASS)

		This function will create convert APASS standards from world coordinate to image coordinate.
		you will need one reference image to determine the transformation and region of interest.
		If 'single_use' is not True, then the template image for filter 'flt' will be used as the reference image;
		otherwise, 'regimg' is required

		INPUTS:
			flt:
			wcs: if True, the wcs information exists in reference image header and will be used
			single_use:
			refimg:
			center_ra:
			center_dec:
			center_X:
			center_Y:
		OUTPUTS:

		'''
		if single_use:
			if refimg is None:
				raise ValueError('the reference image refimg is required...')
			if not os.path.exists(refimg):
				raise IOError('reference image %s not exist'%refimg)
			refimg_key = refimg.split('/')[-1]
			refimg_key_s = refimg_key.split('.')[0]
			if verbose:
				print "image %s will be used to relate the APASS catalog and observed sources"%refimg
		else:
			refimg  = self.templates_after_astrometry[flt]


		#check the input filter
		flts_acceptable = ['B','V','I','R','gp','rp','ip','g','r','i']
		if flt not in flts_acceptable:
			raise IOError,'For now only %s acceptable !' % flts_acceptable


		self.__get_stdref_filenames()
		std_ref_stars_file_abs = os.path.join(self.std_ref_dir,self.std_ref_stars_file['apass'])
		if self.apass_catalog_method == 1:
			refstars = Table.read(std_ref_stars_file_abs, format='ascii')
		elif self.apass_catalog_method == 2:
			refstars = Table.read(std_ref_stars_file_abs, format='csv')
		else:
			raise ValueError('check')

		#you can exclude standard stars with insufficient observation
		if self.apass_mobs_min is not None:
			refstars = refstars[refstars['mobs']>self.apass_mobs_min]
			if verbose:
				print "APASS stds with mobs >%s selected"%mobs_min
		if self.apass_nobs_min is not None:
			refstars = refstars[refstars['nobs']>self.apass_nobs_min]
			if verbose:
				print "APASS stds with nobs >%s selected"%nobs_min

		hdu  = fits.open(refimg)
		header = hdu['PRIMARY'].header
		NX,NY = self.__get_image_size(refimg)

		refstarcolnames = refstars.colnames
		if ('RAJ2000' in refstarcolnames) and ('DEJ2000' in refstarcolnames):
			RA_colname  = 'RAJ2000'
			Dec_colname = 'DEJ2000'
		elif ('RA' in refstarcolnames) and ('Dec' in refstarcolnames):
			RA_colname  = 'RA'
			Dec_colname = 'Dec'
		else:
			raise ValueError("catalog colnames are not supported yet...")
		ra  = refstars[RA_colname]
		dec = refstars[Dec_colname]
		if wcs:
			w = WCS(header)
			world = np.concatenate((ra.reshape(len(ra),1),dec.reshape(len(dec),1)),axis=1)
			pix = w.wcs_world2pix(world,1) # Pixel coordinates of (RA, DEC)
			pix = pix-1 #fits image treat first pixel as (1,1),to make the star pixel position consistent with that of the result from sfind
			xs = pix[:,0]
			ys = pix[:,1]
		else:
			pixscale_deg = self.pixscale/60./60.
			if center_ra is None or center_dec is None:
				center_ra  = self.sn_ra_world_deg
				center_dec = self.sn_dec_world_deg
				print "target coordinate used as image center coordinate"

                        if center_X is None or center_Y is None:
                                if not single_use:
                                        self.__find_template_imagekey()
                                        tpl_imgkey = self.templates[flt]
                                        center_X = self.photometry_info[tpl_imgkey]['x']
                                        center_Y = self.photometry_info[tpl_imgkey]['y']
                                else:
                                        center_X = self.photometry_info[refimg_key]['x']
                                        center_Y = self.photometry_info[refimg_key]['y']
				print "target position (x,y) on image used as reference point"

			xs = -(ra - center_ra)/pixscale_deg  + X_cen
			ys = (dec - center_dec)/pixscale_deg + Y_cen


		xcol = Column(data=xs, name='x') # in pixel
		ycol = Column(data=ys, name='y')
		refstars.add_columns([xcol,ycol], [0,0])

		if refstars.masked:
			refstars = refstars.filled(99.999)

		if flt in ['I','R']: #transformation needed
			refstars = refstars[(refstars['r_mag']<50)*(refstars['i_mag']<50)]
			refstars = refstars[(refstars['r_mag']!=0)*(refstars['i_mag']!=0)]
			rmag = refstars['r_mag'].data
			rmagerr = refstars['e_r_mag'].data
                        imag = refstars['i_mag'].data
                        imagerr = refstars['e_i_mag'].data

			if flt == 'I':
				mag = rmag - 1.2444*(rmag - imag) - 0.3820
				magerr = np.sqrt(0.0078**2 + np.sqrt((1.2444*imagerr)**2+(0.2444*rmagerr)**2))  #uncertainty come from two parts, ri measurements and transformation systematic error (sigma = 0.0078)
			else:
				mag = rmag - 0.2936*(rmag - imag) - 0.1439
				magerr = np.sqrt(0.0072**2 + np.sqrt((0.2936*imagerr)**2+(0.7064*rmagerr)**2))  #uncertainty come from two parts, ri measurements and transformation systematic error (sigma = 0.0072)
		else:
			if flt in ['B','V']:
				magcolname = flt+'mag'
			elif flt in ['g','r','i']:
				magcolname = flt+'_mag'
			else:
				magcolname = flt[0]+'_mag'
			magerrcolname = 'e_'+magcolname

			refstars = refstars[refstars[magcolname]<50]
			refstars = refstars[refstars[magcolname]!=0]
			if self.apass_remove_single_measurement:
				refstars = refstars[refstars[magerrcolname]!=0]
			mag = refstars[magcolname].data
			magerr = refstars[magerrcolname].data

		if not single_use:
			std_ref_big_absfile = os.path.join(self.std_ref_dir, self.flt_std_ref_stars_file_big[flt])
		else:
			std_ref_big_absfile = os.path.join(self.std_ref_dir, 'std_ref_whole_info_'+refimg_key_s+'.txt' )

		self.__delete_file_if_exist(std_ref_big_absfile)
		refstars.write(std_ref_big_absfile,format='csv')


		x = refstars['x'].data
		y = refstars['y'].data
		std = np.array([x,y,mag,magerr]).transpose()

		if self.apass_mag_faint_cut is not None:
			std  = std[std[:,2]<self.apass_mag_faint_cut]
			if verbose:
				print "brighter than %s sources survived"%faint_mag_cut

		if self.apass_mag_saturate_cut is not None:
			std = std[std[:,2]>self.apass_mag_saturate_cut]
			if verbose:
				print "fainter than %s sources survived"%saturate_mag_cut

		if self.apass_remove_single_measurement:
			std = std[std[:,3]>0]
			if verbose:
				print "measurements with magnitude uncertainty = 0 removed"

		mask = ((std[:,0]>0)*(std[:,0]<NX))*((std[:,1]>0)*(std[:,1]<NY))
		std = std[mask,:] #only the standard stars within the input image region are needed

		if not single_use:
			std_ref_flt = self.flt_std_ref_stars_file[flt]
			std_ref_filename_abs = os.path.join(self.std_ref_dir,std_ref_flt)
		else:
			std_ref_filename_abs = os.path.join(self.std_ref_dir, 'std_ref_'+refimg_key_s+'.txt' )

		self.__delete_file_if_exist(std_ref_filename_abs)
		fid = open(std_ref_filename_abs,'w')
		np.savetxt(fid,std,fmt="%8.4f %8.4f %8.3f %8.3f")
		fid.close()

		#save the sources as a ds9  reg file
		if not single_use:
			source_regfile = os.path.join(self.std_ref_dir, flt + '.reg')
		else:
			source_regfile = os.path.join(self.std_ref_dir, refimg_key_s + '.reg')

		create_ds9_region_file(std, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		return std


	def __get_stdref_filenames(self):
		'''
		prepare some filename for standard catalog
		'''
		self.__find_template_imagekey()
		for flt in self.templates.keys():
			tpl_imgkey = self.templates[flt]
			tpl_imgkey_s = tpl_imgkey.split('.')[0]
			self.flt_std_ref_stars_file[flt] = 'std_ref_' +flt+ '.txt'
			self.flt_std_ref_stars_file_big[flt] = 'std_ref_whole_info_' +flt+ '.txt'
			self.flt_std_ref_and_obs_measurement_stars_file[flt] = tpl_imgkey_s +'_std_'+flt+'.match'



	def get_relative_mags(self, tpl_obs_match_method = 'grmatch', updatetable=1):
		'''
		relative magnitude calibration
		'''
		self.__find_template_imagekey(updatetable=updatetable)
		flts_unique = self.templates.keys()
		offset_method = self.cal_offset_method
		for flt in flts_unique:
			print "Now working on %s band"%flt
			self.__get_relative_mags_flt(flt, tpl_obs_match_method = tpl_obs_match_method, offset_method=offset_method, updatetable=updatetable)


	def __get_relative_mags_flt(self,flt, tpl_obs_match_method = 'grmatch', offset_method='median', updatetable=1):

		if updatetable:
			self.__dict2table()
		self.__find_template_imagekey(updatetable=updatetable)

		tpl_imgkey = self.templates[flt]
		images_flt = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)

		for image_key in images_flt['name']:
			print image_key
			if self.photometry_info[image_key]['drop'] > 0:
				continue

			if self.photometry_info[image_key]['relmag'] != 99.99:
				if not self.renew_relative_mag:
					continue


			self.__get_relative_mag_single(image_key,tpl_imgkey, tpl_obs_match_method = tpl_obs_match_method, offset_method=offset_method)



	def standard_calibration_flt_manual_match_method(self, flt, which_dir='raw_image'):
		'''
		standard calibration for situation where automatic matching doesn't work
		'''
		if len(self.templates) == 0:
			self.__find_template_imagekey()

		if flt not in self.templates.keys():
			raise ValueError("template image in filter %s not specified..."%flt)

		tpl_imgkey = self.templates[flt]

		if self.photometry_info[tpl_imgkey]['calmag'] == 99.99 or self.renew_stdcal_mag:
			self.__get_tpl_calibrated_mag_manual_method(tpl_imgkey, flt, which_dir=which_dir)
		else:
			print "old calibration result will be used for the template image"

		self.__dict2table()
		table_flt = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)

		for image_key in table_flt['name']:
			if self.photometry_info[image_key]['template'] == 1:
				continue

			if self.photometry_info[image_key]['drop'] > 0:
				continue

			offset = self.photometry_info[tpl_imgkey]['calmag'] - self.photometry_info[tpl_imgkey]['relmag']
			self.photometry_info[image_key]['calmag'] = self.photometry_info[image_key]['relmag'] + offset

			instmagerr = self.photometry_info[image_key]['instmagerr']
			relmagerr = self.photometry_info[image_key]['relmagerr']
			stdmagerr = np.sqrt(self.photometry_info[tpl_imgkey]['calmagerr']**2 - self.photometry_info[tpl_imgkey]['instmagerr']**2)

			#Improvement needed!!!
			if self.photometry_info[tpl_imgkey]['drop'] >0 or stdmagerr==99.99:
				stdmagerr = 2*relmagerr

			calmagerr = np.sqrt(instmagerr**2 + relmagerr**2 + stdmagerr**2)
			self.photometry_info[image_key]['calmagerr'] = calmagerr


	def __get_tpl_calibrated_mag_manual_method(self, tpl_imgkey, flt, which_dir ='raw_image'):
		'''
		select the standard stars which match the reference stars on the template image and perform the calibration

		What required before this:
			self.get_template_reference_mags_flt
			self.show_standards_for_wcs_not_available
		'''
		input_image = self.__get_internal_image(tpl_imgkey, which_dir=which_dir)
		tpl_skey = tpl_imgkey.split('.')[0]

		std_xymags_file = os.path.join(self.std_ref_dir, "%s_%s_std_xymag.txt"%(flt, tpl_skey))
		if not os.path.exists(std_xymags_file):
			raise IOError("%s not exists; please get it ready first from function show_standards_for_wcs_not_available"%std_xymags_file)

		xysmags = np.loadtxt(std_xymags_file)

		regionfile = os.path.join(self.std_ref_dir, "%s_%s_std_xymag.reg"%(flt, tpl_skey))

		stdcal_idxymag_table = self.__input_xys_through_ds9_get_mags(input_image, xysmags, regionfile=regionfile)

		stdcalmag_filename  = os.path.join(self.template_dir, "%s_%s_std_idxymag.txt"%(flt, tpl_skey))
		stdcal_idxymag_table.write(stdcalmag_filename, format='ascii.fixed_width')


		refmag_filename = os.path.join(self.template_dir, "%s_%s_idxymag.txt"%(flt, tpl_skey))
		if not os.path.exists(refmag_filename):
			raise ValueError("please get the reference mags first; the function will be used get_template_reference_mags_flt ")

		refmag_table = Table.read(refmag_filename, format='ascii.fixed_width')

		stdmags	= []
		stdmagerrs = []
		tplmags = []
		tplmagerrs = []

		for stdstar in stdcal_idxymag_table:
			sid = stdstar['num']

			if sid not in refmag_table['num']:
				raise ValueError("reference star with id number %s not found in reference table"%sid)

			tplstar = refmag_table[refmag_table['num']==sid][0]

			stdmags.append(stdstar['mag'][0])
			stdmagerrs.append(stdstar['magerr'][0])

			tplmags.append(tplstar['mag'])
			tplmagerrs.append(tplstar['magerr'])

		print stdmags

		mags_ref_matched = np.transpose(np.array([stdmags, stdmagerrs]))
		mags_input_matched = np.transpose(np.array([tplmags, tplmagerrs]))

		print mags_ref_matched
		print mags_input_matched

		tpl_instmag = self.photometry_info[tpl_imgkey]['instmag']
		tpl_instmagerr = self.photometry_info[tpl_imgkey]['instmagerr']

		cal_mag,cal_magerr, temp_future_use = self.__relative_calibration_single(mags_ref_matched,mags_input_matched,tpl_instmag)


		self.photometry_info[tpl_imgkey]['calmag'] = cal_mag
		self.photometry_info[tpl_imgkey]['calmagerr'] =  np.sqrt( cal_magerr**2 + tpl_instmagerr**2 )



	def __get_relative_mags_flt_manual_match_method(self, flt, which_dir='raw_image', verbose=1, update_matched_inmags=True):

		self.__dict2table()
		self.__find_template_imagekey()

		tpl_imgkey = self.templates[flt]
		images_flt = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)

		for image_key in images_flt['name']:
			print image_key
			if self.photometry_info[image_key]['drop'] > 0:
				continue

			if self.photometry_info[image_key]['relmag'] != 99.99:
				if not self.renew_relative_mag:
					continue

			if image_key == tpl_imgkey:
				self.photometry_info[image_key]['relmag'] = self.photometry_info[image_key]['instmag']
				self.photometry_info[image_key]['relmagerr'] = 0.0
				continue

			relmag, relmagerr = self.__get_relative_mag_single_manual_match(image_key, flt, which_dir=which_dir, verbose=verbose, update_matched_inmags=update_matched_inmags)
			self.photometry_info[image_key]['relmag'] = relmag
			self.photometry_info[image_key]['relmagerr'] = relmagerr


	def __get_relative_mag_single_manual_match(self, image_key, flt, which_dir ='raw_image', verbose=1, update_matched_inmags = False, refstarnum=None,  refnumin=5, average_method='weighted', update_rettable=False, ds9name=None):
		'''
		relative calibration beween template image and input image. the matching is done by interactively selecting

		INPUTS:
			image_key:
			flt:
			which_dir:
			verbose:
			update_matched_inmags: update the matched mag list if it already exists
			refstarnum: the number of manually seletected reference star
			refnumin: the minimum number of reference stars required for the method of self.__relative_calibration_single; otherwise simple weighted average will be used
			average_method: 'weighted' or 'simple'
		'''

		if len(self.templates) == 0:
			self.__find_template_imagekey()

		if flt not in self.templates.keys():
			raise IOError("template image for %s not available"%flt)

		tpl_imgkey = self.templates[flt]
		tpl_imgkey_s = tpl_imgkey.split('.')[0]

		image_key_s = image_key.split('.')[0]

		refmag_filename = os.path.join(self.template_dir, "%s_%s_idxymag.txt"%(flt, tpl_imgkey_s))
		if not os.path.exists(refmag_filename):
			raise ValueError("please get the reference mags first from function get_template_reference_mags_flt")

		refmag_table = Table.read(refmag_filename, format='ascii.fixed_width')
		if verbose:
			print "mag table for reference stars on template image:"
			print refmag_table

		inmag_table = self.__relative_calibration_prepare_interactive_matching( image_key, which_dir=which_dir, clobber=update_matched_inmags, refstarnum=refstarnum, ds9name=ds9name)
		if verbose:
			print "mag table for reference stars on input image:"
			print inmag_table

		refmags	= []
		refmagerrs = []
		inmags = []
		inmagerrs = []
		for instar in inmag_table:
			sid = instar['num']

			if sid not in refmag_table['num']:
				raise ValueError("reference star with id number %s not found in reference table"%sid)

			refstar = refmag_table[refmag_table['num']==sid][0]

			inmags.append(instar['mag'])
			inmagerrs.append(instar['magerr'])

			refmags.append(refstar['mag'])
			refmagerrs.append(refstar['magerr'])

		refmags_match = np.transpose(np.array([refmags, refmagerrs]))
		inmags_match  = np.transpose(np.array([inmags, inmagerrs]))
		instmag = self.photometry_info[image_key]['instmag']

		if len(refmags)<refnumin:
			if verbose:
				print "less than %s reference stars will be used for the relative calibration. The %s average of the offsets will be applied and std of the offsets will be used as uncertainty"%(refnumin, average_method)

			if average_method == 'weighted':
				offsets = np.array(refmags) - np.array(inmags)
				offsetserr = np.sqrt(np.array(refmagerrs)**2 + np.array(inmagerrs)**2)
				print "The offsets and uncertainties on offsets are:"
				print offsets
				print offsetserr

				weights = 1./offsetserr**2
                        	offset_ave = np.sum(weights*offsets) / np.sum(weights)
                       		#offseterr_ave = 1./ np.sqrt(np.sum(weights))
				relmag = instmag + offset_ave
				relmagerr = np.std(offsets)

			elif average_method == 'simple':
				offset = np.mean(np.array(refmags) - np.array(inmags))
				relmag = instmag + offset
				relmagerr = np.std(offset)
			else:
				raise IOError("invalid input for average method")
		else:
			relmag, relmagerr, temp_future_use = self.__relative_calibration_single(refmags_match, inmags_match, instmag)

		if verbose:
			print "relmag=%s, relmagerr=%s"%(relmag, relmagerr)

		if update_rettable:
			self.photometry_info[image_key]['relmag'] = relmag
			self.photometry_info[image_key]['relmagerr'] = relmagerr

		return relmag, relmagerr


	def __get_relative_mag_single(self,image_key,tpl_imgkey,  tpl_obs_match_method = 'grmatch', offset_method='median'):
		'''
		relative calibration between template image and input image
		'''
		instmag = self.photometry_info[image_key]['instmag']
                instmagerr = self.photometry_info[image_key]['instmagerr']

                if image_key == tpl_imgkey:
                	relmag = instmag
                	relmagerr = 0
		else:
         	       	mags_ref_matched,mags_input_matched = self.__relative_calibration_prepare(image_key,tpl_imgkey, tpl_obs_match_method = tpl_obs_match_method)

			if self.photometry_info[image_key]['drop']>0:
				print "image %s has been dropped... please check"%image_key
				return

			if self.photometry_method   == 'apphot':
				match_prereject_filename  = os.path.join(self.aperture_photometry_dir,image_key.split('.')[0]+ '_tpl_tstrmd.match')
				match_caldata_filename  = os.path.join(self.aperture_photometry_dir,image_key.split('.')[0]+ '_tpl_caldata.match')
			elif self.photometry_method == 'psfphot':
        	       	        match_prereject_filename  = os.path.join(self.psf_photometry_dir,image_key.split('.')[0]+ '_tpl_tstrmd.match')
        	       	        match_caldata_filename  = os.path.join(self.psf_photometry_dir,image_key.split('.')[0]+ '_tpl_caldata.match')
			else:
				raise IOError("Invalid input for photometry_method")

			relmag,relmagerr, caldata_index  = self.__relative_calibration_single(mags_ref_matched,mags_input_matched,instmag, offset_method=offset_method)
			match_alldata = np.loadtxt(match_prereject_filename)
			match_caldata = match_alldata[caldata_index,:]
			np.savetxt(match_caldata_filename, match_caldata, fmt="%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f")

            	print relmag,relmagerr
                self.photometry_info[image_key]['relmag'] = relmag
                self.photometry_info[image_key]['relmagerr'] = relmagerr





	def __relative_calibration_prepare_interactive_matching(self,image_key, which_dir='raw_image', clobber=True, refstarnum=None, ds9name=None):
		'''
		prepare magnitudes list for relative calibration: input magnitudes with corresponding id against the reference mags

		INPUTs:
			image_key: input image key, eg. 000.fits
		OUTPUTS:
			mags_input_matched:

		PRODUCTS:
			matched list between input magnitudes and reference magnitudes, i.e 'xxx_tpl.match'
		'''
		if self.photometry_method is None:
			self.photometry_method = raw_input("please input the photometry method used(apphot/psfphot):")

		image_key_s = image_key.split('.')[0]

		if self.photometry_method   == 'apphot':
			input_list_filename    = os.path.join(self.aperture_photometry_dir, image_key_s +self.apphot_ret_file_suffix)
			match_output_filename  = os.path.join(self.aperture_photometry_dir,image_key_s + '_tpl.match')
			photdir = self.aperture_photometry_dir
		elif self.photometry_method == 'psfphot':
                        input_list_filename    = os.path.join(self.psf_photometry_dir, image_key_s +self.psfphot_ret_file_suffix)
                        match_output_filename  = os.path.join(self.psf_photometry_dir,image_key_s + '_tpl.match')
			photdir = self.psf_photometry_dir
		else:
			raise IOError("Invalid input for photometry_method")

		if clobber:
			self.__delete_file_if_exist(match_output_filename)


		if not os.path.exists(match_output_filename):
			if not os.path.exists(input_list_filename):
				raise IOError("PSF photometry result file %s not exists"%input_list_filename)

			psfphot_ret = np.loadtxt(input_list_filename)
			xys = psfphot_ret[:,0:2]
			mags = psfphot_ret[:,2]
			magerrs = psfphot_ret[:,3]

			input_image = self.__get_internal_image(image_key, which_dir=which_dir)
			regionfile = os.path.join(self.psf_photometry_dir, image_key_s +"_xymag.reg")

			if not os.path.exists(regionfile):
				create_ds9_region_file(psfphot_ret, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='red', width=1, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)


			calmatch_table = self.__input_xys_through_ds9_get_mags(input_image, psfphot_ret, regionfile=regionfile, wantnum=refstarnum, ds9name=ds9name)
			calmatch_table.write(match_output_filename, format='ascii.fixed_width')

		calmatch_table = Table.read(match_output_filename, format='ascii.fixed_width') #the table read from previouly save table

		return calmatch_table


	def __relative_calibration_prepare(self,image_key,tpl_imgkey, remove_transient = True, tpl_obs_match_method = 'grmatch'):
		'''
		prepare magnitudes list for relative calibration: input magnitudes and reference magnitudes

		INPUTs:
			image_key: input image key, eg. 000.fits
			tpl_imgkey: reference image key, eg. 001.fits
			remove_transient: remove the transient from the matched result of reflist and input list
			tpl_obs_match_method: 'grmatch' or 'surrounding_search'

		OUTPUTS:
			mags_ref_matched:
			mags_input_matched:

		PRODUCTS:
			matched list between input magnitudes and reference magnitudes, i.e 'xxx_tpl.match'
			transform coefficient, i.e 'xxx_tpl.coef'

		Notes:
			the matched magnitudes will be filterd according to self.reference_mag_min and self.reference_mag_max

		'''
		if self.photometry_method is None:
			self.photometry_method = raw_input("please input the photometry method used(apphot/psfphot):")

		if self.photometry_method   == 'apphot':
			ref_list_filename      = os.path.join(self.aperture_photometry_dir,tpl_imgkey.split('.')[0]+self.apphot_ret_file_suffix)
			input_list_filename    = os.path.join(self.aperture_photometry_dir, image_key.split('.')[0]+self.apphot_ret_file_suffix)
			match_output_filename  = os.path.join(self.aperture_photometry_dir,image_key.split('.')[0]+ '_tpl.match')
			match_output_filename_tstrmd  = os.path.join(self.aperture_photometry_dir,image_key.split('.')[0]+ '_tpl_tstrmd.match')  # transient removed; too faint, too bright and too large uncertainty stars are rmeoved
			trans_fitting_filename = os.path.join(self.aperture_photometry_dir,image_key.split('.')[0] + '_tpl.coef')
		elif self.photometry_method == 'psfphot':
			ref_list_filename      = os.path.join(self.psf_photometry_dir,tpl_imgkey.split('.')[0]+self.psfphot_ret_file_suffix)
                        input_list_filename    = os.path.join(self.psf_photometry_dir, image_key.split('.')[0]+self.psfphot_ret_file_suffix)
                        match_output_filename  = os.path.join(self.psf_photometry_dir,image_key.split('.')[0]+ '_tpl.match')
			match_output_filename_tstrmd  = os.path.join(self.psf_photometry_dir,image_key.split('.')[0]+ '_tpl_tstrmd.match')  # transient removed; too faint, too bright and too large uncertainty stars are rmeoved
			trans_fitting_filename = os.path.join(self.psf_photometry_dir,image_key.split('.')[0] + '_tpl.coef')
		else:
			raise IOError("Invalid input for photometry_method")

		self.__delete_file_if_exist(match_output_filename)
		self.__delete_file_if_exist(trans_fitting_filename)

		if self.trim_before_grmatch:
			if self.trim_left is None or self.trim_right is None or self.trim_bottom is None or self.trim_up is None:
				raise IOError("trim boundary required!")

			ref_list_data = np.loadtxt(ref_list_filename)
			input_list_data = np.loadtxt(input_list_filename)

			print len(input_list_data)
			nl = self.trim_left
			nb = self.trim_bottom

			if self.trim_right <0 or self.trim_up <0:

				img_template = self.images[tpl_imgkey]
				hdu_tpl = fits.open(img_template)
				header_tpl = hdu_tpl['PRIMARY'].header
				try:
					NX_tpl = header_tpl['NAXIS1']
					NY_tpl = header_tpl['NAXIS2']
				except:
					print "image size not found in the tempalte image"


				img_input = self.images[image_key]
				hdu_input = fits.open(img_input)
				header_input = hdu_input['PRIMARY'].header
				try:
					NX_input = header_input['NAXIS1']
					NY_input = header_input['NAXIS2']
				except:
					print "image size not found in the input image"

				nr_tpl = NX_tpl + self.trim_right
				nu_tpl = NY_tpl + self.trim_up
				nr_input = NX_input + self.trim_right
				nu_input = NY_input + self.trim_up

			else:
				nr_tpl = self.trim_right
				nu_tpl = self.trim_up
				nr_input = self.trim_right
				nu_input = self.trim_up

			print "input image boundaries:",   nl,nr_input,nb,nu_input
			print "template image boundaries:",nl,nr_tpl,  nb,nu_tpl

			ref_list_data_trim_x = ref_list_data[np.logical_and(ref_list_data[:,0]>nl,ref_list_data[:,0]<nr_tpl),:]
			ref_list_data_trim = ref_list_data_trim_x[np.logical_and(ref_list_data_trim_x[:,1]>nb,ref_list_data_trim_x[:,1]<nu_tpl),:]

			input_list_data_trim_x = input_list_data[np.logical_and(input_list_data[:,0]>nl,input_list_data[:,0]<nr_input),:]
			input_list_data_trim = input_list_data_trim_x[np.logical_and(input_list_data_trim_x[:,1]>nb,input_list_data_trim_x[:,1]<nu_input),:]
			print len(input_list_data_trim)
			ref_list_filename = ref_list_filename + '.trim'
			input_list_filename = input_list_filename + '.trim'

			np.savetxt(ref_list_filename,  ref_list_data_trim)
			np.savetxt(input_list_filename,input_list_data_trim)

		if self.magnitude_cut_before_grmatch:
			if self.reference_mag_min is None and self.reference_mag_max is None:
				raise IOError("magnitude cut required... self.reference_mag_min or self.refence_mag_max or both are required")

			ref_list_data = np.loadtxt(ref_list_filename)
			input_list_data = np.loadtxt(input_list_filename)

			if self.reference_mag_min is not None:
				ref_list_data_magcut = ref_list_data[ref_list_data[:,2]<self.reference_mag_min]
				input_list_data_magcut = input_list_data[input_list_data[:,2]<self.reference_mag_min]
			else:
				ref_list_data_magcut = ref_list_data
				input_list_data_magcut = input_list_data

			if self.reference_mag_max is not None:
				ref_list_data_magcut = ref_list_data_magcut[ref_list_data_magcut[:,2]>self.reference_mag_max]
				input_list_data_magcut = input_list_data_magcut[input_list_data_magcut[:,2]>self.reference_mag_max]

			ref_list_filename = ref_list_filename + '.magcut'
			input_list_filename = input_list_filename + '.magcut'

			np.savetxt(ref_list_filename,  ref_list_data_magcut)
			np.savetxt(input_list_filename,input_list_data_magcut)

		visbfmatch = self.tpl_obs_before_match_display #show the data before matching
		if visbfmatch:
			refdata_plot = np.loadtxt(ref_list_filename)
			inputdata_plot = np.loadtxt(input_list_filename)

			plt.plot(refdata_plot[:,0], refdata_plot[:,1],'ro', label='reference')
			plt.plot(inputdata_plot[:,0], inputdata_plot[:,1],'bo', label='input')
			plt.xlabel('x')
			plt.ylabel('y')
			plt.legend()
			plt.show()

		if tpl_obs_match_method == 'grmatch':
			matched_table = self.__fitsh_grmatch(ref_list_filename, input_list_filename, match_output_filename, trans_fitting_filename)
		elif tpl_obs_match_method == 'surrounding_search':
			input_list_data = np.loadtxt(input_list_filename)
			ref_list_data   = np.loadtxt(ref_list_filename)
			tg_xys  = input_list_data[:,0:2]
	                ref_xys = ref_list_data[:,0:2]
	                criteria = self.tpl_obs_match_criteria
        	        target_stars,indexmat,ref_stars,indexref = self.__crude_match(tg_xys,ref_xys,criteria)

			if len(indexmat) != 0:
		                match_ret = np.hstack((ref_list_data[indexref],input_list_data[indexmat]))
        		        np.savetxt(match_output_filename,match_ret)


		match_result_plot = self.tpl_obs_match_result_display
		if match_result_plot:
			if tpl_obs_match_method == 'grmatch':
				if not os.path.exists(match_output_filename):
					raise IOError("grmatch failure...")

				if self.photometry_method == 'apphot':
					trans_output = os.path.join(self.aperture_photometry_dir, image_key.split('.')[0] + '_tpl.trans')
				elif self.photometry_method == 'psfphot':
					trans_output = os.path.join(self.psf_photometry_dir, image_key.split('.')[0] + '_tpl.trans')
				else:
					raise IOError("Invalid input for photometry_method")

				self.__fitsh_grtrans(match_output_filename, trans_output, trans_fitting_filename)
				trans_result = np.loadtxt(trans_output)
				ref_xys = trans_result[:,[0,1]]
				input_xys  = trans_result[:,[4,5]]

				plt.plot(ref_xys[:,0]-input_xys[:,0],ref_xys[:,1]-input_xys[:,1],'o')
				plt.xlabel('x_ref - x_input')
				plt.ylabel('y_ref - y_input')
				plt.show()

			elif tpl_obs_match_method == 'surrounding_search':
				if os.path.exists(match_output_filename):
					match_result = np.loadtxt(match_output_filename)
					ref_xys   = match_result[:,[0,1]]
					input_xys = match_result[:,[4,5]]

					plt.plot(ref_xys[:,0], ref_xys[:,1], 'bo')
					plt.plot(input_xys[:,0], input_xys[:,1], 'rd')
					plt.xlabel('X')
					plt.ylabel('Y')
					plt.legend(['ref','input'])
					plt.show()
				else:
					print "no valid match obtained..."
		try:
			matched_ret   = np.loadtxt(match_output_filename)

			if self.reference_mag_min is not None:
				mag_faint_cut = self.reference_mag_min
				matched_ret = matched_ret[matched_ret[:,2]<mag_faint_cut,:]

			if self.reference_mag_max is not None:
				mag_bright_cut  = self.reference_mag_max
				matched_ret = matched_ret[matched_ret[:,2]>mag_bright_cut,:]


			if self.reference_magerr_cut is not None:
				magerr_cut_ref = self.reference_magerr_cut
				matched_ret = matched_ret[matched_ret[:,3]<magerr_cut_ref,:]

			if self.input_magerr_cut is not None:
				magerr_cut_input = self.input_magerr_cut
				matched_ret = matched_ret[matched_ret[:,7]<magerr_cut_input,:]

			#get rid of the transient
			if remove_transient:
				x = self.photometry_info[tpl_imgkey]['x']
				y = self.photometry_info[tpl_imgkey]['y']

				transient_exist, corresp_xy, index_transient = self.__find_corresponding(matched_ret[:,0:2], [x,y], 2)
				if transient_exist:
					indexs = range(len(matched_ret))
					for index_remove in index_transient:
						indexs.remove(index_remove)
					matched_ret = matched_ret[indexs,:]

				np.savetxt(match_output_filename_tstrmd, matched_ret, fmt="%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f")


			ref_matched   = matched_ret[:,[0,1,2,3]]
			input_matched = matched_ret[:,[4,5,6,7]]
		except:
			self.photometry_info[image_key]['drop'] = 5
			return 0,0

		xys_ref_matched   = ref_matched[:,[0,1]]
		xys_input_matched = input_matched[:,[0,1]]

		mags_ref_matched = ref_matched[:,[2,3]]
		mags_input_matched = input_matched[:,[2,3]]


		return mags_ref_matched,mags_input_matched


	def __relative_calibration_single(self,ref_mags,input_mags,instrument_mag, xdata=None, offset_method='median', sigclipping_first  = True):
		'''

		INPUTS:
			ref_mags:
			input_mags:
			instrument_mag:
			offset_method: median, mean, or funcfit.
		'''
		if xdata is None:
			xdata = ref_mags

		if sigclipping_first:
			offsets_first = ref_mags[:,0] - input_mags[:,0]
			survive_mask,discard_mask = self.__sigma_clipping(offsets_first,variability_method='mad')
			mags_ref = ref_mags[survive_mask,0]
			mags_input = input_mags[survive_mask,0]
			magserr_ref = ref_mags[survive_mask,1]
			magserr_input = input_mags[survive_mask,1]
			xs = xdata[survive_mask,0]
			xerrs = xdata[survive_mask,1]
		else:
			mags_ref = ref_mags[:,0]
			mags_input = input_mags[:,0]
			magerr_ref = ref_mags[:,1]
			magerr_input = input_mags[:,1]

		mags_offset = mags_ref - mags_input
		magserr_offset = np.sqrt(magserr_ref**2+magserr_input**2)


		if offset_method == "funcfit":
			magoffset0 = np.median(mags_offset)
			if self.cal_offset_funcfit_type == 'o1poly':
				def func(x,a,b):
					return a*x+b
				def yerr(x,aerr,berr):
					return np.sqrt((x*aerr)**2+berr**2)
				p0 = np.array([0,magoffset0])
			elif self.cal_offset_funcfit_type == 'constant':
				def func(x,c):
					return x*0+c
				def yerr(x,cerr):
					return x*0+cerr
				p0 = np.array([magoffset0])
			else:
				raise ValueError('%s for self.cal_offset_funcfit_type not supported yet'%self.cal_offset_funcfit_type)

			fitdata, popt, perr = self.__funcfit_remove_outlier(xs,mags_offset,magserr_offset,func, p0, nsig=6, rm_mode = 'all')
			offset = func(fitdata[:,0],*popt)
			#offset_err =  yerr(fitdata[:,0],*perr)
			offset_ret = func(instrument_mag+np.mean(offset), *popt)
		elif offset_method in ['median', 'mean']:
			N = len(xs)
			fitdata = np.hstack((xs.reshape((N,1)), mags_offset.reshape((N,1)), magserr_offset.reshape((N,1))))
			if offset_method == 'median':
				offset = np.median(mags_offset)
			else:
				offset = np.mean(mags_offset)
			offset_ret = offset
		else:
			raise IOError("invalid input for offset_method...")

		cut = int(np.floor((len(fitdata)*0.688)))
		offset_dev = np.abs(fitdata[:,1]-offset)
		indice = np.argsort(offset_dev)[cut]
		offset_err = offset_dev[indice]

		if self.cal_plot:
			from matplotlib import gridspec

			fig = plt.figure(figsize=(9, 6))
			gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
			ax0 = plt.subplot(gs[0])
			plt.xlabel('reference magnitudes')
			plt.ylabel('magnitude offset (reference - input)')
			ax1 = plt.subplot(gs[1],sharey=ax0)

			ax0.errorbar(xs, mags_offset, yerr=magserr_offset, fmt='gx',label='survived')

			if sigclipping_first:
				xs_discard = xdata[discard_mask,0]
				offset_discard_mag = ref_mags[discard_mask,0] - input_mags[discard_mask,0]
				offset_discard_magerr = np.sqrt(ref_mags[discard_mask,1]**2 + input_mags[discard_mask,1]**2)
				if len(xs_discard) != 0:
					ax0.errorbar(xs_discard, offset_discard_mag, yerr=offset_discard_magerr,fmt= 'g.',label='discarded')

			ax0.errorbar(fitdata[:,0],fitdata[:,1],yerr=fitdata[:,2],fmt='ro',label='fitdata')

			if offset_method == "funcfit":
				ax0.plot(fitdata[:,0], offset,'k')
				ax0.plot(fitdata[:,0], offset+offset_err,'b',alpha=0.4)
				ax0.plot(fitdata[:,0], offset-offset_err,'b',alpha=0.4)
			else:
				ax0.plot(fitdata[:,0],np.ones(len(fitdata))*offset,'k')
				ax0.plot(fitdata[:,0],np.ones(len(fitdata))*offset+offset_err, 'b',alpha=0.4)
				ax0.plot(fitdata[:,0],np.ones(len(fitdata))*offset-offset_err, 'b',alpha=0.4)


			hist,edges = np.histogram(fitdata[:,1],bins=15)
			bin_centers = edges[:-1] + (edges[1:]-edges[:-1])/2.
			ax1.step(hist,bin_centers)

			if instrument_mag != 0:
				upper_limit = offset_ret + 5*offset_err
				lower_limit = offset_ret - 5*offset_err
				ax0.plot([instrument_mag+offset_ret,instrument_mag+offset_ret],[lower_limit,upper_limit])
			plt.show()

		relative_mag = instrument_mag + offset_ret
		relative_magerr = offset_err
		caldata_index = survive_mask  #the index of input data which are used in the relative calibration

		return relative_mag,relative_magerr, caldata_index


	def __sigma_clipping(self,input_data, sig=3, meanfunc=np.median, variability_method="stddev", verbose=0):
		"""
		Identify outliers considering the mean (if meanfunc=np.mean) or median (if meanfunc=np.median) value and 3 sigma (3*stdev),
		iterating until convergence.

		20180318: add new parameter variability_method which is used to input the choice of how to measure the variability of a univariate sample of quantitative data.
			  Two options are supported. If 'stddev', stddev =  np.sqrt(np.sum((data-meanfunc(data))**2)/last_total)
						     If 'mad', stddev = 1.48*np.median(np.abs(data-np.median(data)))

		INPUTS:
			input_data:

		"""
		data = input_data.copy()
		last_total = len(data)

		if variability_method == "stddev":
			stdev = np.sqrt(np.sum((data-meanfunc(data))**2)/last_total)
		elif variability_method == "mad":
			stdev = 1.48*np.median(np.abs(data-np.median(data)))
		else:
			raise IOError("invalid input for variability_method")

		if verbose:
			print "The stddev in the first iteration:", stdev

		diff = data - meanfunc(data)
		sfilter = np.abs(diff) < sig*stdev
		current_total = len(data[sfilter])
		# Continue iterating until convergence (no more points are removed)
		i = 1
		while last_total > current_total:
			data = data[sfilter]
			last_total = current_total

			if variability_method == "stddev":
				stdev = np.sqrt(np.sum((data-meanfunc(data))**2)/last_total)
			elif variability_method == "mad":
				stdev = 1.48*np.median(np.abs(data-np.median(data)))
			else:
				raise IOError("invalid input for variability_method")

			if verbose:
				i = i+1
				print "The stddev in iteration %s:"%i, stdev

			diff = data - meanfunc(data)
			sfilter = np.abs(diff) < sig*stdev
			current_total = len(data[sfilter])

		if variability_method == "stddev":
			stdev = np.sqrt(np.sum((data-meanfunc(data))**2)/last_total)
		elif variability_method == "mad":
			stdev = 1.48*np.median(np.abs(data-np.median(data)))
		else:
			raise IOError("invalid input for variability_method")

		meanv = meanfunc(data)
		diff = input_data - meanv
		survive_mask = np.abs(diff) < sig*stdev
		discard_mask = 	np.abs(diff) > sig*stdev

		return survive_mask,discard_mask



	def __funcfit_remove_outlier(self, x,y,yerr,func,p0, nsig=3, criteria_mode= 'individual', rm_mode = 'all', stddev_method = 'mad'):
		"""
		fitting the func to input data (x, y+/-yerr) and get the best fitting parameter.
		During the fitting process, the iterative rejection on the input data will be applied

		INPUT:
			x:
			y:
			yerr:
			func: function object takes x and parameter list and output y
			p0:
			nsig:
			criteria_mode: 'whole': the fixed single criteria for the sample; 'indivual': yerr value used for each point in criteria
			rm_mode: 'all': remove all beyond the criteria for one iteration, 'worst': only remove the worst for one iteration
			stddev_method: 'direct' or 'mad';  np.std() for direct, and 1.48*MAD() for mad
		"""

		def residuals(p, y, x,yerr):
			return (y - func(x, p))/yerr

		def get_stddev(y, stddev_method):
			if stddev_method == 'direct':
				stddev = np.std(y)
			elif stddev_method == 'mad':
				stddev = np.median(np.abs(y-np.median(y)))
			else:
				raise ValueError("invalid input for stddev_method...")
			return stddev

		def remove_outliers(x, y, yerr, abs_diffs, nsig, stddev, rm_mode, criteria_mode):
			if criteria_mode == 'whole':
				goodcriteria = nsig*stddev
			elif criteria_mode == 'individual':
				goodcriteria = nsig*yerr
			else:
				raise ValueError('invalid input for criteria_mode')
			if rm_mode == 'all':
				goodflt = abs_diffs < goodcriteria
				x    = x[goodflt]
				y    = y[goodflt]
				yerr = yerr[goodflt]
			elif rm_mode == 'worst':
				indice = np.argsort(abs_diffs)
				if abs_diffs[indice[-1]] > goodcriteria:
					goodflt = indice[:-1]
					x    = x[goodflt]
					y    = y[goodflt]
					yerr = yerr[goodflt]
			else:
				raise ValueError('invalid input for rm_mode...')
			return x,y,yerr

		popt, pcov = curve_fit(func, x, y, p0=p0, sigma=yerr)
		p0 = popt
		perr = np.sqrt(np.diag(pcov))
		y_model = func(x,*popt)
		abs_diffs = np.abs(y-y_model)
		last_total = len(y)

		stddev = get_stddev(y, stddev_method)
		x,y,yerr = remove_outliers(x,y,yerr,abs_diffs,nsig,stddev, rm_mode, criteria_mode)
		current_total = len(y)

		while last_total > current_total: #Continue iterating until convergence (no more points are removed)
			last_total = current_total
			popt, pcov = curve_fit(func, x, y, p0=p0, sigma=yerr)
			p0 = popt
			perr = np.sqrt(np.diag(pcov))
			y_model = func(x,*popt)
			abs_diffs = np.abs(y-y_model)

			stddev = get_stddev(y, stddev_method)
			x,y,yerr = remove_outliers(x,y,yerr,abs_diffs,nsig,stddev, rm_mode, criteria_mode)
			current_total = len(y)

		gooddata = np.hstack((x.reshape((len(x),1)), y.reshape((len(y),1)), yerr.reshape((len(yerr),1))))
		return gooddata,popt,perr



	def __crude_match(self,array1,array2,criteria,verbose=False):
	    '''
	    Conduct crude match between two given list of positions of stars

	    INPUTS:
	    array1: the first star array with first and second column array1[:,0:2] are star position(x,y) and array1.shape = (N1,m)
	    array2: the second star array with first and second column array1[:,0:2] are star position(x,y)and array1.shape = (N2,m)
	            array1 and array have the same width(columns) and ONE of N1 and N2 can equal 1
	    criteria: = c(PIXES) the criteria for match is that one star is within c pixels from the corrresponding star's central position
	    '''

	    if len(array1) < len(array2):
	        refarray = array2
	        matarray = array1
	        alter = 1
	    else:
	        refarray = array1
	        matarray = array2
	        alter = 0

	    N1 = len(matarray)
	    N2 = len(refarray)

	    matmask=[]
	    refmask=[]
	    if matarray.shape == (2,) or matarray.shape == (1,2):
	        matarray = matarray.copy().ravel()
	        diffarray = refarray - matarray
	        temp = []
	        for j,diff in enumerate(diffarray):
	            if np.abs(diff[0])<criteria and np.abs(diff[1])< criteria:
	                temp.append(j)
	                #print diff

	        if len(temp)==1:
	            matched_ref = refarray[temp[0]]
	            matched_mat = matarray
	            matmask.append(1)
	            refmask.append(temp[0])
	        else:
	            if len(temp)>1:
	                print 'more than one objects fall into the criteria region of reference star. This object is droped!'
	            stars1 = []
	            stars2 = []
	            mask1 =[]
	            mask2 =[]
	            return stars1,mask1,stars2,mask2
	    else:
	        for i,current in enumerate(matarray):
	            diffarray = refarray - current
	            if verbose:
	                print diffarray
	            temp = []
	            for j,diff in enumerate(diffarray):
	                try:
	                    diff1 = diff[0]
	                    diff2 = diff[1]
	                except:
	                    diff1 = diff[0,0]
	                    diff2 = diff[0,1]

	                if np.abs(diff1)<criteria and np.abs(diff2)< criteria:
	                    temp.append(j)
	                    #print diff
	            if len(temp)==1:
	                matmask.append(i)
	                refmask.append(temp[0])
	            elif len(temp)>1:
	                print 'more than one objects fall into the criteria region of reference star. This object is droped!'

	        if verbose:
	            print refmask
	            print matmask

	        matched_ref = refarray[refmask]
	        matched_mat = matarray[matmask]

	    #make sure the output stars1,stars2 correspond to array1,array2 respectively

	    if alter:
	        stars1 = matched_mat
	        stars2 = matched_ref
	        mask1 = matmask
	        mask2 = refmask
	    else:
	        stars1 = matched_ref
	        stars2 = matched_mat
	        mask1 = refmask
	        mask2 = matmask

	    if verbose:
	        print "From array1:\n"
	        print stars1
	        print "From array2:\n"
	        print stars2

	    return stars1,mask1,stars2,mask2



	def __find_corresponding(self,base_array,target,criteria,verbose=False):
		'''
		find the location of 'target' in 'base_array' when it meet the 'criteria'
		Input:
			base_array: (N,2) array
			target:     two element array or list
			criteria: scaler, the distance to the target in both x and y direction

		Output:
			exist: True or False
			base_array[ret_mask]: (N,2) shape array and N can be 1
			ret_mask: list, indice of rows meeting requirement
		'''
		if isinstance(target,list):
			target = np.array(target)

		target = target.reshape((1,2))
		diffs = base_array-target
		ret_mask= []
		for i,diff in enumerate(diffs):
		    if np.abs(diff[0])<criteria and np.abs(diff[1])<criteria:
			ret_mask.append(i)
		if ret_mask == []:
		    exist = False
		    return exist,0,0
		else:
		    exist = True
		    return exist,base_array[ret_mask],ret_mask


	def psf_photometry_dophot(self, which_dir = 'raw_image'):
		'''
		See self.__psf_photometry_dophot_flt
		'''
		self.__find_template_imagekey()
		for flt in self.templates.keys():
			self.__psf_photometry_dophot_flt(flt, which_dir = which_dir)
			self.__psf_photometry_dophot_flt_target(flt)



	def __psf_photometry_dophot_flt(self,flt, which_dir = 'raw_image'):
		'''
		See self.__psf_photometry_dophot_single
		'''
		self.__dict2table()
		flt_table = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)
		for image_key in flt_table['name']:
			print image_key
			if self.photometry_info[image_key]['drop'] > 0:
				continue
			psfret_savefile = os.path.join(self.psf_photometry_dir,image_key.split('.')[0]+self.psfphot_ret_file_suffix)
			if not os.path.isfile(psfret_savefile) or  self.renew_psf_photometry_retfile:
				self.__psf_photometry_dophot_single(image_key,psfret_savefile,output_residuals=True, which_dir=which_dir)



	def __psf_photometry_dophot_flt_target(self,flt):
		'''
		See self.__psf_photometry_dophot_single_target
		'''
		self.__dict2table()
		flt_table = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)
		for image_key in flt_table['name']:
			if self.photometry_info[image_key]['drop'] > 0:
				continue
			if self.photometry_info[image_key]['instmag'] != 99.99 and not self.renew_psf_photometry:
				continue
			match_criteria = self.psfphot_sn_match_criteria
			self.__psf_photometry_dophot_single_target(image_key,match_criteria = match_criteria)


	def __psf_photometry_dophot_single_target_ds9_pick(self, image_key, match_criteria=3, which_dir='raw_image',  verbose=1 , update_xy=True, newds9=True, box_length=100, box_xc=None,box_yc=None):
		'''
		this is created for the initial purpose of SMARTS-NIR photometry where the source can be very sparse on the image. The PSF photometry is obtained first then select the target mag from the
		magnitude list.
		'''
		imgkey_s = image_key.split('.')[0]
		psfphot_ret_file = os.path.join(self.psf_photometry_dir,imgkey_s+self.psfphot_ret_file_suffix)
		if not os.path.exists(psfphot_ret_file):
			print "PSF photometry result file %s not exists"%psfphot_ret_file
			self.photometry_info[image_key]['drop'] = 4
		else:
			psfphot_ret = np.loadtxt(psfphot_ret_file)
			xys = psfphot_ret[:,0:2]
			mags = psfphot_ret[:,2]
			magerrs = psfphot_ret[:,3]

			input_image = self.__get_internal_image(image_key, which_dir=which_dir)
			regionfile = os.path.join(self.psf_photometry_dir, imgkey_s+"_xymag.reg")

			if (box_xc is not None) and (box_yc is not None):
				tmpmask = (psfphot_ret[:,0]>(box_xc-box_length/2))*(psfphot_ret[:,0]<(box_xc+box_length/2))*(psfphot_ret[:,1]<(box_yc+box_length/2))*(psfphot_ret[:,1]>(box_yc-box_length/2))
				psfphot_reg = psfphot_ret[tmpmask,:]
			else:
				psfphot_reg = psfphot_ret

			create_ds9_region_file(psfphot_reg, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='red', width=1, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)

                	xy = self.__get_xy_on_image_from_ds9(input_image=input_image, regionfile=regionfile, newds9=newds9)
                	print "xy from mouse pick up:", xy

                	x = xy[0] - 1
                	y = xy[1] - 1
                	xy_target = [x,y]
                	yesfound, xy_match,index = self.__find_corresponding(xys, xy_target, match_criteria)

                	if (not yesfound):
                	        print "no object found within criteeia..."
                	        mag = 99.99
                	        magerr = 99.99
				sigdrop = 4
                	elif len(index)>1:
                	        print 'more than one objects fall into the criteria region of reference star. This object is droped!'
                	        mag = 99.99
                	        magerr = 99.99
				sigdrop = 4
                	else:
                	        print "the corresponding object found at (%s)"%xy_match

				if update_xy:
					snx = xy_match[0][0]
					sny = xy_match[0][1]
					print "the (x,y) coordinate will be updated to (%s,%s)"%(snx, sny)
					self.photometry_info[image_key]['x'] = snx
					self.photometry_info[image_key]['y'] = sny

                	        mag = psfphot_ret[index,2]
                	        magerr = psfphot_ret[index,3]
				sigdrop = 0
				if verbose:
					print "the mag: %s+/-%s"%(mag, magerr)

			self.photometry_info[image_key]['instmag'] = mag
			self.photometry_info[image_key]['instmagerr'] = magerr
			self.photometry_info[image_key]['drop'] = sigdrop




	def __psf_photometry_dophot_single_target(self,image_key,ds9_display = False,match_criteria = 3, drop_nstar_less_than_this = 3):
		'''
		get the mag for the target from mag list obtained from PSF photometry
		'''
		psfphot_ret_file = os.path.join(self.psf_photometry_dir,image_key.split('.')[0]+self.psfphot_ret_file_suffix)
		psfphot_ret = np.loadtxt(psfphot_ret_file)

		if len(psfphot_ret) < drop_nstar_less_than_this:
			self.photometry_info[image_key]['drop'] = 8
			return

		xys = psfphot_ret[:,0:2]
		mags = psfphot_ret[:,2]
		magerrs = psfphot_ret[:,3]
		x = self.photometry_info[image_key]['x']
		y = self.photometry_info[image_key]['y']
		xy = [x-1, y-1]

		exist,xy_match,indice_match = self.__find_corresponding(xys, xy, match_criteria)

		if exist:

			if len(indice_match)>1:
				print "more than one source detected at the SN position... drop this one"
				self.photometry_info[image_key]['drop'] = 9
			else:
				self.photometry_info[image_key]['instmag'] = mags[indice_match]
				self.photometry_info[image_key]['instmagerr'] = magerrs[indice_match]

			if ds9_display:
				input_image = self.images[image_key]
				d = self.__display_image_with_ds9(input_image, newds9=True)
				print "Red circle indicate the corresponding position from PSF photometry results"
				if len(indice_match) >1:
					for xy_reg in xy_match:
						self.__load_source_regions_to_ds9(d, xy_reg, radius=5, color='red', width=1)
				else:
					xy_reg = xy_match.ravel()
					self.__load_source_regions_to_ds9(d, xy_reg, radius=5, color='red', width=1)
				print "Green circle indicate the SN position previously provided"
				xy_reg = [x,y]
				self.__load_source_regions_to_ds9(d, xy_reg, radius=5, color='green',width=1)
		else:
			if ds9_display:
				input_image = self.images[image_key]
				d = self.__display_image_with_ds9(input_image, newds9=True)
				print "Green circle indicate the SN position previously provided"
				xy_reg = [x,y]
				self.__load_source_regions_to_ds9(d, xy_reg, radius=5, color='green',width=1)
			print "photometry result not found for %s at position %s"%(image_key,xy)
			self.photometry_info[image_key]['drop'] = 4



	def __psf_photometry_dophot_single(self,image_key,ret_save_filename=None, output_residuals=False, which_dir ='raw_image', pmfile=None, renew_pmfile=True, verbose=1):
		'''
		perform PSF photometry with DoPhot package

		INPUTS:
			image_key:
			ret_save_filename: the output photometry file containing the well formated x,y,mag,magerr
			output_residuals: save the residual image after model subtraction?
			which_dir: where the input image is from
			pmfile: the input parameter file; None if not specified
			renew_pmfile: you can provide explicitly the parameter file under /path/to/sndir/warehouse/psfphotometry/xxx.pm, or the parameter file can be from last run. The above existing parameter file will be replaced if renew_pmfile is True

		'''
		if verbose:
			print image_key

		img_key_s  = image_key.split('.')[0]

		if pmfile is not None: #if parameter file is specified through input variable, the pmfile will be used
			if not os.path.exists(pmfile):
				raise ValueError("the input parameter file %s not exists"%pmfile)
			pmfile_external = True
		else:
			pmfile_tmp = img_key_s + '.pm'
			pmfile_external = False
			cwd = os.getcwd()
			cwdod = os.path.join(cwd, self.current_sn)#oject dir under cwd
			cwdodid = os.path.join(cwdod, self.current_telescope)#instrument dir under cwdod
			if self.insulate_pmfile:
				if not os.path.exists(cwdod):
					os.mkdir(cwdod)
				if not os.path.exists(cwdodid):
					os.mkdir(cwdodid)
				pmfile  = os.path.join(self.current_sn, os.path.join(self.current_telescope, pmfile_tmp))
				photimg = os.path.join(self.current_sn, os.path.join(self.current_telescope, image_key))
			else:
				pmfile = pmfile_tmp
				photimg = image_key

			pmfile_src = os.path.join(cwd,pmfile)  #this is the pm file which will be the direct input of dophot command execution
			pmfile_dst = os.path.join(self.psf_photometry_dir,pmfile_tmp) #this is the archived pm file for specified image stored under event specific directory

			if os.path.exists(pmfile_dst) and (not renew_pmfile):
				self.__delete_file_if_exist(pmfile_src)
				shutil.copy(pmfile_dst, pmfile_src)
				savepmfile = False #need to save/archive the newly created pmfile
			else:
				savepmfile = True
				if self.dophot_version == 'fortran':
					dophot_fortran_pm = self.__prepare_pm_data_dophot_fortran(image_key, photimg)
					self.__prepare_pm_file_dophot_fortran(pmfile, **dophot_fortran_pm)
				elif self.dophot_version == 'C':
					dophot_C_pm = self.__prepare_pm_data_dophot_C(image_key, photimg)
					self.__prepare_pm_file_dophot_C(pmfile, **dophot_C_pm)
				else:
					raise IOError("dophot version %s is not supported..."%self.dophot_version)
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		input_image_here = os.path.join(os.getcwd(), photimg)
		self.__delete_file_if_exist(input_image_here)

		if self.dophot_image_prepare == 'softlink':
			ln_command = "ln -s %s %s"%(input_image,input_image_here)
			if verbose:
				print ln_command
			os.system(ln_command)
		elif self.dophot_image_prepare == 'imcopy':
			imcopy_iraf(input_image+'[0]', input_image_here)
			if verbose:
				print "imcopy %s[0] %s"%(input_image, input_image_here)

		if self.dophot_version == 'fortran':
			dophot = os.path.join(self.base.dophot_fortran_dir, 'dophot')
		elif self.dophot_version == 'C':
			dophot = os.path.join(self.base.dophot_C_dir, 'dophot')
		else:
			raise IOError("dophot version %s is not supported..."%self.dophot_version)

		if verbose:
			psfphot_command = '%s %s'%(dophot,pmfile)
		else:
			psfphot_command = '%s %s >/dev/null'%(dophot,pmfile)

		print psfphot_command
		os.system(psfphot_command)

		if not pmfile_external:
			if savepmfile:
				shutil.move(pmfile_src,pmfile_dst)

		self.__clean_dophot_output_files(image_key)
		os.remove(input_image_here)

		ret_good_small = self.__extract_dophot_phot_results(image_key)
		if ret_save_filename is None:
			ret_save_filename = os.path.join(self.psf_photometry_dir,image_key.split('.')[0]+self.psfphot_ret_file_suffix)
		np.savetxt(ret_save_filename,ret_good_small)


	def __extract_dophot_phot_results(self, image_key):
		'''
		extract and output the clean result in format of x,y,mag,magerr
		'''
		img_key_s = image_key.split('.')[0]
		output_obj = img_key_s + '.out'	#the psf photometry result directly from the dophot  with COMPLETE data format
		psfphot_ret_dst = os.path.join(self.psf_photometry_dir,output_obj)
		#output result file of dophot photometry with complete format has 15 columns; but elements in some rows stick toghther
		#we need to eliminate wrong rows
		temp_file = os.path.join(self.psf_photometry_dir, img_key_s + '_temp.out')
		fid_out = open(temp_file, 'w')
		renew_required = False
		for line in open(psfphot_ret_dst).readlines():
			line_segs = line.split()
			if len(line_segs) < 15:
				renew_required = True
				continue
			fid_out.write(line)
		fid_out.close()

		if renew_required:
			shutil.move(temp_file, psfphot_ret_dst)

		psfret_whole = np.loadtxt(psfphot_ret_dst)
		mags = psfret_whole[:,4]
		mask_good = np.where(mags != -99.999)[0]
		ret_good  = psfret_whole[mask_good,:]

		if self.dophot_select_stars_method == 'eq': 	#source type identation 1(star); see docs/dophot_notes.txt for summary of object types
			ret_good_small = ret_good[ret_good[:,1]==self.dophot_select_star_value]	 #only x,y,mag,magerr
		elif self.dophot_select_stars_method == 'lt':
			ret_good_small = ret_good[ret_good[:,1]<self.dophot_select_star_value]
		elif self.dophot_select_stars_method == 'gt':
			ret_good_small = ret_good[ret_good[:,1]>self.dophot_select_star_value]
		else:
			raise ValueError("invalid input for self.dophot_select_stars_method")

		xymags = ret_good_small[:,[2,3,4,5]]
		#save the source regions as a ds9 region file,
		phot_source_reg =  os.path.join(self.psf_photometry_dir, image_key.split('.')[0] + '_all.reg')
		create_ds9_region_file(ret_good, phot_source_reg, clobber = True,  x_col=2, y_col=3, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		#save the sources as a ds9  reg file, only sources meet the requirement of identification and photometry quality
		phot_source_reg =  os.path.join(self.psf_photometry_dir, image_key.split('.')[0] + '.reg')
		create_ds9_region_file(xymags, phot_source_reg, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		return xymags

	def __clean_dophot_output_files(self, image_key):
		'''
		Clean the miscellaneous output from dophot task
		'''
		img_key_s = image_key.split('.')[0]
		finish_file   = img_key_s +'.finish'
		finish_src = os.path.join(os.getcwd(),finish_file)
		finish_dst = os.path.join(self.psf_photometry_dir,finish_file)
		if os.path.exists(finish_src):
			shutil.move(finish_src,finish_dst)

		psf_file   = img_key_s +'_psf.fits'
		psffile_src = os.path.join(os.getcwd(),psf_file)
		psffile_dst = os.path.join(self.psf_photometry_dir,psf_file)
		if os.path.exists(psffile_src):
			shutil.move(psffile_src,psffile_dst)

		logfile    = img_key_s +'_dophot.log'
		logfile_src = os.path.join(os.getcwd(),logfile)
		logfile_dst = os.path.join(self.psf_photometry_dir,logfile)
		if os.path.exists(logfile_src):
			shutil.move(logfile_src,logfile_dst)

		pmfile_out = img_key_s +'out.pm'
		pmfile_out_src = os.path.join(os.getcwd(),pmfile_out)
		pmfile_out_dst = os.path.join(self.psf_photometry_dir,pmfile_out)
		shutil.move(pmfile_out_src,pmfile_out_dst)

		output_obj = img_key_s + '.out'	#the psf photometry result directly from the dophot  with COMPLETE data format
		output_image = 'out_' + image_key
		output_obj_src = os.path.join(os.getcwd(), output_obj)
		psfphot_ret_dst = os.path.join(self.psf_photometry_dir,output_obj)
		shutil.move(output_obj_src,psfphot_ret_dst)

		output_img_src = os.path.join(os.getcwd(), output_image)
		psfphot_residual_img_dst = os.path.join(self.psf_photometry_dir,output_image)
		if os.path.isfile(output_img_src):
			shutil.move(output_img_src,psfphot_residual_img_dst)


	def __prepare_pm_data_dophot_fortran(self, image_key, photimg):
		'''
		prepare the input data for the pm(parameter modification) file
		'''
		img_key_s = image_key.split('.')[0]

		fwhm = self.photometry_info[image_key]['fwhm']
		skyvalue = self.photometry_info[image_key]['bkg']
		eperdn   = float(self.base.telescopes_info[self.current_telescope]['gain'])
		rdnoise  = float(self.base.telescopes_info[self.current_telescope]['readnoise'])
		output_obj = img_key_s + '.out'	#the psf photometry result directly from the dophot  with COMPLETE data format
		output_image = 'out_' + image_key
		pmfile_out = img_key_s +'out.pm'
		finishfile = img_key_s + '.finish'

		autopmdata_dict = {}
		autopmdata_dict['FWHM']    = fwhm
		autopmdata_dict['SKY']     = skyvalue
		autopmdata_dict['EPERDN']  = eperdn
		autopmdata_dict['RDNOISE'] = rdnoise
		autopmdata_dict['IMAGE_IN'] = photimg
		autopmdata_dict['IMAGE_OUT']=output_image
		autopmdata_dict['OBJECTS_OUT']=output_obj
		autopmdata_dict['PARAMS_OUT']=pmfile_out
		autopmdata_dict['FINISHFILE'] =finishfile

		dophot_fortran_pm = OrderedDict()
		for kw in self.dophot_fortran_pm.keys():
			dophot_fortran_pm[kw] = self.dophot_fortran_pm[kw]

		for kw in autopmdata_dict.keys():
			if kw not in dophot_fortran_pm:
				dophot_fortran_pm[kw] = autopmdata_dict[kw]
			elif dophot_fortran_pm[kw] is None:
				dophot_fortran_pm[kw] = autopmdata_dict[kw]
			else:
				print "%s explicitly set"%pm
		return dophot_fortran_pm


	def __prepare_pm_file_dophot_fortran(self,pm_file, **dophot_fortran_pm):
		'''
		prepare dophot input file
		'''
        	fid_out = open(pm_file,'wt')
		for kw in dophot_fortran_pm.keys():
			if dophot_fortran_pm[kw] is None:
				continue
			elif isinstance(dophot_fortran_pm[kw],str):
				fid_out.write("%s = '%s' \n"%(kw, dophot_fortran_pm[kw]))
			else:
				fid_out.write("%s = %s \n"%(kw, str(dophot_fortran_pm[kw])))
		fid_out.write('END')
		fid_out.close()


	def __prepare_pm_data_dophot_C(self, image_key, photimg):
		'''
		prepare the input data for the pm(parameter modification) file
		'''
		img_key_s = image_key.split('.')[0]

		fwhm = self.photometry_info[image_key]['fwhm']
		skyvalue = self.photometry_info[image_key]['bkg']
		eperdn = float(self.base.telescopes_info[self.current_telescope]['gain'])
		rdnoise = float(self.base.telescopes_info[self.current_telescope]['readnoise'])
		saturation_value = float(self.base.telescopes_info[self.current_telescope]['satval'])

		output_obj = img_key_s + '.out'	#the psf photometry result directly from the dophot  with COMPLETE data format
		output_image = 'out_' + image_key
		pmfile_out = img_key_s +'out.pm'
		logfile    = img_key_s +'_dophot.log'
		psf_file   = img_key_s +'_psf.fits'
		centintmax = saturation_value

		autopmdata_dict = {}
		autopmdata_dict['FWHM']    = fwhm
		autopmdata_dict['SKY']     = skyvalue
		autopmdata_dict['EPERDN']  = eperdn
		autopmdata_dict['RDNOISE'] = rdnoise
		autopmdata_dict['IMAGE_IN'] = photimg
		autopmdata_dict['IMAGE_OUT']=output_image
		autopmdata_dict['OBJECTS_IN'] = ' '
		autopmdata_dict['OBJECTS_OUT']=output_obj
		autopmdata_dict['PARAMS_OUT'] = pmfile_out
		autopmdata_dict['EMP_SUBRAS_OUT'] = psf_file
		autopmdata_dict['SHADOWFILE_OUT'] = ' '
		autopmdata_dict['ERRORS_OUT']     = ' '
		autopmdata_dict['LOGFILE'] = logfile

		dophot_C_pm = OrderedDict()
		for kw in self.dophot_C_pm.keys():
			dophot_C_pm[kw] = self.dophot_C_pm[kw]

		for kw in autopmdata_dict.keys():
			if kw not in dophot_C_pm:
				dophot_C_pm[kw] = autopmdata_dict[kw]
			elif dophot_C_pm[kw] is None:
				dophot_C_pm[kw] = autopmdata_dict[kw]
			else:
				print "%s explicitly set"%kw

		return dophot_C_pm


	def __prepare_pm_file_dophot_C(self,pm_file, **dophot_C_pm):
		'''
		prepare pm file
		'''
		fid_out = open(pm_file,'wt')
		for kw in dophot_C_pm.keys():
			if dophot_C_pm[kw] is None:
				continue
			elif isinstance(dophot_C_pm[kw],str):
				fid_out.write("%s = '%s' \n"%(kw, dophot_C_pm[kw]))
			else:
				fid_out.write("%s = %s \n"%(kw, str(dophot_C_pm[kw])))
		fid_out.write('END')
		fid_out.close()



	def get_apphot_iraf_parameters(self, image_key, saveparfile = True):
		'''
		Get the default config parameters for apphot.daophot
		'''
		imgkey = image_key.split('.')[0]
		parfile = os.path.join(self.aperture_photometry_dir, imgkey+'.par')

		config = ConfigParser.ConfigParser()
		if os.path.exists(parfile) and (not self.renew_apphot_parfile):
			with open(parfile,'rw') as cfgfile:
				config.readfp(cfgfile)
				app     = config.getfloat('param','app')
				skyin   = config.getfloat('param','skyin')
				skywidth  = config.getfloat('param','skywidth')
				fwhmpsf = config.getfloat('param','fwhmpsf')
				sky_sigma_down  = config.getfloat('param','sky_sigma_down')
				sky_sigma_up    = config.getfloat('param','sky_sigma_up')
				def_zeropt      = config.getfloat('param','def_zeropt')

				self.apphot_iraf_options['fwhmpsf']       = fwhmpsf
				self.apphot_iraf_options['app']           = app
				self.apphot_iraf_options['skyin']         = skyin
				self.apphot_iraf_options['skywidth']      = skywidth
				self.apphot_iraf_options['sky_sigma_down']= sky_sigma_down
				self.apphot_iraf_options['sky_sigma_up']  = sky_sigma_up
				self.apphot_iraf_options['def_zeropt']    = def_zeropt

		else:
			fwhmpsf = self.photometry_info[image_key]['fwhm']
			self.apphot_iraf_options['fwhmpsf'] = fwhmpsf

			if self.apphot_iraf_autoscale_aperture:
				self.__autoscale_apphot_iraf_aperture_pars(fwhmpsf)

			if saveparfile:
				self.__delete_file_if_exist(parfile)
				config.add_section('param')
				config.set('param', 'fwhmpsf', self.apphot_iraf_options['fwhmpsf'])
				config.set('param', 'app', self.apphot_iraf_options['app'])
				config.set('param', 'skyin', self.apphot_iraf_options['skyin'])
				config.set('param', 'skywidth', self.apphot_iraf_options['skywidth'])
				config.set('param', 'sky_sigma_down', self.apphot_iraf_options['sky_sigma_down'])
				config.set('param', 'sky_sigma_up', self.apphot_iraf_options['sky_sigma_up'])
				config.set('param', 'def_zeropt', self.apphot_iraf_options['def_zeropt'])

				fid = open(parfile,'w')
				config.write(fid)
				fid.close()


	def __autoscale_apphot_iraf_aperture_pars(self, fwhmpsf):
		'''
		Auto-scale the aperture sizes (source aperture and background region)
		'''
		if fwhmpsf <=0 or fwhmpsf>50: #check to prevent radicuous results
			raise ValueError("The star profile with negative fwhm or fwhm value larger than 50 pixels. Really??")

		self.apphot_iraf_options['app']     = self.apphot_iraf_app_nfwhm * fwhmpsf
		self.apphot_iraf_options['skyin']   = self.apphot_iraf_skyin_nfwhm * fwhmpsf
		self.apphot_iraf_options['skywidth']= self.apphot_iraf_skywidth_nfwhm * fwhmpsf



	def aperture_photometry_apphot_iraf(self,flts = None, which_dir = 'raw_image', aperture_size_fixed = False, centering_target = True, centering_ref_stars = True, updatetable=1):
		'''
		Aperture photometry on images with filter given in 'flts', which default is all filters for
		which the template image is given.

		Two functions involved in for this purpose:
		self.__aperture_photometry_apphot_iraf(flt): aperture photometry for all given objects in the image

		self.__aperture_photometry_apphot_iraf_target(flt): aperture photometry for target


		Inputs:
			flts: default None; allowed input for example ['B','V','R','I']
			which_dir: perform aperture photometry on images stored in self.raw_image_dir or self.modified_image_dir, indicated by
			'raw_image' or 'modified_image'
		'''
		if updatetable:
			self.__dict2table()
		self.__find_template_imagekey(updatetable=updatetable)

		if flts is None:
			flts = self.templates.keys()

		for flt in flts:
			print "Working on %s band"%flt
			self.__aperture_photometry_apphot_iraf(flt, which_dir = which_dir, aperture_size_fixed=aperture_size_fixed, centering = centering_ref_stars)
			self.__aperture_photometry_apphot_iraf_target(flt, which_dir = which_dir, aperture_size_fixed=aperture_size_fixed, centering = centering_target, updatetable=updatetable)



	def __aperture_photometry_apphot_iraf_target(self,flt, which_dir = 'raw_image', aperture_size_fixed = False, imagefwhm=None, centering=True, x=None, y=None, updatetable=1):
		'''
		pipeline script for aperture photometry on images with filter 'flt'
		which_dir: 'raw_image', 'modified_image', or 'subtracted_image'
		See also:
			self.__aperture_photometry_apphot_iraf_target_single
		'''
		if updatetable:
			self.__dict2table()
		info_flt = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)
		images = info_flt['name']

		for image_key in images:
			if self.photometry_info[image_key]['drop'] >0:
				continue
			if self.photometry_info[image_key]['instmag'] != 99.99 and not self.renew_aperture_photometry:
				continue
			photret_singlept = self.__aperture_photometry_apphot_iraf_target_single(image_key, which_dir =  which_dir, aperture_size_fixed=aperture_size_fixed, fwhm_img=imagefwhm, centering = centering, x=x, y=y)

			if len(photret_singlept) == 0:
				self.photometry_info[image_key]['drop'] = 4
				continue

			mag = photret_singlept[0,7]
			magerr = photret_singlept[0,8]

			self.photometry_info[image_key]['instmag'] = mag
			self.photometry_info[image_key]['instmagerr'] = magerr



	def __aperture_photometry_apphot_iraf_target_single(self, image_key, which_dir = 'raw_image', aperture_size_fixed = False, fwhm_img=None, centering= True, x=None, y=None):
		'''
                For each image: get target xy; get apphot input parameters; perform aperture photometry
		Involved functions:
                        self.get_apphot_iraf_parameters
                        self.__aperture_photometry_apphot_iraf_single_image
		'''
		print image_key
		if x is None:
			x = self.photometry_info[image_key]['x'] - 1
		if y is None:
			y = self.photometry_info[image_key]['y'] - 1

                if x == -1 or y==-1:
                        raise ValueError("no valid souce position is given...")

		xy = np.array([x,y]).reshape((1,2))
		image = self.__get_internal_image(image_key, which_dir=which_dir)
		self.get_apphot_iraf_parameters(image_key, saveparfile=True)
		options = self.apphot_iraf_options.copy()

		print options
		photret_singlept = self.__aperture_photometry_apphot_iraf_single_image(image,xy,options,centering = centering)
		#xpos, ypos, flux, bkg,bkg_stdev, sky_counts, target_area,mag,mag_err

		photoret_singlept = photret_singlept.ravel()
		print photret_singlept

		return photret_singlept


	def __aperture_photometry_single_target_ds9_pick(self, image_key, match_criteria=3, which_dir='raw_image',  verbose=1 , update_xy=True, newds9=True, box_length=100, box_xc=None,box_yc=None):
		'''
		this is created for the initial purpose of SMARTS-NIR photometry where the source can be very sparse on the image. The PSF photometry is obtained first then select the target mag from the
		magnitude list.
		'''
		imgkey_s = image_key.split('.')[0]
		apphot_ret_file = os.path.join(self.aperture_photometry_dir,imgkey_s+self.apphot_ret_file_suffix)
		if not os.path.exists(apphot_ret_file):
			print "aperture photometry result file %s not exists"%apphot_ret_file
			self.photometry_info[image_key]['drop'] = 4
		else:

			apphot_ret = np.loadtxt(apphot_ret_file)
			xys = apphot_ret[:,0:2]
			mags = apphot_ret[:,2]
			magerrs = apphot_ret[:,3]

			input_image = self.__get_internal_image(image_key, which_dir=which_dir)
			regionfile = os.path.join(self.aperture_photometry_dir, imgkey_s+"_xymag.reg")

			if (box_xc is not None) and (box_yc is not None):
				tmpmask = (apphot_ret[:,0]>(box_xc-box_length/2))*(apphot_ret[:,0]<(box_xc+box_length/2))*(apphot_ret[:,1]<(box_yc+box_length/2))*(apphot_ret[:,1]>(box_yc-box_length/2))
				apphot_reg = apphot_ret[tmpmask,:]
			else:
				apphot_reg = apphot_ret

			create_ds9_region_file(apphot_reg, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='red', width=1, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)

                	xy = self.__get_xy_on_image_from_ds9(input_image=input_image, regionfile=regionfile, newds9=newds9)
                	print "xy from mouse pick up:", xy

                	x = xy[0] - 1
                	y = xy[1] - 1
                	xy_target = [x,y]
                	yesfound, xy_match,index = self.__find_corresponding(xys, xy_target, match_criteria)

                	if (not yesfound):
                	        print "no object found within criteeia..."
                	        mag = 99.99
                	        magerr = 99.99
				sigdrop = 4
                	elif len(index)>1:
                	        print 'more than one objects fall into the criteria region of reference star. This object is droped!'
                	        mag = 99.99
                	        magerr = 99.99
				sigdrop = 4
                	else:
                	        print "the corresponding object found at (%s)"%xy_match

				if update_xy:
					snx = xy_match[0][0]
					sny = xy_match[0][1]
					print "the (x,y) coordinate will be updated to (%s,%s)"%(snx, sny)
					self.photometry_info[image_key]['x'] = snx
					self.photometry_info[image_key]['y'] = sny

                	        mag = apphot_ret[index,2]
                	        magerr = apphot_ret[index,3]
				sigdrop = 0
				if verbose:
					print "the mag: %s+/-%s"%(mag, magerr)

			self.photometry_info[image_key]['instmag'] = mag
			self.photometry_info[image_key]['instmagerr'] = magerr
			self.photometry_info[image_key]['drop'] = sigdrop



	def aperture_photometry_apphot_iraf_single_image(self, image_key, overwrite=False, which_dir='modified_image', aperture_size_fixed=False, centering=True):
		'''
		Aperture photometry for individual image. This is not intended to be used in pack mode (not for pipeline usage)
		'''
		photret_file = image_key.split('.')[0]  + self.apphot_ret_file_suffix
		photret_file_abs = os.path.join(self.aperture_photometry_dir,photret_file)

		if os.path.exists(photret_file_abs):
			if not overwrite:
				print "photometry result file for %s exists: %s; please specify overwrite=True if you want to perform aperture photometry again"%(image_key, photret_file_abs)
				return

		xys = self.__get_xys_on_image(image_key)
		image = self.__get_internal_image(image_key, which_dir=which_dir)
		self.get_apphot_iraf_parameters(image_key, saveparfile=True)
		options = self.apphot_iraf_options.copy()

		photret = self.__aperture_photometry_apphot_iraf_single_image(image,xys,options, centering=centering)
		#xpos, ypos, flux, bkg,bkg_stdev, sky_counts, target_area,mag,mag_err

		photret_file_big = image_key.split('.')[0] +'_whole' + self.apphot_ret_file_suffix
		photret_file_big_abs = os.path.join(self.aperture_photometry_dir,photret_file_big)
		self.__delete_file_if_exist(photret_file_big_abs)
		np.savetxt(photret_file_big_abs,photret,fmt="%6.2f %6.2f %10.3f %8.3f %8.3f %8.3f %8.3f %6.3f %6.3f")

		photret_small = photret[:,[0,1,7,8]]
		self.__delete_file_if_exist(photret_file_abs)
		np.savetxt(photret_file_abs,photret_small,fmt="%6.2f %6.2f %6.3f %6.3f")

                phot_source_reg =  os.path.join(self.aperture_photometry_dir, image_key.split('.')[0] + '.reg')
		create_ds9_region_file(photret_small, phot_source_reg, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)




	def __aperture_photometry_apphot_iraf(self,flt, which_dir = 'raw_image', aperture_size_fixed = False, centering = True):
		'''
		pipeline script for aperture photometry on images with filter 'flt'
		For each image: get objects xys; get apphot input parameters; perform aperture photometry

		Involved functions:
			self.__get_xys_on_image
			self.get_apphot_iraf_parameters
			self.__aperture_photometry_apphot_iraf_single_image
		'''

		info_flt = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)
		images = info_flt['name']

		for image_key in images:
			print image_key
			if self.photometry_info[image_key]['drop'] > 0:
				continue

			photret_file = image_key.split('.')[0]  + self.apphot_ret_file_suffix
			photret_file_abs = os.path.join(self.aperture_photometry_dir,photret_file)

			if os.path.exists(photret_file_abs):
				if not self.renew_aperture_photometry_retfile:
					continue

			xys = self.__get_xys_on_image(image_key)
			image = self.__get_internal_image(image_key, which_dir=which_dir)
			self.get_apphot_iraf_parameters(image_key, saveparfile=True)
			options = self.apphot_iraf_options.copy()

			photret = self.__aperture_photometry_apphot_iraf_single_image(image,xys,options, centering=centering)
			#xpos, ypos, flux, bkg,bkg_stdev, sky_counts, target_area,mag,mag_err

			photret_file_big = image_key.split('.')[0] +'_whole' + self.apphot_ret_file_suffix
			photret_file_big_abs = os.path.join(self.aperture_photometry_dir,photret_file_big)
			self.__delete_file_if_exist(photret_file_big_abs)
			np.savetxt(photret_file_big_abs,photret,fmt="%6.2f %6.2f %10.3f %8.3f %8.3f %8.3f %8.3f %6.3f %6.3f")

			photret_small = photret[:,[0,1,7,8]]

			self.__delete_file_if_exist(photret_file_abs)
			np.savetxt(photret_file_abs,photret_small,fmt="%6.2f %6.2f %6.3f %6.3f")

                	phot_source_reg =  os.path.join(self.aperture_photometry_dir, image_key.split('.')[0] + '.reg')
			create_ds9_region_file(photret_small, phot_source_reg, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)


	def __aperture_photometry_apphot_iraf_single_point(self,image,xy,options=None, imagefwhm=None, centering = True, aperture_size_fixed=1):
		'''
		The same as self.__aperture_photometry_apphot_iraf_single_image
		but no region file will be saved

		Please see self.__aperture_photometry_apphot_iraf_single_image for details

		INPUTS:
			image: the input image not image key
			xy:
			options: iraf appot parameter dictory: app, def_zeropoint, fwhmpsf, sky_sigma_down, sky_sigma_up, skyin, skywidth
			imagefwhm: the fwhm of the image, if provided this will overwrite the one in options
			centering: if True, then centering algorithm in self.apphot_iraf_calgorithm will be used
		'''

		coo_list =  'txtfiles/temp.coo'
		ap_list = 'txtfiles/temp.ap'

		datamin = self.apphot_iraf_datamin
		datamax = self.apphot_iraf_datamax

		if centering:
			calgorithm = self.apphot_iraf_calgorithm
		else:
			calgorithm = 'none'

		epadu = float(self.base.telescopes_info[self.current_telescope]['gain'])
		readnoise = float(self.base.telescopes_info[self.current_telescope]['readnoise'])

		if self.current_ccd_readnoise is not None:
			readnoise = self.current_ccd_readnoise
		if self.current_ccd_epadu is not None:
			epadu = self.current_ccd_epadu
		if options is None:
			options = self.apphot_iraf_options.copy()

		if imagefwhm is not None:
			options['fwhmpsf'] = imagefwhm
		else:
			print "Attention!!!: default fwhm value of %s will be used"%options['fwhmpsf']

		if not aperture_size_fixed:
			fwhm_img = options['fwhmpsf']
			if fwhm_img <0 or fwhm_img>50:
				raise ValueError("The star profile with negative fwhm or fwhm value larger than 50 pixels. Really??")
			options['app'] = 2*fwhm_img
			options['skyin'] = 3*fwhm_img
			options['skywidth']= 3*fwhm_img

		ret  = apphot_iraf(image, xy, calgorithm, coo_list, ap_list,options, datamin=datamin, datamax = datamax, epadu=epadu, readnoise=readnoise)
		print "xpos, ypos, flux, bkg,bkg_stdev, sky_counts, target_area,mag,mag_err"

		return ret



	def __aperture_photometry_apphot_iraf_single_image(self,image,xys,options, centering = True,source_reg_filename=None):
		'''
		underlying aperture photometry task which utilize 'apphot_iraf'.

		ds9 region file will be saved in 'source_reg_filename'

		Inputs:
			image:  absolute name of the image in interest
			xys:	the objects position on the images at which the aperture
				photometry will take place
			options: dict containing input parameters for aperture photometry
				{'app': xx,
				 'def_zeropt': xx,
				 'fwhmpsf': xx,
				 'sky_sigma_down': xx,
				 'sky_sigma_up': xx,
				 'skyin': xx,
				 'skywidth': xx}

			centroid_option: required by apphot_iraf function; default 'centroid'

		See photometry_collection.apphot_iraf for details of aperture photometry task
		'''

		coo_list =  'temp.coo'
		imgkey = image.split(',')[0]
		apphot_ret_filename = '%s_%s'%(imgkey, self.phot_iraf_apphot_suffix)
		ap_list = os.path.join(self.aperture_photometry_dir, apphot_ret_filename)

		image = image+'[0]'
		datamin = self.apphot_iraf_datamin
		datamax = self.apphot_iraf_datamax

		if centering:
			calgorithm = self.apphot_iraf_calgorithm
		else:
			calgorithm = 'none'

		epadu = float(self.base.telescopes_info[self.current_telescope]['gain'])
		readnoise = float(self.base.telescopes_info[self.current_telescope]['readnoise'])

		if self.current_ccd_readnoise is not None:
			readnoise = self.current_ccd_readnoise
		if self.current_ccd_epadu is not None:
			epadu = self.current_ccd_epadu

		ret  = apphot_iraf(image, xys, calgorithm, coo_list, ap_list,options, datamin = datamin, datamax = datamax, epadu=epadu, readnoise=readnoise)
		#xpos, ypos, flux, bkg,bkg_stdev, sky_counts, target_area,mag,mag_err
		print ret

		mags = ret[:,7]
		success = np.where(mags != 99.999)[0]
		ret_good =ret[success,:]

		if source_reg_filename is None:
			source_reg_filename =  os.path.join(self.ds9_region_dir, image.split('/')[-1].split('.')[0] + '.reg')
		create_ds9_region_file(ret_good, source_reg_filename, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		return ret_good


	def __transform_xys(self,dxfit_values,dyfit_values,input_xys):
		'''
		allowd input_xys: (x,y) or (N,2) array
		dxfit_values = [a0,a1,a2]
		dyfit_values = [b0,b1,b2]
		(x,y) ?==> (x',y'): x'= a0+a1*x+a2*y; y'=b0+b1*x+b2*y
		'''
		a0,a1,a2 = dxfit_values
		b0,b1,b2 = dyfit_values

		if isinstance(input_xys,list):
			x = input_xys[0]
			y = input_xys[1]
			xout= a0+a1*x+a2*y
			yout= b0+b1*x+b2*y
			output_xys = [xout,yout]
		else:
			output_xys = None
			for xy in input_xys:
				x = xy[0]
				y = xy[1]
				xout= a0+a1*x+a2*y
				yout= b0+b1*x+b2*y
				output_xy = np.array([xout,yout]).reshape((1,2))
				if output_xys is None:
					output_xys = output_xy
				else:
					output_xys = np.append(output_xys,output_xy,axis=0)

		return output_xys



	def __extract_transformation_fitting(self,trans_file):
		'''
		Read in data from transformation coefficient file
		'''
		fid = open(trans_file,'r')
		data= fid.readlines()
		dxfit = data[-2].strip()
		dyfit = data[-1].strip()

		dxfit_strs = dxfit.split('=')[1].split(',')
		dyfit_strs = dyfit.split('=')[1].split(',')
		dxfit_values = [float(xv) for xv in dxfit_strs]
		dyfit_values = [float(yv) for yv in dyfit_strs]

		return dxfit_values,dyfit_values



	def __fitsh_grtrans(self,input_list,output_list,input_trans,):
		'''
		For details on grtrans, please refer to http://fitsh.szofi.net/task/grtrans

		grtransh [options] <input> [-o <output>]
		'''
		self.__delete_file_if_exist(output_list)
		grtrans = os.path.join(self.base.fitsh_dir, 'grtrans')
		command = "%s -i %s -o %s --input-transformation %s --col-xy 1,2"%(grtrans,input_list,output_list,input_trans)
		try:
			os.system(command)
		except:
			print "grtrans failure..."



	def __fitsh_grmatch(self,ref_list,input_list, match_output, trans_output, mode ='file'):
		'''
		For details on grmatch refer to http://fitsh.szofi.net/task/grmatch

		Basics for grmatch:
		grmatch [options] -r <reference> -i <input> [-o <output>]

		The program 'grmatch' matches lines read from two input files,
		namely from a reference and from an input file.
		All implemented algorithms are symmetric, in the manner that
		the result should be the same if these two files
		are swapped. The only case when the order of these files is
		important is when a geometrical transformation is
		also returned (see point matching below), in this case
		the swapping of the files results the inverse form of the
		original transformation. The lines (rows) can be matched using various criteria.


		If there is no negative sign before the column index, the data are sorted
		in descending(!) order, therefore the lines with the lines with the highest(!)
		values are selected for triangulation. If there is a negative sign before the index,
		the data are sorted in ascending order by these values, therefore the lines with
		the smallest(!) values are selected for triangulation.

		Inputs:
			ref_list: filename or data depending on 'mode'
			input_list:
			match_output: the table file containing match sources
			trans_output: the transformation coefficient file
			mode: 'file' or 'data'
		'''

		self.__delete_file_if_exist(match_output)
		self.__delete_file_if_exist(trans_output)

		if mode == 'file':
			ref_list_filename = ref_list
			input_list_filename = input_list
		elif mode == 'data':
			ref_list_filename = "grmatch_ref_list_temp.txt"
			input_list_filename = "grmatch_input_list_temp.txt"
			np.savetxt(ref_list_filename,ref_list)
			np.savetxt(input_list_filename, input_list)

		grmatch = os.path.join(self.base.fitsh_dir, 'grmatch')
		command = "%s -r %s -i %s -o %s --output-transformation %s"%(grmatch,ref_list_filename,input_list_filename,match_output,trans_output)
		if self.fitsh_grmatch_type == 'point':
			command = command + ' --match-points'
			for kw in self.fitsh_grmatch_pointmatch_pars.keys():
				if self.fitsh_grmatch_pointmatch_pars[kw] is not None:
					command = command + ' ' + kw + ' ' + str(self.fitsh_grmatch_pointmatch_pars[kw])
		elif self.fitsh_grmatch_type == 'coord':
			command = command + ' --match-coord'
			for kw in self.fitsh_grmatch_coordmatch_pars.keys():
				if self.fitsh_grmatch_coordmatch_pars[kw] is not None:
					command = command + ' ' + kw + ' ' + str(self.fitsh_grmatch_coordmatch_pars[kw])
		elif self.fitsh_grmatch_type == 'id':
			command = command + ' --match-id'
			for kw in self.fitsh_grmatch_idmatch_pars.keys():
				if self.fitsh_grmatch_idmatch_pars[kw] is not None:
					command = command + ' ' + kw + ' ' + str(self.fitsh_grmatch_idmatch_pars[kw])
		else:
			raise ValueError('invalid input for self.fitsh_grmatch_type')

		try:
			os.system(command)
		except:
			print "grmatch failure"

		if mode == 'data':
			os.remove(ref_list_filename)
			os.remove(input_list_filename)

		if os.path.exists(match_output):
			matched_ret = np.loadtxt(match_output)
		else:
			matched_ret = None

		return matched_ret



	def __get_trans_matrix_diapl_cross(self,list1,list2):

		print "on the way"


	def __get_rough_cross_shift_by_imaga_correlation(self,ref_img,target_img):
		'''
		cross cross.par instrument.par output_file reference_image target_image

		the function give the relative shift of target image against the reference image
		suppose sources on refence image (x,y); sources on the target image (x',y')
		in the ideal circumstance, x = x'+dx; y = y'+dy

		'''

		#prepare the images(cut the image to 2048*2048 for cross correlation)
		in_img_name = ref_img
		cutted_ref_name = 'cutted_ref_img.fits'
		self.__cut_image(in_img_name,cutted_ref_name,1,1,2048,2048)

		in_img_name = target_img
		cutted_target_name = 'cutted_target_img.fits'
		self.__cut_image(in_img_name,cutted_target_name,1,1,2048,2048)

		shifts_file = "temp_shift.dat"
		self.__delete_file_if_exist(shifts_file)

		#cross correlate two images to get the shift amout;
		cross_par = os.path.join(self.parafile_dir,'cross.par')
		instrument_par = os.path.join(self.parafile_dir,'instrument.par')

		cross = os.path.join(self.base.diapl_dir, 'cross')
		cross_command = "%s %s %s %s %s %s >/dev/null" %(cross, cross_par,instrument_par,shifts_file, cutted_ref_name, cutted_target_name)
		if os.system(cross_command):
		    print "Error: sorry dude, cannot run", command
		    sys.exit(1)


		os.remove(cutted_ref_name)
		os.remove(cutted_target_name)

		fid = open(shifts_file)
		line = fid.read()
		prepared_target_image, dx, dy = tuple(line.split())

		os.remove(shifts_file)

		return dx,dy


	def image_cutout(self, x0, y0, xe, ye, which_dir='raw_image'):
		'''
		cut out the wanted portion
		'''
		for imgkey in self.images.keys():
			self.image_cutout_single(imgkey, x0,y0,xe,ye, which_dir=which_dir)


	def image_cutout_single(self, imgkey, x0, y0, xe, ye, which_dir='raw_image'):
		'''
		image cut-out
		'''
		input_image = self.__get_internal_image(imgkey, which_dir=which_dir)
		output_image = os.path.join(self.modified_image_dir, imgkey)
		temp_image = os.path.join(self.modified_image_dir, 'cuttemp.fits')
		if os.path.exists(temp_image):
			os.remove(temp_image)
		self.__cut_image(input_image, temp_image, x0, y0, xe, ye)
		shutil.move(temp_image, output_image)



	def __cut_image(self,in_img_name,cutted_img_name,x0,y0,xe,ye):
		'''
		cut the rectangle region defined by low left corner(x0,y0) and upper right corner (xe,ye) of the 'in_img_name'
		and save as 'out_img_name'
		'''
		cut = os.path.join(self.base.diapl_dir, 'cutfitsim')
		cut_command = "%s %s %s %s %s %s %s" % (cut,in_img_name,cutted_img_name,x0,xe,y0,ye)

		try:
			os.system(cut_command)
		except:
			print "cut image failure on %s"%in_img_name



	def __get_xys_on_image(self,image_key):
		'''
		read the star file and get the xys of found stars
		return (N,2) array
		'''
		stars_file = image_key.split('.')[0]+self.starfile_suffix

		stars_file_abs = os.path.join(self.stars_dir,stars_file)

		if os.path.isfile(stars_file_abs):
			data = np.loadtxt(stars_file_abs)
			xys = data[:,0:2]
		else:
			print "The star file %s does not exist... you can perform sources finding through 'sfind'"%stars_file_abs
			xys = None

		return xys



	def get_target_xys(self,tpl_method='astrometry', which_dir = 'raw_image', updatetable=1):
		'''
		get the target position on image
		INPUTS:
			tpl_method: 'astrometry', 'ds9' or 'manual'
			which_dir: the directory containing the image to work on
		'''
		if updatetable:
			self.__dict2table()
		self.__find_template_imagekey(updatetable=updatetable)

		names = self.result_table_newsaved['name']
#		flts = self.result_table_newsaved['flt']
#		flts_unique = np.unique(flts)

		for flt in self.templates.keys():
			self.__get_templates_xy(flt,method = tpl_method, which_dir = which_dir)

		for image_key in names:
			if self.photometry_info[image_key]['drop']>0:
				continue

			if self.photometry_info[image_key]['template'] == 1:
				continue

			if self.photometry_info[image_key]['x'] == 0.0 or self.photometry_info[image_key]['y'] == 0:

				self.__get_target_xy_single(image_key, which_dir=which_dir)

				if self.photometry_info[image_key]['x'] == 0.0 or self.photometry_info[image_key]['y'] == 0:
					self.photometry_info[image_key]['drop'] = 5



	def __get_templates_xy(self,flt,method = 'astrometry', which_dir='raw_image', updatetable=1):
		'''
		Input:
			method: 'astrometry' or 'ds9' or 'manual'
		'''

		print "working on %s band template"%flt
		self.__find_template_imagekey(updatetable=updatetable)
		tpl_imgkey = self.templates[flt]
		tpl_image_cal = self.templates_after_astrometry[flt]

		if self.photometry_info[tpl_imgkey]['x'] == 0.0 or self.photometry_info[tpl_imgkey]['y'] == 0 or self.renew_template_xys:
			if method == 'astrometry':
				xy = self.__get_tpl_xy_on_image_from_astrometry(flt, which_dir=which_dir)
			elif method == 'ds9':
				xy = self.__get_xy_on_image_from_ds9(image_key = tpl_imgkey, which_dir = which_dir)
			elif method == 'manual':
				x = raw_input("Enter the x position of the current supernova on %s band template image:"%flt)
				y = raw_input("Enter the y position of the current supernova on %s band template image:"%flt)
				xy = [float(x),float(y)]
			else:
				print "Invalid input for method"
				raise IOError("Sorry...")

			x_image = round(xy[0],1)
			y_image = round(xy[1],1)
			self.photometry_info[tpl_imgkey]['x'] = x_image
			self.photometry_info[tpl_imgkey]['y'] = y_image




	def __get_target_xy_single(self,image_key,depend_on_tpl = True, non_tpl_method = 'ds9', which_dir = 'raw_image'):
		'''
		Get the supernova position on images
		Input:

			depend_on_tpl: if True, the supernova position on the template image required and working-on image
					get its position relative the template image
				       if False, the working on image gets its position independently
			non_tpl_method: if depend_on_tpl is False, then individual image get sn position with method assigned by this input
		'''
		print image_key
		if depend_on_tpl:
			flt = self.photometry_info[image_key]['flt']
			tpl_imgkey = self.templates[flt]
			ref_list   = os.path.join(self.stars_dir,tpl_imgkey.split('.')[0]+self.starfile_suffix)
			input_list = os.path.join(self.stars_dir,image_key.split('.')[0]+self.starfile_suffix)
			match_output = os.path.join(self.stars_dir,image_key.split('.')[0]+ '_tpl.match')
			self.__delete_file_if_exist(match_output)
			trans_fitting = os.path.join(self.stars_dir,image_key.split('.')[0] + '_tpl.coef')
			self.__delete_file_if_exist(trans_fitting)
			matched_table = self.__fitsh_grmatch(ref_list, input_list, match_output, trans_fitting)
			if matched_table is None:
				self.photometry_info[image_key]['drop'] = 5
				x_image = 0
				y_image = 0
			else:
				trans_output = os.path.join(self.stars_dir,image_key.split('.')[0]+'_tpl.trans')
				self.__delete_file_if_exist(trans_output)
				self.__fitsh_grtrans(match_output,trans_output,trans_fitting)
				dxfit,dyfit = self.__extract_transformation_fitting(trans_fitting)
				tpl_xy = [self.photometry_info[tpl_imgkey]['x'],self.photometry_info[tpl_imgkey]['y']]
				target_xy = self.__transform_xys(dxfit,dyfit,tpl_xy)
				print target_xy
				x_image = round(target_xy[0],1)
				y_image = round(target_xy[1],1)

				if x_image < 0 or y_image <0:
					self.photometry_info[image_key]['drop'] = 5
		else:
			if non_tpl_method == 'ds9':
				xy = self.__get_xy_on_image_from_ds9(image_key = image_key, which_dir=which_dir)
			elif non_tpl_method == 'manual':
				x = raw_input("Enter the x position of the current supernova on image %s:"%image_key)
				y = raw_input("Enter the y position of the current supernova on image %s:"%image_key)
				xy = [float(x),float(y)]
			else:
				print "Invalid input for method"
				raise IOError("Sorry...")

			x_image = round(xy[0],1)
			y_image = round(xy[1],1)

		self.photometry_info[image_key]['x'] = x_image
		self.photometry_info[image_key]['y'] = y_image




	def __get_xy_on_image_from_wcsinfo(self, inputimage, verbose=0):
		'''
		world2image conversion: wcs info in fits header required
		'''
		if self.sn_ra_world_deg is None or self.sn_dec_world_deg is None:
			if self.current_sn in self.base.registerd_photometry_objects_info.keys():
				self.sn_ra_world_deg  = float(self.base.registerd_photometry_objects_info[self.current_sn]['ra'])
				self.sn_dec_world_deg = float(self.base.registerd_photometry_objects_info[self.current_sn]['dec'])

		if self.sn_ra_world_deg is None or self.sn_dec_world_deg is None:
			sn_ra = raw_input("Enter the RA of the current supernova:")
			sn_dec = raw_input("Enter the Dec of the current supernova:")
			self.sn_ra_world_deg = float(sn_ra)
			self.sn_dec_world_deg = float(sn_dec)

		world_radec = np.array([self.sn_ra_world_deg,self.sn_dec_world_deg])
		world_radec = world_radec.reshape(1,2)
		image_xy = self.__world2image(inputimage,world_radec)
		if verbose>0:
			print image_xy

		if verbose>1:
			image = self.images[self.templates[flt]]
			self.label_single_xy_region(image,image_xy)


		return image_xy[0]


	def __get_tpl_xy_on_image_from_astrometry(self, flt, verbose=1):
		'''
		Output:
			image_xy: array
		'''
		self.__find_template_imagekey()
		inputimage = self.templates_after_astrometry[flt]
		tplxy = self.__get_xy_on_image_from_wcsinfo(inputimage)
		return tplxy



	def __world2image_fitshead(self, hdr, worldcoor):
		w = WCS(hdr)
		image_coor = w.wcs_world2pix(worldcoor,1) # Pixel coordinates of (RA, DEC)
		return image_coor

	def __world2image(self,image,world, hduindex=None):
		hdu  = fits.open(image)
		if hduindex is None:
			header = hdu['PRIMARY'].header
		else:
			header = hdu[hduindex].header

		image_coor = self.__world2image_fitshead(header, world)

		return image_coor

	def __image2world_fitshead(self, hdr, imagecoor):
		'''
		convert from image coordinate to physical world coordinate

		'''
		w = WCS(hdr)
		world_coor = w.wcs_pix2world(imagecoor,1)

		return world_coor

	def __image2world(self,image,image_xys, hduindex=None):
		hdu = fits.open(image)

		if hduindex is None:
			header = hdu['PRIMARY'].header
		else:
			header = hdu[hduindex].header

		world_coor = self.__image2world_fitshead(header, image_xys)

		return world_coor



	def __get_xy_on_image_from_ds9(self,image_key=None,input_image=None, which_dir='raw_image', regionfile =None, newds9=True):
		'''
		pick up xy coordinate from ds9 image
		'''
		if input_image is None:
			if image_key is not None:
				input_image = self.__get_internal_image(image_key, which_dir=which_dir)
			else:
				raise IOError("No valid input...")
		else:
			if os.path.exists(input_image):
				print input_image
			else:
				raise IOError("input image doesn't exist!")

		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		if regionfile is not None:
			print regionfile
                	d.set('regions load %s'%regionfile)

		xys = d.get('imexam coordinate image')
		x,y = xys.split()
		xy = np.array([float(x),float(y)])
		print xy

		return xy


	def get_standards(self,catalog = 'apass', center_ra = None, center_dec = None, distance = None):
		'''
		search the standard star catalog to extract subset of the data for future match with stars in working image field
		currently supported catalogs: apass and 2mass

		INPUTS:
			catalog: 'apass', '2mass' or 'panstarrs'
			center_ra: catalog search center ra; if None, self.sn_ra_world_deg will be used
			center_dec: catalog search center dec; if None, self.sn_dec_world_deg will be used
			distance: search radius in degree; if None, self.std_region_radius will be used
		'''

		std_ref_stars_file = os.path.join(self.std_ref_dir, self.std_ref_stars_file[catalog])

		renew_stds = self.renew_standards
		if os.path.exists(std_ref_stars_file) and renew_stds:
			os.remove(std_ref_stars_file)


		if not os.path.exists(std_ref_stars_file):
			if center_ra is None or center_dec is None:
				center_ra, center_dec = self.__get_standards_prepare()
			if catalog == 'panstarrs' and center_dec<-30:
				raise ValueError('PanSTARRS catalog does not cover region beyond -30 degree for the south')
			if distance is None:
				distance = self.std_region_radius
			self.__get_standards(catalog=catalog, center_ra = center_ra, center_dec =center_dec, distance = distance)
		else:
			if catalog == 'apass':
				self.standards = Table.read(std_ref_stars_file, format='ascii')
			elif catalog == '2mass':
				self.standards = Table.read(std_ref_stars_file, format='ascii.csv')
			elif catalog == 'panstarrs':
				self.standards = Table.read(std_ref_stars_file, format='ascii.csv')
			else:
				print "Warning: the catalog %s not supported yet"%catalog


	def __get_standards(self, catalog='apass', center_ra=None, center_dec=None, distance=None):
		'''
		standard reference data will be saved in self.std_ref_stars_file[catalog]
		and self.standards will be updated

		INPUTS:
			Refer to self.get_standards

		'''

		if center_ra is None or center_dec is None:
			raise ValueError('subset region center (center_ra, center_dec) required...')

		if distance is None:
			raise ValueError('standard region radius required...')

		std_ref_dir = self.std_ref_dir
		std_ref_stars_file = self.std_ref_stars_file[catalog]
		output_filename = os.path.join(std_ref_dir, std_ref_stars_file)

		if catalog == '2mass':
			twomassmethod = self.twomass_catalog_method
			if twomassmethod == 1:
				query_VO_SCS(center_ra, center_dec,distance,table='fp_psc',out_format='csv', output_file = output_filename)
			elif twomassmethod == 2:
				twomass_query_Vizier(center_ra, center_dec, distance, outfile=output_filename)
			else:
				raise IOError('2MASS catalog download method %s not available yet'%str(twomassmethod))
			self.standards = Table.read(output_filename, format='ascii.csv')
		elif catalog == 'apass':
			apassmethod  = self.apass_catalog_method
			if apassmethod == 1:
				stdref_database_dir  = self.local_apass_dir
				if not os.path.exists(stdref_database_dir):
					raise ValueError("the local apass data expected in %s not available"%stdref_database_dir)
				query_local_APASS(center_ra, center_dec,distance,output_file=output_filename, local_database_dir=stdref_database_dir)
				self.standards = Table.read(output_filename, format='ascii')
			elif apassmethod == 2:
				apass_query_Vizier(center_ra, center_dec, distance, outfile=output_filename)
				self.standards = Table.read(output_filename, format='csv')
			else:
				raise IOError('PS1 catalog download method %s not available yet'%str(PS1method))
		elif catalog == 'panstarrs':
			PS1method = self.panstarrs_catalog_method
			print "The standard stars will be saved at %s"%output_filename
			if PS1method ==1:
				if self.std_region_radius > 0.25:
					raise ValueError('sorry, radius too large...')
				query_General_MAST(center_ra, center_dec, distance, FORMAT='csv', catalog='PS1V3OBJECTS', filename=output_filename)
			elif PS1method == 2:
				panstarrs_query_Vizier(center_ra, center_dec, distance, outfile=output_filename)
			else:
				raise IOError('PS1 catalog download method %s not available yet'%str(PS1method))

			self.standards = Table.read(output_filename, format='ascii.csv')

		else:
			print "Warning: the catalog %s not supported yet"%catalog


	def __get_standards_prepare(self, verbose=1):
		'''
		get the coordinate of the target in interest
		'''
		if self.sn_ra_world_deg is not None and self.sn_dec_world_deg is not None:
			if verbose:
				print "target (RA,Dec): (%s, %s)"%(self.sn_ra_world_deg, self.sn_dec_world_deg)
		else:
			sn_ra = raw_input("Enter the RA of the current supernova:")
			sn_dec = raw_input("Enter the Dec of the current supernova:")
			self.sn_ra_world_deg = float(sn_ra)
			self.sn_dec_world_deg = float(sn_dec)

		ra = self.sn_ra_world_deg
		dec = self.sn_dec_world_deg

		return ra,dec


	def load_2mass_standards_to_ds9_image(self, input_image,text1col=None, text2col=None,newds9=True):
		'''
		load 2mass catalog in DS9 to given image
		The standard stars are loaded from catalog file standard_reference_star_2mass.txt

		INPUTS:
		'''
		twomassstdfile = os.path.join(self.std_ref_dir, self.std_ref_stars_file['2mass'])
		if not os.path.exists(twomassstdfile):
			raise ValueError('standard star catalog %s not exists... you can prepare that with self.get_standards()'%twomassstdfile)

		standards = Table.read(twomassstdfile, format='ascii.csv')
		self.__load_standards_region_to_ds9_image(input_image, standards, color='red', text1col=text1col, text2col=text2col, newds9=newds9)



	def load_panstarrs_standards_to_ds9_image(self, input_image, text1col=None, text2col=None, magcol=None, magbcut=None, magfcut=None, nbright=None, fbright=None, newds9=True):
		'''
		load PS1 catalog in DS9 to given image
		The standard stars are loaded from catalog file standard_reference_star_panstarrs.txt
		INPUTS:
			input_image:
			magcol: the magnitude column providing data for filtering
			magbcut: the bright end cut of the loaded sources; bright corresponds to small value of magnitudes
			magfcut: the faint end cut of the loaded sources
			nbright: the brightest number of targets to be loaded;
			fbright: the fraction of the brightest sources to be loaded
		'''

		ps1stdfile = os.path.join(self.std_ref_dir, self.std_ref_stars_file['panstarrs'])
		if not os.path.exists(ps1stdfile):
			raise ValueError('standard star catalog %s not exists... you can prepare that with self.get_standards()'%ps1stdfile)

		standards = Table.read(ps1stdfile, format='ascii.csv')
		if standards.masked:
			standards = standards.filled(99.99)

		if np.any(np.array([magbcut, magfcut, nbright, fbright])):
			if magcol is None:
				raise ValueError('the magnitude column required to provide data for filtering ')
			if magbcut is not None and magfcut is not None:
				if magbcut > magfcut:
					raise ValueError("bright end cut must brighter than the faint end cut")
				standards = standards[standards[magcol]>magbcut*standards[magcol]<magfcut]
			elif magbcut is not None:
				standards = standards[standards[magcol]>magbcut]
			elif magfcut is not None:
				standards = standards[standards[magcol]<magfcut]
			else:
				standards.sort(magcol)
				Nn = nbright
				if nbright is not None and fbright is not None:
					Nf = len(standards)*fbright
					if Nn > Nf:
						standards = standards[:Nf]
					else:
						standards = standards[:Nn]
				else:
					if nbright is not None:
						standards = standards[:Nn]
					if fbright is not None:
						Nf = len(standards)*fbright
						standards = standards[:Nf]

		self.__load_standards_region_to_ds9_image(input_image, standards, color='red', text1col=text1col, text2col=text2col, newds9=newds9)


	def load_apass_standards_to_ds9_image(self, img, text1col=None, text2col=None,newds9=True):
		'''
		The standard stars are loaded from catalog file standard_reference_star_apass.txt

		See self.__load_standards_region_to_ds9_image
		'''
		apassstdfile = os.path.join(self.std_ref_dir, self.std_ref_stars_file['apass'])
		if not os.path.exists(apassstdfile):
			raise ValueError('standard star catalog %s not exists... you can prepare that with self.get_standards()'%apassstdfile)

		standards = Table.read(apassstdfile, format='ascii.csv')
		self.__load_standards_region_to_ds9_image(img, standards, color='red', text1col=text1col, text2col=text2col, newds9=newds9)

	def stat_apass_standards(self, outfile=None):
		'''
		display the mag VS. e_mag plot
		'''
		if self.standards is None:
			raise ValueError('prepare the standard catalog first')

		fig = plt.figure(figsize=(8,7))
		ax1 = fig.add_subplot(221)
		ax2 = fig.add_subplot(222)
		ax3 = fig.add_subplot(223)
		ax4 = fig.add_subplot(224)

		for ax,magcol in zip([ax1,ax2,ax3,ax4],['Bmag','Vmag','r_mag','i_mag']):
			ax.plot(self.standards[magcol], self.standards['e_'+magcol],'.')
			ax.set_xlabel(magcol)
			ax.set_ylabel('e_'+magcol)
			ax.set_ylim([-0.02, 0.3])

		if outfile is not None:
			savefile = os.path.join(self.std_ref_dir, outfile)
			plt.savefig(savefile)

		plt.show()

	def load_current_standards_to_ds9_image(self, img, text1col=None, text2col=None, newds9=True):
		'''
		The standard stars are from self.standards
		INPUTS:
			img: the input image (with absolute path)
			text1col: the column name in the standard star table
			text2col: the column name in the standard star table
		'''
		if not self.standards:
			raise ValueError("please load standard catalog to self.standards first")

		self.__load_standards_region_to_ds9_image(img, self.standards, color='red', text1col=text1col, text2col=text2col, newds9=newds9)



	def __load_standards_region_to_ds9_image(self, img, standards, color='red', text1col=None, text2col=None, newds9=True):
		'''
		load the input standards to ds9 image; WCS info required in the image header

		INPUTS:
			img:
			standards: the input standards catalog
			newds9: whether start a new pyds9 new instance
		'''
		if standards.masked:
			standards = standards.filled(99.99)

		region_filename = os.path.join(self.std_ref_dir, 'standards_temp.reg')
		stdcolnames = standards.colnames
		if ('RAJ2000' in stdcolnames) and ('DEJ2000' in stdcolnames):
			RA_colname  = 'RAJ2000'
			Dec_colname = 'DEJ2000'
		elif ('RA' in stdcolnames) and ('Dec' in stdcolnames):
			RA_colname  = 'RA'
			Dec_colname = 'Dec'
		elif  ('ra' in stdcolnames) and ('dec' in stdcolnames):
			RA_colname  = 'ra'
			Dec_colname = 'dec'
		else:
			raise ValueError("catalog columns names are not supported yet...")

		try:
                        RAs  = standards[RA_colname].data
                        Decs = standards[Dec_colname].data
		except:
			print "please check whether input parameter std_catalog is correct..."
		if text1col is not None and text1col not in stdcolnames:
			raise ValueError("the required colume %s not exist"%text1col)
		if text2col is not None and text2col not in stdcolnames:
			raise ValueError("the required colume %s not exist"%text2col)

		if text1col is not None and text2col is not None:
			loadtext = True
			textcol1 = 2
			textcol2 = 3
			input_sources = np.array([RAs, Decs, standards[text1col].data, standards[text2col].data]).transpose()
		elif text1col is not None:
			loadtext = True
			textcol1 = 2
			textcol2 = None
			input_sources = np.array([RAs, Decs, standards[text1col].data]).transpose()
		elif text2col is not None:
			loadtext = True
			textcol1 = None
			textcol2 = 2
			input_sources = np.array([RAs, Decs, standards[text2col].data]).transpose()
		else:
			loadtext = False
			textcol1 = None
			textcol2 = None
			input_sources = np.array([RAs, Decs]).transpose()

		create_ds9_region_file(input_sources, region_filename, x_col=0, y_col=1, coordinate_type = 'fk5', radec_deg = True, circle_radius = 0.0042, color=color, load_text=loadtext, text_uoffset=0.005, text_loffset=0.005, textcol1=textcol1, textcol2=textcol2) #circle radius around 15 arcsec = APASS aperture size for photometry

		if not os.path.exists(img):
			raise IOError("image %s doesn't exist"%img)

		d = self.__display_image_with_ds9(img, newds9=newds9)
		self.__load_source_regions_to_ds9(d, region_filename)




	def show_reference_stars_with_id_on_template_image_flt(self, flt, which_dir='raw_image'):
		'''
		This is used to display the reference stars with labeled ID which aids selecting corresponding stars in input images or from standards stars
		'''
		if len(self.templates) == 0:
			self.__find_template_imagekey()

		if flt not in self.templates.keys():
			raise ValueError("no template image avaiable for filter %s"%flt)

		tpl_imgkey = self.templates[flt]
		tpl_imgkey_s = tpl_imgkey.split('.')[0]
		refmag_filename = os.path.join(self.template_dir, "%s_%s_idxymag.txt"%(flt, tpl_imgkey_s))
		refmag_table = Table.read(refmag_filename, format='ascii.fixed_width')

		xymagids = np.transpose(np.array([refmag_table['x'], refmag_table['y'], refmag_table['mag'], refmag_table['num']]))
		input_image = self.__get_internal_image(tpl_imgkey, which_dir=which_dir)
		regionfile = os.path.join(self.template_dir, "%s_%s_idxymag.reg"%(flt, tpl_imgkey_s))

		self.__load_stars_xy_mag_to_ds9(input_image, xymagids, regionfile, mag_offset = 0, coortype = 'image', color = 'red', width=1, radius=10)


	def show_standards_for_wcs_not_available(self,flt, std_catalog='apass', radius_arcmin = 16, rotate_deg = 0, brightest_N = None, invert_x=False, invert_y = False, x_mid = None, y_mid = None, display_mags = False, which_dir='raw_image'):
		'''
		prepare and display standard stars on the template image for the purpose of 'mannually' matching standard stars with calibration stars on the template image

		INPUTS:
			flt:
			std_catalog: 'apass', '2mass'
			radius_arcmin:
			rotate_deg:
			brightest_N:
			invert_x:
			invert_y:
			x_mid:
			y_mid:
			display_mags:
		'''

		if std_catalog == 'apass':
			corresmap = {'RA':"RA", 'Dec':"Dec", 'V':"V", 'B':"B",'gp':"gp", 'rp':"rp", 'ip':"ip", 'Verr':"Verr", 'Berr':"Berr", 'gperr':"gperr", 'rperr':"rperr", 'iperr':"iperr" }
		elif std_catalog == '2mass':
			corresmap = {'RA':"ra", 'Dec':"dec", 'J':"j_m", 'Jerr':"j_msig", 'H':"h_m", 'Herr':"h_msig",'K':"k_m", 'Kerr':"k_msig", }
		else:
			raise ValueError("%s not supported..."%std_catalog)


		racen  = self.sn_ra_world_deg
		deccen = self.sn_dec_world_deg

		if self.standards is None:
			self.get_standards(catalog=std_catalog)

		ras  = self.standards[corresmap['RA']]
		decs = self.standards[corresmap['Dec']]

		#
		if flt == 'I':
			flt = 'ip'
			print "SDSS-i band data used for I band"
		if flt == 'R':
			flt = 'rp'
			print "SDSS-r band data used for R band"

		#mag_col = corresmap[flt]
		mags    = self.standards[corresmap[flt]]
		magerrs = self.standards[corresmap[flt+'err']]

		distances = np.sqrt((ras-racen)**2 + (decs-deccen)**2)

		radius_deg = radius_arcmin/60.

		ra_yes    = ras[distances<radius_deg]
		dec_yes  = decs[distances<radius_deg]
		mags_yes = mags[distances<radius_deg]
		magerrs_yes = magerrs[distances<radius_deg]

		self.__find_template_imagekey()
		tpl_imgkey = self.templates[flt]
		x0 = self.photometry_info[tpl_imgkey]['x']
		y0 = self.photometry_info[tpl_imgkey]['y']

		if x0 == 0 or y0==0:
			raise ValueError("the image coordinate for the target is required")

		if self.pixscale is None:
			raise ValueError("please specify the pixel scale")

		pixscale_deg = self.pixscale / 3600.

		x = (ra_yes - racen)/pixscale_deg
		y = (dec_yes - deccen)/pixscale_deg

		rotate_rad = rotate_deg/180.*np.pi

		xp =  x*np.cos(rotate_rad) + y*np.sin(rotate_rad)
		yp = -x*np.sin(rotate_rad) + y*np.cos(rotate_rad)

		X = x0-xp	#left east and right west
		Y = y0+yp

		if invert_x:
			if x_mid is None:
				x_mid = (np.max(X) + np.min(X))/2
			X = 2*x_mid - X
		if invert_y:
			if y_mid is None:
				y_mid = (np.max(Y) + np.min(Y))/2
			Y = 2*y_mid - Y

		if brightest_N is None:
			N = len(mags)
		else:
			N = brightest_N

		mags_sort = np.argsort(mags_yes)


		X_plot = X[mags_sort[:N]]
		Y_plot = Y[mags_sort[:N]]
		mags_plot = mags_yes[mags_sort[:N]]
		magerrs_plot = magerrs_yes[mags_sort[:N]]

		xymags = np.transpose(np.array([X_plot, Y_plot, mags_plot, magerrs_plot]))
		input_image = self.__get_internal_image(tpl_imgkey, which_dir=which_dir)

		tpl_skey = tpl_imgkey.split('.')[0]
		std_xymags_file = os.path.join(self.std_ref_dir, "%s_%s_std_xymag.txt"%(flt, tpl_skey))
		np.savetxt(std_xymags_file, xymags, fmt="%8.2f %8.2f %8.3f %8.3f")

		regionfile = os.path.join(self.std_ref_dir, "%s_%s_std_xymag.reg"%(flt, tpl_skey))
		self.__load_stars_xy_mag_to_ds9(input_image, xymags, regionfile, mag_offset = 0, coortype = 'image', color = 'red', width=1, radius=10)


		return xymags



	def __find_template_imagekey(self, updatetable=1):
		'''
		find template images
		'''
		if updatetable:
			self.__dict2table()

		info_notdrop = self.__select_rows_from_table(self.result_table_newsaved,'drop',0)

		info_template  = self.__select_rows_from_table(info_notdrop,'template',1)

		#print info_template

		for img_key in info_template['name']:
			flt = self.photometry_info[img_key]['flt']
			self.templates[flt] = img_key
			self.templates_after_astrometry[flt] = os.path.join(self.template_dir,'cal_' + img_key)


	def get_template_image(self):

		self.__get_template_image()


	def renew_template_image(self):

		for image_key in self.photometry_info.keys():
			self.photometry_info[image_key]['template'] = 0

		self.__get_template_image()


	def __get_template_image(self, updatetable=1):
		'''
		assign template image for each filter

		This function will not remove old template image
		refer to self.renew_template_image, if you want to assign new template image
		'''
		self.__dict2table()

		imgs_info = self.result_table_newsaved
		imgs_info = self.__select_rows_from_table(imgs_info,'drop',0)

		#print imgs_info

		names = imgs_info['name']
		flts = imgs_info['flt']
		flts_unique = np.unique(flts)

		print flts_unique

		self.__find_template_imagekey(updatetable=updatetable)
		flts_tpl_old = self.templates.keys()


		for flt in flts_unique:

			if flt in flts_tpl_old:
				continue

			flt_mask = flts==flt
			if np.sum(flt_mask)==1:
				image_key = names[flt_mask][0]
				self.photometry_info[image_key]['template'] = 1
			else:
				self.__get_template_image_single(flt)




	def __get_template_image_single(self,flt,verbose=1):
		'''
		if verbose==2 the template will be displayed and saved by aplpy
		'''
		self.__dict2table()
		info_sf = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt,)
		info_sf = self.__select_rows_from_table(info_sf,'drop',0)

		bkg     = info_sf['bkg']
		fwhm    = info_sf['fwhm']
		starnum = info_sf['nstar']

		bkg_mask = bkg != self.anomaly_signs['bkg']
		fwhm_mask = fwhm != self.anomaly_signs['fwhm']
		starnum_mask = starnum != self.anomaly_signs['nstar']
		valid_mask = bkg_mask*fwhm_mask*starnum_mask

		info_sf = info_sf[valid_mask]
		bkg     = info_sf['bkg']
		fwhm    = info_sf['fwhm']
		starnum = info_sf['nstar']
		coadd_template = 0 	# use the best single image as the template
		if coadd_template == 0:
		    num = 1		#how many best images you want for the given filter
		    inputpara = np.hstack((bkg.reshape(len(bkg),1),fwhm.reshape(len(fwhm),1),starnum.reshape(len(starnum),1)))
		    if verbose:
			print inputpara
		    ret,wantedindice=get_best(inputpara,self.template_factors,num)
		    image_key = info_sf['name'][wantedindice[0]]
		    self.photometry_info[image_key]['template'] = 1

		    if verbose == 2:
			import aplpy
			tpl_image_name = 'template_'+flt+'.eps'
			tpl_image_save = os.path.join(self.template_dir,tpl_image_name)
			self.__delete_file_if_exist(tpl_image_save)
			img = aplpy.FITSFigure(self.images[image_key])
			img.show_grayscale()
			img.save(tpl_image_save)



	def __delete_file_if_exist(self,filename):
		'''
		if file 'filename' exist, remove the existing file
		'''
		if os.path.isfile(filename):
			os.remove(filename)



	def assign_template_image(self,flt,image_key):
		'''
		delete old if exist and add new one assigned
		'''
		self.__dict2table()
		info_sf = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt,)
		names = info_sf['name']
		template_flags     = info_sf['template']

		for img in names:
			self.photometry_info[img]['template'] = 0
		if image_key not in names:
			raise ValueError("%s not in %s band"%(image_key,flt))
		self.photometry_info[image_key]['template'] = 1



	def __select_rows_from_table(self,table,colname,value,mode='eq',criteria = 0):
		'''
		select rows in table which meet specific requirements given by mode and criteria
		available mode 'eq','lt','gt','lgt'

		'eq':X == value; 'lt': X<= value; 'gt': X>=value; 'lgt': value-criteria  <= X <= value+criteria
		if mode is 'lgt' then non-zero criteria is needed, otherwise 'lgt' is the same as 'eq' with default criteria

		return the table by rows which meet the criteria
		'''
		colvalues = table[colname]
		if mode   == 'eq':
			mask  = colvalues == value
		elif mode == 'lt':
			mask  = colvalues <= value
		elif mode == 'gt':
			mask  = colvalues >= value
		elif mode == 'lgt':
			if criteria == 0:
				print "Attention: if mode is 'lgt' then non-zero criteria is needed, otherwise 'lgt' is the same as 'eq' with default criteria!!!"
			mask = np.logical_and(colvalues>=value-criteria,colvalues<=value+criteria)
		else:
			raise KeyError("the criteria %s is not available"%criteria)

		selected = table[mask]

		return selected


	def __select_columns_from_table_to_ascii(self, table, colnames, outfile):
		'''
		output data in colnames to pure ascii data file
		'''
		outdata = table[colnames]
		self.__delete_file_if_exist(outfile)
		outdata.write(outfile, format='ascii.fast_commented_header')


	def get_star_number_rough(self,method='sfind', which_dir = None, fistar_ntime_above_bkg = 3, daofind_nsigma_thresh = 3):
		'''
		get the approximate number of stars in the observation images
		INPUTS:
			method: 'sfind' or 'fistar' or 'daofind'
			which_dir:
			fistar_ntime_above_bkg:
		'''
		for key in self.photometry_info.keys():
			if which_dir is None:
				if self.photometry_info[key]['bitpix'] == -32:
					imgdir  = 'modified_image'
				else:
					imgdir = 'raw_image'
			else:
				imgdir = which_dir

			if self.photometry_info[key]['nstar'] == 0 and self.photometry_info[key]['drop'] == 0:
				self.__get_star_number_rough(key,method=method,ntime_above_bkg = fistar_ntime_above_bkg,  daofind_nsigma_thresh = daofind_nsigma_thresh, which_dir = imgdir)


	def renew_star_number_rough(self,method='sfind', which_dir = 'modified_image', fistar_ntime_above_bkg = 3,  daofind_nsigma_thresh = 3):
		'''
		get the approximate number of stars in the observation images
		'''
		for key in self.photometry_info.keys():
			if self.photometry_info[key]['drop'] == 0:
				self.__get_star_number_rough(key,method=method,ntime_above_bkg = fistar_ntime_above_bkg,  daofind_nsigma_thresh = daofind_nsigma_thresh, which_dir = which_dir)



	def __get_star_number_rough(self,key,method = 'sfind',ntime_above_bkg =3,  daofind_nsigma_thresh = 3, drop_nstar_less_than_this = 10, which_dir = 'modified_image'):

		if method == 'sfind':
			try:
				print key
				starnum = self.__get_star_number_rough_single_sfind_diapl(key, which_dir = which_dir)
				self.photometry_info[key]['nstar'] = starnum
				if starnum < drop_nstar_less_than_this:
					self.photometry_info[key]['drop'] = 2
			except Exception,e:
				nstar = self.anomaly_signs['nstar']
				print "source finding or region file saving error on %s"%key
				self.photometry_info[key]['nstar'] = nstar
				self.photometry_info[key]['drop'] = 2
		elif method == 'fistar':
			try:
				print key
				starnum = self.__get_star_number_rough_single_fistar_fitsh(key,ntime_above_bkg = ntime_above_bkg, which_dir = which_dir)
				self.photometry_info[key]['nstar'] = starnum
				if starnum < drop_nstar_less_than_this:
					self.photometry_info[key]['drop'] = 2
			except Exception,e:
				nstar = self.anomaly_signs['nstar']
				print "source finding or region file saving error on %s"%key
				self.photometry_info[key]['nstar'] = nstar
				self.photometry_info[key]['drop'] = 2
		elif method == 'daofind':
			try:
				print key
				starnum = self.__get_star_number_rough_single_apphot_daofind(key, which_dir = which_dir, threshold=daofind_nsigma_thresh )
				self.photometry_info[key]['nstar'] = starnum
				if starnum < drop_nstar_less_than_this:
					print "%s stars detected on %s, and drop value 2 will be assigned"%(starnum, key)
					self.photometry_info[key]['drop'] = 2
			except Exception,e:
				nstar = self.anomaly_signs['nstar']
				print "source finding or region file saving error on %s"%key
				self.photometry_info[key]['nstar'] = nstar
				self.photometry_info[key]['drop'] = 2

		else:
			raise IOError("Invalid input for method...")



	def __get_star_number_rough_single_apphot_daofind(self, image_key,which_dir='raw_image', fwhmpsf=None, emission=True, sigma=None, datamin=None, datamax=None, readnoise=None, epadu=None, threshold=3, deleteINDEF = True, nsigma=1.5,  psffwhm_hc=30):
		'''
		source detection using daofind from iraf apphot package

		'''
		imgkey= image_key.split('.')[0]
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		stars_filename = os.path.join(self.stars_dir,imgkey+self.starfile_suffix)
		self.__delete_file_if_exist(stars_filename)

		daofind_outfilename = '%s%s'%(imgkey, self.iraf_daofind_output_suffix)
		iraf_daofind_outfile = os.path.join(self.stars_dir, daofind_outfilename)

		if fwhmpsf is None:
			fwhmpsf= self.photometry_info[image_key]['fwhm']

		if fwhmpsf <= 0 or fwhmpsf >psffwhm_hc: #reasonability check
			raise ValueError("unrealistic fwhm of stellar profile, please check...")

		if epadu is None:
			try:
				epadu = float(self.base.telescopes_info[self.current_telescope]['gain'])
			except:
				print "gain value of 1.0 will be used for %s"%image_key
				epadu = 1.0

		if sigma is None:
			stat_retdict = imstat_iraf(input_image,fields = "image,npix,mode, midpt,mean,stddev,min,max",lower = 'INDEF',upper = 'INDEF',nclip = 5,lsigma = 3.0,usigma = 3.0,binwidth = 1.0, verbose=0)
			sigma = stat_retdict['stddev']
		if datamin is None:
			datamin = 'INDEF'
		if datamax is None:
			datamax = 'INDEF'
		if readnoise is None:
			try:
				readnoise = float(self.base.telescopes_info[self.current_telescope]['readnoise'])
			except:
				print "readnoise value of 10.0 will be used for %s"%image_key
				readnoise = 10.0

		daofind_iraf(input_image, output=stars_filename, fwhmpsf=fwhmpsf, emission=emission, sigma=sigma, datamin=datamin, datamax=datamax, readnoise=readnoise, epadu=epadu, threshold=threshold, nsigma=nsigma)
		if deleteINDEF:
			savelines_index = []
			source_lines = np.array(open(stars_filename).readlines())
			for i,line in enumerate(source_lines):
				if i<41:
					continue
				linesegs= line.split()
				if 'INDEF' not in linesegs:
					savelines_index.append(i)

			savelines = source_lines[savelines_index]
#			print savelines

			if os.path.exists(stars_filename):
				os.remove(stars_filename)

			fid = open(stars_filename,'awt')
			for line in savelines:
				fid.write(line)
			fid.close()

		shutil.copy(stars_filename, iraf_daofind_outfile)

		#save ds9 region of found sources
                regfile = os.path.join(self.stars_dir,image_key.split('.')[0]+self.regfile_suffix)
		stars = np.loadtxt(stars_filename)
		create_ds9_region_file(stars, regfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', circle_radius = 10)
                starnum_this = len(open(stars_filename).readlines())

		return starnum_this


	def __get_star_number_rough_single_fistar_fitsh(self,image_key,ntime_above_bkg=3, which_dir = 'modified_image'):
		'''
		get the source list with fistar task from fitsh package
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		stars_filename = os.path.join(self.stars_dir,image_key.split('.')[0]+self.starfile_suffix)
		self.__delete_file_if_exist(stars_filename)

		peakthres = self.photometry_info[image_key]['bkg'] + ntime_above_bkg*np.sqrt(self.photometry_info[image_key]['bkg'])
		output_psf = os.path.join(self.stars_dir,image_key.split('.')[0]+'.psf')

		fistar = os.path.join(self.base.fitsh_dir, 'fistar')
		command = "%s -i %s -o %s -t %s --model elliptic -F x,y,magnitude,bg"%(fistar, input_image,stars_filename,peakthres)
		try:
                	os.system(command)
                except:
                	raise IOError("source finding on image %s fails..."%image_key)

		#fistar source detection algorithm find false targets in the edge region of the image if the image has bright edge
		stars = np.loadtxt(stars_filename)
		if self.filter_after_fistar:
			if self.edge_xl is not None:
				xmin = self.edge_xl
				stars = stars[stars[:,0]> xmin]
			if self.edge_xu is not None:
				xmax = self.edge_xu
				stars = stars[stars[:,0]< xmax]
			if self.edge_yl is not None:
				ymin = self.edge_yl
				stars = stars[stars[:,1]> ymin]
			if self.edge_yu is not None:
				ymax = self.edge_yu
				stars = stars[stars[:,1]< ymax]
			np.savetxt(stars_filename, stars, fmt='%6.2f %6.2f %5.2f %6.2f')

		#save ds9 region of found sources
                regfile = os.path.join(self.stars_dir,image_key.split('.')[0]+self.regfile_suffix)
		stars = np.loadtxt(stars_filename)
		create_ds9_region_file(stars, regfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', circle_radius = 10)

                starnum_this = len(open(stars_filename).readlines())

		return starnum_this




	def __get_star_number_rough_single_sfind_diapl(self,image_key, which_dir = 'modified_image'):
		'''
		get star number through sfind and save the star list file
		'''
		sfind_par = os.path.join(self.parafile_dir,'sfind.par')
		instrument_par = os.path.join(self.parafile_dir,'instrument.par')
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)

		stars_filename = os.path.join(self.stars_dir,image_key.split('.')[0]+self.starfile_suffix)
		self.__delete_file_if_exist(stars_filename)
		sfind = os.path.join(self.base.diapl_dir, 'sfind')
		command = "%s %s %s %s %s" %(sfind, sfind_par,instrument_par,input_image,stars_filename)
		print command
		try:
			os.system(command)
			print "1"
		except:
			raise IOError("source finding on image %s fails..."%image_key)

		regfile = os.path.join(self.stars_dir,image_key.split('.')[0]+self.regfile_suffix)
		stars = np.loadtxt(stars_filename)
		create_ds9_region_file(stars, regfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', circle_radius = 10)

		print "2"

		starnum_this = len(open(stars_filename).readlines())

		return starnum_this




	def __get_fwhm_bkg_single_fistar(self,key, which_dir = 'modified_image', time_out_max = 60):

		success, info = self.__fistar_action(key, which_dir = which_dir, time_out_max = time_out_max)
		return success, info


	def __fistar_action(self,key, which_dir = 'modified_image', time_out_max = 60):

		imgkey = key.split('.')[0]

		if which_dir == 'modified_image':
			img_file = os.path.join(self.modified_image_dir, key)
		else:
			img_file = self.images[key]

		print img_file
		hdu = fits.open(img_file)
		img_data = hdu[0].data
		print img_data.shape
		N1,N2 = img_data.shape
		print N1,N2


		from estimate_background import estimate_bkg
		bkgmean, bkgmedian, bkgstd = estimate_bkg(img_data)
		peakthreshold = bkgmedian+3*bkgstd
		source_file = os.path.join(self.stars_dir,imgkey + '.temp')

		fistar = os.path.join(self.base.fitsh_dir, 'fistar')
		fistar_command = "%s -i %s -o %s -t %s -s x -F id,x,y,bg,amp,fwhm"%(fistar, img_file,source_file,peakthreshold)
		command = Command(fistar_command)

		failure_signal = command.run(timeout=time_out_max)
		if not failure_signal:
			source_info = np.loadtxt(source_file)
			bkg = np.mean(source_info[:,3])
			fwhm = np.median(source_info[:,5])
			info = [bkg,fwhm]
			success = True
		else:
			print "dear friend, we can't get the result within %s s"%time_out_max
			info = None
			success = False

		return success, info


	def get_fwhm_bkg(self, method = 'diapl', which_dir = 'modified_image', time_out_max = 60):
		'''
		get fwhm and background

		INPUTs:
			method: 'diapl' or 'fistar' or 'borrow'
			which_dir: 'modified_image' or 'raw_image'
			time_out_max: the maximum time before time out for single image
		'''
		for key in self.photometry_info.keys():
			if self.photometry_info[key]['fwhm'] == 0:
				self.__get_fwhm_bkg(key,method =method, which_dir = which_dir, time_out_max=time_out_max)



	def renew_fwhm_bkg(self, method ='diapl', which_dir = 'modified_image', time_out_max = 60):
		'''
		get fwhm and background
		'''
		for key in self.photometry_info.keys():
                        self.__get_fwhm_bkg(key, method =method, which_dir = which_dir, time_out_max = time_out_max)



	def __get_fwhm_bkg(self,key,method='diapl', which_dir = 'modified_image', time_out_max = 60):
		'''
		get the fwhm and background of given image with imagekey 'key'

		INPUTS:
			key:
			method: 'fwhm' from diapl, or 'fistar' from fitsh, 'borrow' is using results from external reductions for example from WFCAM CASU reduction pipeline which give SEEING and SKYLEVEL
			which_dir: 'raw_image' or 'modified_image'
			time_out_max: time out threshold for this task
		'''

		if method == 'diapl':
			success,result = self.__get_fwhm_bkg_single_diapl(key, which_dir = which_dir, time_out_max = time_out_max)
			if success:
				fwhm  = round(np.float(result[-2]),2)
				bkg   = round(np.float(result[-3]),2)
				print key,fwhm,bkg
			else:
				fwhm  = self.anomaly_signs['fwhm']
				bkg   = self.anomaly_signs['bkg']
				self.photometry_info[key]['drop'] = 3
		elif method == 'fistar':
			success, result = self.__get_fwhm_bkg_single_fistar(key, which_dir = which_dir, time_out_max = time_out_max)
			if success:
				bkg  = round(result[0],2)
				fwhm = round(result[1],2)
				print key,fwhm,bkg
			else:
				fwhm  = self.anomaly_signs['fwhm']
				bkg   = self.anomaly_signs['bkg']
				self.photometry_info[key]['drop'] = 3
		elif method == 'borrow':
			result = self.__get_fwhm_bkg_single_borrow(key, which_dir =which_dir)
			bkg  = round(result[0],2)
			fwhm = round(result[1],2)
			print key,fwhm,bkg
		else:
			raise IOError("method %s not supported"%method)

		self.photometry_info[key]['fwhm'] = fwhm
		self.photometry_info[key]['bkg'] = bkg

		if self.photometry_info[key]['fwhm'] == -99.99:
			self.photometry_info[key]['drop'] = 3


	def __get_fwhm_bkg_single_borrow(self, image_key, which_dir='modified_image'):
		'''
		get fwhm and background by extracting external reduction results
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		fitsinfo = get_fits_info(input_image,['SEEING', 'SKYLEVEL'], extension=None)
		fwhm = fitsinfo['SEEING']
		bkg = fitsinfo['SKYLEVEL']
		outdata = [bkg, fwhm]

		return outdata

	def __get_fwhm_bkg_single_diapl(self,image_key, which_dir = 'modified_image', time_out_max = 60, monitorplot=False, deletemfile=True):
		'''
		get fwhm and background with function 'fwhm' from  diapl package
		'''
		templog = image_key.split('.')[0] + '_fwhm.log'
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)

		input_image_int16 = 'tmpimage_delete_after.fits'
		if os.path.exists(input_image_int16):
			os.remove(input_image_int16)

		flt2int = os.path.join(self.base.ccdproc_dir, 'flt2int')
		command = "%s %s %s"%(flt2int, input_image, input_image_int16)
		os.system(command)

		print input_image

		telcode = self.current_telescope
		fwhm_parfile_telspec = os.path.join(self.parafile_dir, 'fwhm_%s.par'%telcode)
		if os.path.exists(fwhm_parfile_telspec):
			fwhm_parfile = fwhm_parfile_telspec
		else:
			fwhm_parfile = os.path.join(self.parafile_dir, 'fwhm.par')

		fwhmstarsfile= os.path.join(self.stars_dir, image_key.split('.')[0]+'_fwhm.stars')

		fwhm = os.path.join(self.base.diapl_dir, 'fwhm_ping2')
		command_line = '%s %s %s %s>>%s' %(fwhm, fwhm_parfile, input_image_int16, fwhmstarsfile, templog)
		print command_line
		command = Command(command_line)
		failure_signal = command.run(timeout=time_out_max)

		if os.path.exists(input_image_int16) and deletemfile:
			os.remove(input_image_int16)

		if failure_signal:
			print "Sorry friend, we can't get the result within %s s"%time_out_max
			success = False
			data = None
		else:
			loglines = open(templog).readlines()
			fwhmline = loglines[-1]
			print fwhmline
			data = fwhmline.strip().split()[-5:]
			success = True

			if monitorplot:
				fwhmdata = np.loadtxt(fwhmstarsfile)
				fig = plt.figure(figsize=(5,5))
				plt.hist(fwhmdata[:,4], bins=30)
				plt.show()

		if os.path.exists(templog) and deletemfile:
			os.remove(templog)


		return success, data


	def remove_cosmic_ray(self, which_dir = 'modified_image', method='ccdproc',sigclip=5):
		'''
		remove cosmic ray
		INPUTS:
			which_dir: the directory containing the input images
			method: 'ccdproc' or 'cosmics'
		'''

		for key in self.photometry_info.keys():
			try:
				if not self.photometry_info[key]['rmCR']:
					self.__remove_cosmic_ray(key, which_dir = which_dir, method=method, sigclip=sigclip)
					self.photometry_info[key]['rmCR'] = 1
			except Exception as e:
				print e
				print "Remove cosmic ray failure on image %s" %key
				self.photometry_info[key]['rmCR'] = 0


	def __remove_cosmic_ray(self,image_key,which_dir = 'modified_image', method='ccdproc', output_image = None,single_use = False, image_gain = 1.0, image_readnoise=10.0, sigclip =5):
		'''
		INPUTS:
			image_key:
			which_dir:
			method: cosmics or ccdproc
			output_image:
			single_use:
			image_gain:
			image_read_noise:
			sigclip:
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		if output_image is None:
			output_image = os.path.join(self.modified_image_dir,image_key)

		# Build the object:
		if len(self.base.telescopes_info[self.current_telescope]['gain']):
			gain = float(self.base.telescopes_info[self.current_telescope]['gain'])
		else:
			gain = image_gain

		if len(self.base.telescopes_info[self.current_telescope]['readnoise']):
			readnoise = float(self.base.telescopes_info[self.current_telescope]['readnoise'])
		else:
			readnoise = image_readnoise

		if method == 'cosmics':
			# Read the FITS :
			array, header = cosmics.fromfits(input_image)
			# array is a 2D numpy array
			c = cosmics.cosmicsimage(array, gain=gain, readnoise=readnoise, sigclip = 5.0, sigfrac = 0.3, objlim = 5.0)
			# There are other options, check the manual...

			# Run the full artillery :
			c.run(maxiter = 4)

			# Write the cleaned image into a new FITS file, conserving the original header :
			if os.path.isfile(output_image):
				os.remove(output_image)
			cosmics.tofits(output_image, c.cleanarray, header)

			# If you want the mask, here it is :
			#mask_filename = ''
			#cosmics.tofits(mask_filename, c.mask, header)
			# (c.mask is a boolean numpy array, that gets converted here to an integer array)
		elif method == 'ccdproc':
			print input_image, output_image, gain, readnoise
			if os.path.exists(output_image):
				output_temp = "temp_clean_image.fits"
				if os.path.exists(output_temp):
					os.remove(output_temp)
				movedata = 1
			else:
				output_temp = output_image
				movedata = 0

			remove_CR_ccdproc_cosmicray_lacosmic(input_image, output_temp, gain, readnoise, sigclip=sigclip)
			if movedata:
				shutil.move(output_temp, output_image)


		else:
			raise IOError("method of %s not supported"%method)


		if single_use:
			self.photometry_info[image_key]['rmCR'] = 1



	def fix_bad_pixels(self,mode='vicinity', which_dir='raw_image', bad_threshold = 60000):
		'''
		replace the 'bad' pixels with 'good' values which determined by the 'fix_mode'

		Inputs:
			mode: how you want to fix the 'bad' pixels, available options are the followings,
				'vicinity', replace the 'bad' values with the values from the vicinity
				'bkg', replace the 'bad' values with the rough background
				'zero', replace the 'bad' values with zeros

		'''
		print "Fix bad pixels on image..."
		for key in self.photometry_info.keys():
			print key
			try:
				self.__fix_bad_pixels(key,mode=mode, which_dir = which_dir, bad_threshold = bad_threshold )
				self.photometry_info[key]['fixbad'] = 1
			except:
				print "Fix bad pixels failure on image %s" %key
				self.photometry_info[key]['fixbad'] = 0
		print "Fix bad pixels on image... Done!"


	def fix_negative_values_pixels(self, which_dir ='raw_image', fillmethod='bkg', verbose=1):
		'''
		fix the negative pixel value issue
		'''
		for img in self.images.keys():
			if fillmethod == 'bkg':
				bkgvalue = self.photometry_info[img]['bkg']
			self.__replace_negative_value(img, which_dir=which_dir, fillvalue=bkgvalue)
			if verbose:
				print "negative values in %s will be replaced with %s"%(img, bkgvalue)


	def __replace_negative_value(self, image_key, which_dir = 'raw_image', fillvalue=None, fillmethod='bkg', output_image=None):
		'''
		replace the negative values with fillvalue or the value derived from fillmethod
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		if output_image is None:
			output_image = os.path.join(self.modified_image_dir,image_key)

		hdu = fits.open(input_image)
		header = hdu[0].header
		data = hdu[0].data
		hdu.close()
		bad_pixels = self.__get_bad_pixels(data,0, mode='negative')
		data_new = self.__fill_bad_pixels(data,bad_pixels, fill_value=fillvalue, mode=fillmethod)

		hdu_new = fits.PrimaryHDU(data_new,header)
		hdulist = fits.HDUList([hdu_new])

		if os.path.isfile(output_image):
			os.remove(output_image)
		hdulist.writeto(output_image)


	def __fix_bad_pixels(self,image_key,mode='bkg',which_dir = 'raw_image',bad_threshold = 60000,  output_image = None,single_use = False):
		'''
		    Replace the saturated pixels with normal value which determined by the 'mode'
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)

		if output_image is None:
			output_image = os.path.join(self.modified_image_dir,image_key)

		hdu = fits.open(input_image)
		header = hdu[0].header
		data = hdu[0].data
		hdu.close()
		bad_pixels = self.__get_bad_pixels(data,bad_threshold)
		data_new = self.__fill_bad_pixels(data,bad_pixels,mode=mode)

		hdu_new = fits.PrimaryHDU(data_new,header)
		hdulist = fits.HDUList([hdu_new])
		if os.path.isfile(output_image):
			os.remove(output_image)
		hdulist.writeto(output_image)

		if single_use:
			self.photometry_info[image_key]['fixbad'] = 1



	def __get_bad_pixels(self,image,threshold, mode='saturation'):
		'''
		bad pixels have value higher than saturation threshold in saturation mode or have negative values in negative mode
		INPUTS:
			image:
			threshold:
			mode: saturation or negative
		'''
		mask = np.zeros(image.shape)
		if mode == 'saturation':
			mask[np.where(image > threshold)]=1
		elif mode == 'negative':
			mask[np.where(image<threshold)] =1
		else:
			raise ValueError("invalid input for mode...")

		return mask



	def __fill_bad_pixels(self,image,bad_pixels,fill_value=None, mode='bkg'):
		'''
		fill the bad pixels with the desired value which is determined by option 'mode'
		mode = 'bkg' or 'vicinity'
		when mode = 'bkg' the bad_pix will be filled by the value of the image background
		when mode = 'vicinity' the bad_pix will be filled by the average value ofthe surrounding pixels(time consuming!!)
		when mode = 'zero' the bad value will be replaced with zeros
		'''
		M,N = image.shape
		X,Y = np.where(bad_pixels ==1)

		if fill_value is not None:
			image[X,Y] = fill_value
		else:
			if mode=='vicinity':
				shifts = self.__get_surrounding(radius = 1)
				#print shifts
				temp = np.ma.zeros((len(shifts),M,N))

				for i,[x,y] in enumerate(shifts):
				    xlo = 0+x
				    ylo = 0+y
				    xho = M+x
				    yho = N+y
				    img_shift = self.__cut_or_extend_imagedata(image,xlo,ylo,xho,yho)
				    badpix_shift = self.__cut_or_extend_imagedata(bad_pixels,xlo,ylo,xho,yho)
				    masked_image = np.ma.array(img_shift,mask = badpix_shift)
				    #print type(masked_image)
				    temp[i] = masked_image
				#print type(temp)
				mean_image=np.mean(temp,axis = 0)
				for x,y in zip(X,Y):
				    #print image[x,y],mean_image[x,y]
				    image[x,y]=mean_image[x,y]
				    #print i
			if mode=='bkg':
			    	bkg,std =self.__get_rough_background(image)
			        image[X,Y]=bkg

			if mode == 'zero':
			    	image[X,Y] = 0

		return image


	def __get_rough_background(self,image):
		'''
		get the background intensity of the given image
		'''
		imgbkg = image.reshape(1,np.size(image))
		cri = np.mean(imgbkg) + np.std(imgbkg)
		imgbkg = imgbkg[imgbkg<cri]
		#print imgbkg.shape
		background,bkgfilter,std = sigma_clipping(imgbkg, sig=3, meanfunc=np.median)
		#show_histogram(background,nbins = 100)
		bkg = np.mean(background)

		return bkg,std



	def __get_surrounding(self,radius = 1):
		'''
		get the pixels which are inside the circle with center at [0,0] and radius = radius
		'''
		shifts = []

		for x in np.arange(-int(radius),int(radius)+1):
		    for y in np.arange(-int(radius),int(radius)+1):
		        if x**2+y**2 <= radius**2:
		    		shifts.append([x,y])

		return shifts


	def __cut_or_extend_imagedata(self,imagedata,xlo,ylo,xho,yho):
		'''
		cut or extend the input image data
		we can get the size of the image from the image coordinate in which the left bottom pix has (0,0)
		and the right top pix has (M-1,N-1) where (M,N) = imagedata.shape
		the image coordinates of the left bottom pix and right top pix in the output image are (xlo,ylo),(xho,yho) respectively
		(xlo,ylo),(xho,yho) are give in the input image frame
		note: the record of data on CCD is not the same as the normal array. the x axis is  along the bottom to top direction and
		the y axis is along the left to right direction
		'''
		data = imagedata
		temp = np.zeros((xho-xlo,yho-ylo),dtype=data.dtype.name)
		xtemp1,ytemp1 = temp.shape

		x0,y0=tuple([0,0])
		x1,y1=data.shape

		if xlo<=0:
		    xf0 = x0
		    xt0 = x0-xlo
		else:
		    xf0 = xlo
		    xt0 = x0

		if ylo<=0:
		    yf0 = y0
		    yt0 = y0 - ylo
		else:
		    yf0 = ylo
		    yt0 = y0

		if xho>=x1:
		    xf1 = x1
		    xt1 = x1-xlo
		else:
		    xf1 = xho
		    xt1 = xtemp1

		if yho>=y1:
        	    yf1 = y1
        	    yt1 = y1-ylo
		else:
        	    yf1 = yho
        	    yt1 = ytemp1

            	temp[xt0:xt1,yt0:yt1]=data[xf0:xf1,yf0:yf1]
	        new_data = temp
        	return new_data



	def wfcam_prepare_image(self, extension_want=None, verbose=0):
		'''
		WFCAM reduced image has 4 science images in 4 seperate extensions. Select the one containing the target and write to individual image
		'''
		for img in self.photometry_info.keys():
			input_image = self.images[img]
			output_image =  os.path.join(self.modified_image_dir, img)
			if os.path.exists(output_image):
				continue

			self.__get_wfcam_target_image(input_image, output_image,extension_want=extension_want, verbose=verbose)



	def __get_wfcam_target_image(self, input_image, output_image, extension_want=None, verbose=0):
		'''
		the image containing the target should be closest to the target coordinate
		'''
		if extension_want is not None:
			self.__create_new_image_from_selected_extension(input_image, output_image, extension_want)
		else:
			if self.sn_ra_world_deg is None or self.sn_dec_world_deg is None:
				print "the target coordinate"

			snrarad = self.sn_ra_world_deg / 180.0 * np.pi
			sndecrad = self.sn_dec_world_deg / 180.0 * np.pi

			distances = np.array([])
			hdu = fits.open(input_image)
			for i in (np.arange(4)+1):
				hdr = hdu[i].header
				data = hdu[i].data
				NX,NY = data.shape
				imagecoor = np.array([NX/2, NY/2]).reshape(1,2)
				center_coor = self.__image2world_fitshead(hdr, imagecoor)
				cradeg,cdecdeg = center_coor[0]
				crarad = cradeg/180.0*np.pi
				cdecrad = cdecdeg/180.0*np.pi
				raraddiff = np.abs(snrarad - crarad)
				distance = np.arccos(np.sin(sndecrad)*np.sin(cdecrad) + np.cos(sndecrad)*np.cos(cdecrad)*np.cos(raraddiff))	#compute the great-circle distance
				distances = np.append(distances, distance)
				if verbose:
					print i, distance

			extension_want = np.where(distances==np.min(distances))[0][0] + 1
			if verbose:
				print "the image containing the target for %s is in extension %s"%(input_image, extension_want)
			self.__create_new_image_from_selected_extension(input_image, output_image, extension_want)



	def __create_new_image_from_selected_extension(self, input_image,output_image, want_hduindex, hdr_kw_del = ['XTENSION'], output_verify='ignore'):
		'''
		get the wanted extension data and header, delete specified header keywords, and write to output image
		'''
		hdu = fits.open(input_image)
		hdudata = hdu[want_hduindex].data
		hduhdr  = hdu[want_hduindex].header

		for delkw in hdr_kw_del:
			hduhdr.pop(delkw)

		newimg = fits.PrimaryHDU(data=hdudata, header=hduhdr)
		newimg.writeto(output_image, output_verify=output_verify)


	def image_flt2int(self, which_dir = 'modified_image', nocheck=False):
		'''
		convert data type from float to integer for those with bitpix == -32
		'''
		try:
			for key in self.photometry_info.keys():
				if nocheck:
					self.__flt2int(key, which_dir = which_dir )
				else:
					if self.photometry_info[key]['bitpix'] == -32 and (not self.photometry_info[key]['flt2int']):
						self.__flt2int(key, which_dir = which_dir )
		except:
			print "flt2int not available..."

		print "CCDProc flt2int... Done!"




	def __flt2int(self, image_key, which_dir='modified_image'):
        	'''
        	convert the data type from float32 to int to meet the requirement of some functions in diapl
        	'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		output_image_temp = os.path.join(self.modified_image_dir, 'int16_temp.fits')
		output_image = os.path.join(self.modified_image_dir,image_key)

		flt2int = os.path.join(self.base.ccdproc_dir, 'flt2int')
        	command = '%s %s %s'%(flt2int,input_image,output_image_temp)
		try:
			os.system(command)
			status = 1
			if os.path.isfile(output_image):
				os.remove(output_image)
			os.rename(output_image_temp, output_image)

		except:
			print "skip image %s in flt2int"%image_key
			status = 0

		self.photometry_info[image_key]['flt2int'] = status



	def __renew_bitpix_info(self,image_key):
		img_abs = self.images[image_key]
                bitpix_key = self.base.telescopes_info[self.current_telescope]['bitpixkey']
                self.photometry_info[image_key]['bitpix'] = get_bitpix(img_abs,bitpix_key,self.current_telescope)

	def __renew_flt_info(self,image_key):
		img_abs = self.images[image_key]
		flt_key = self.base.telescopes_info[self.current_telescope]['fltkey']
		try:
			flt = get_filter(img_abs,flt_key,self.current_telescope)
			self.photometry_info[image_key]['flt'] = flt
		except:
			self.photometry_info[image_key]['drop'] = 12

	def __renew_airmass_info(self, image_key):
		img_abs = self.images[image_key]
		airmass_key = self.base.telescopes_info[self.current_telescope]['airmasskey']
                self.photometry_info[image_key]['airmass'] = get_airmass(img_abs,airmass_key,self.current_telescope)


	def __get_airmass_astroobs_single(self, imgkey):
		obsJD = self.photometry_info[imgkey]['obstime']
		airmass_ret = get_airmass_given_time_target_observatory(obsJD, self.current_telescope, self.sn_ra_str, self.sn_dec_str)
		return airmass_ret

	def get_airmass(self):
		airmass_key = self.base.telescopes_info[self.current_telescope]['airmasskey']
		if airmass_key != '':
			for img in self.images.keys():
				img_abs = self.images[img]
				self.photometry_info[img]['airmass'] = get_airmass_telescope(img_abs,airmass_key,self.current_telescope)
		else:
			for img	in self.images.keys():
				airmass_ret = self.__get_airmass_astroobs_single(img)
				self.photometry_info[img]['airmass'] = airmass_ret

	def __renew_obstime_info(self,image_key):
		img_abs = self.images[image_key]
		obstime_key = self.base.telescopes_info[self.current_telescope]['obstimekey']
                self.photometry_info[image_key]['obstime'] = get_obstime(img_abs,obstime_key,self.current_telescope)

	def __renew_exptime_info(self,image_key):
		img_abs = self.images[image_key]
		exptime_key = self.base.telescopes_info[self.current_telescope]['exptimekey']
                self.photometry_info[image_key]['exptime'] = get_exptime(img_abs,exptime_key,self.current_telescope)


	def renew_image_info(self,image_key=None,renewall=False,):
		'''
		renew image information: filter,obstime and bixpix
		'''
		if renewall:
			for key in self.photometry_info.keys():
				self.__renew_flt_info(key)
                                self.__renew_obstime_info(key)
				self.__renew_bitpix_info(key)
				self.__renew_exptime_info(key)

		elif image_key is not None and image_key in self.images.keys():
			self.__renew_flt_info(image_key)
			self.__renew_obstime_info(image_key)
			self.__renew_bitpix_info(image_key)
			self.__renew_exptime_info(key)
		else:
			print "no valid input for this action..."

	def save(self, outfile=None):
		'''
		pickle the current photometry instance to self.result_dir
		'''
		if outfile is None:
			outfile = 'photometry.pickle'
		savetofile = os.path.join(self.result_dir, outfile)

		if os.path.exists(savetofile):
			print "%s already exists; leave now and let you check!"
			return

		print "the photometry object will be pickled to %s"%savetofile
		f = open(savetofile, 'w')
		pickle.dump(self, f)
		f.close()


	def save_results(self,save_filename = None):
		'''
		convert photometry_info dictionary to table and save.
		'''
		self.__dict2table()
		self.result_table_newsaved.sort('obstime')

		if save_filename is None:
			save_filename = self.result_table_file
		self.__delete_file_if_exist(save_filename)
		Table.write(self.result_table_newsaved,save_filename,format='ascii.fixed_width')


	def __dict2table(self):
		'''
		self.photometry_info --> self.result_table_newsaved
		'''
		self.result_table_newsaved = Table(names=self.photometry_info_keys,dtype=self.photometry_info_dtypes)
        	for img in self.photometry_info.keys():
            		self.result_table_newsaved.add_row(self.photometry_info[img].values())
		for colname in ['relmag', 'relmagerr','calmag','calmagerr']:
			self.result_table_newsaved[colname] = np.round(self.result_table_newsaved[colname], 3)


	def __load_old_results(self,info_table,info_dict):
		for key in info_dict.keys():
			for colname in info_table.colnames:
				name_key = np.where(info_table['name'] == key)[0]
				if len(name_key):
					name_key_indice = name_key[0]
					info_dict[key][colname] = info_table[name_key_indice][colname]


	def __load_source_regions_to_ds9(self,ds9_object,input_reg,radius=10,color='green',width=1):
		'''
		Inputs:
			ds9_object:
			input_reg: ds9 region file or position list with two element (x,y)
			radius:
			color:
			width:
		'''
		if isinstance(input_reg,str) and os.path.isfile(input_reg):
			ds9_object.set('regions load %s'%input_reg)
		elif isinstance(input_reg,list):
			x = input_reg[0]
			y = input_reg[1]
			region_set = 'regions command {circle %s %s %s #color=%s width=%s}'%(x,y,radius,color,width)
			print region_set

			ds9_object.set(region_set)
		else:
			try:
				xy = input_reg.ravel()
				x = xy[0]
                		y = xy[1]
                		ds9_object.set('regions command {circle %s %s %s #color=%s width=%s}'%(x,y,radius,color,width))
			except:
				raise IOError("invalid input for ds9 region")


	def load_dophot_photret_to_ds9(self, image_key, which_dir='raw_image', xcol=2, ycol=3, text1col=4, text2col=1, filtercol=1, filtertype='eq',filtervalue=1, coortype = 'image', color = 'green', width=1, radius=15, newds9=True):
		'''
		load phot object with different caption choice
		INPUTS:
			xcol, ycol: define the object coordinate
			text1col, text2col: the text columns to show
			filtercol: if None then no filtering
			filtertype: 'lt', less than; 'eq', equal to; 'gt', greater than
			filtervalue: the criteria for filtering
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
                photdir = self.psf_photometry_dir
                imgkey = image_key.split('.')[0]
                photret_file = os.path.join(photdir, imgkey+'.out')
                photret = np.loadtxt(photret_file)

		if filtercol is not None:
			if filtertype == 'eq':
				mask = photret[:,filtercol] == filtervalue
			elif filtertype == 'lt':
				mask = photret[:,filtercol] < filtervalue
			elif filtertype == 'gt':
				mask = photret[:,filtercol] > filtervalue
			else:
				raise ValueError('not support filter type %s'%filtertype)
			photret = photret[mask]

		regionfile = os.path.join(photdir, imgkey+"_photret.reg")
		create_ds9_region_file(photret, regionfile, clobber = True,  x_col=xcol, y_col=ycol, x_offset=0, y_offset=0, coordinate_type = coortype, radec_deg = True, circle_radius = radius, color=color, width=width, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=text1col, textcol2=text2col)

		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		d.set('regions load %s'%regionfile)



        def load_stars_xy_mag_to_ds9(self,image_key,photometry_method = 'apphot', which_dir='raw_image', offset= 0, mag_cut = 1, xl = None, xh=None, yl=None, yh = None, coortype = 'image', color = 'green', width=1, radius=15, image_physical_diff_x = 0, image_physical_diff_y = 0, newds9=True):
		'''
		load star region ans magnitude measurements to the image in DS9

		Inputs:
			image_key:
			photometry_method:
			offset: add an offset to magnitudes measurements
			mag_cut: bright %(mag_cut*100) will be displayed
			image_physical_diff_x:
			image_physical_diff_y:
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
                if photometry_method == 'apphot':
                    photdir = self.aperture_photometry_dir
		    photret_suffix = self.apphot_ret_file_suffix
                elif photometry_method == 'psfphot':
                    photdir = self.psf_photometry_dir
		    photret_suffix = self.psfphot_ret_file_suffix
                else:
                    raise IOError("input for photometry_method is not supported...")

                imgkey = image_key.split('.')[0]
                photret_file = os.path.join(photdir, imgkey+ photret_suffix)
                xymags = np.loadtxt(photret_file)

		if xl is not None:
			xymags = xymags[xymags[:,0] > xl,:]
		if xh is not None:
			xymags = xymags[xymags[:,0] < xh,:]
		if yl is not None:
			xymags = xymags[xymags[:,1] > yl,:]
		if yh is not None:
			xymags = xymags[xymags[:,1] < yh,:]

		xymags = xymags[np.argsort(xymags[:,2]),:]
		N = len(xymags)
		N_want = int(np.ceil(N*mag_cut))
		mask = range(N_want)
		xymags_show = xymags[mask,:]

		regionfile = os.path.join(photdir, imgkey+"_xymag.reg")
		self.ds9D  = self.__load_stars_xy_mag_to_ds9(input_image, xymags_show, regionfile, mag_offset = offset, coortype = coortype, color =color, width=width, radius=radius, newds9=newds9)


	def __load_stars_xy_mag_to_ds9(self, input_image, xymags, regionfile, mag_offset  = 0, coortype = 'image', color = 'green', width=1, radius=15, newds9=True):
		'''
		load object list to ds9
		'''
		if mag_offset != 0:
			xymags[:,2] = xymags[:,2] + mag_offset
		create_ds9_region_file(xymags, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = coortype, radec_deg = True, circle_radius = radius, color=color, width=width, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)

		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		d.set('regions load %s'%regionfile)

                return d



	def load_matched_stars_to_ds9(self, img, matched_source_file = None, ref_reg_file = None, input_reg_file = None):
		'''
		load two set sources which are matched one to one.

		INPUTS:
			matched_source_file:	first and second columns are images coordinates of ref sources; and fifth and sixth columns are images coordinates of input sources
			ref_reg_file:
			input_reg_file:
		'''

		#check inputs
		if ref_reg_file is not None and input_reg_file is not None:
			if os.path.exists(ref_reg_file) and os.path.exists(input_reg_file):
				regionfile_provided = True
			else:
				raise IOError("input region file %s and/or %s doesn't exist..."%(ref_reg_file, input_reg_file))
		else:
			regionfile_provided = False

		if matched_source_file is None and (not regionfile_provided):
			raise IOError("No valid region file provided...")

		#explicit ds9  region files have higher priority
		if regionfile_provided:
			self.__load_matched_stars_to_ds9(img, ref_reg_file, input_reg_file)
		else:
			ref_reg_file_temp = 'temp_ref_sources.reg'
			input_reg_file_temp = 'temp_input_source.reg'
			input_sources = np.loadtxt(matched_source_file)

			create_ds9_region_file(input_sources, ref_reg_file_temp,   clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 20)
			create_ds9_region_file(input_sources, input_reg_file_temp, clobber = True,  x_col=4, y_col=5, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

			self.__load_matched_stars_to_ds9(img, ref_reg_file_temp, input_reg_file_temp)




	def __load_matched_stars_to_ds9(self, input_image, ref_reg_file, input_reg_file, newds9=True):
		'''
		load ds9 region files 'ref_reg_file' and 'input_reg_file' to image 'img'
		'''
		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		print "reference with large circle"
		self.__load_source_regions_to_ds9(d, ref_reg_file)
		self.__load_source_regions_to_ds9(d, input_reg_file)


	def load_std_stars_region_to_ds9_flt(self,flt,load_mags=False, on_image_after_astrometry = False, newds9=True, updatetable=1):
		'''
		load standard reference stars to the template image
		Here standard stars are all available ones before match with sources measured on the science image
		INPUTs:
			flt:
			load_mags: if True, then load the magnitudes
			on_image_after_astrometry: if True, then the image in ./warehouse/template/cal_xxx.fits will be used
		'''
		self.__find_template_imagekey(updatetable=updatetable)
		self.__get_stdref_filenames()
		if not on_image_after_astrometry:
			image_key = self.templates[flt]
        		input_image = self.images[image_key]
		else:
			input_image = self.templates_after_astrometry[flt]

		d = self.__display_image_with_ds9(input_image, newds9=newds9)

		regfile = os.path.join(self.std_ref_dir,flt+self.regfile_suffix)
		if not os.path.exists(regfile):
			raise ValueError("region file %s not exist"%regfile)
		self.__load_source_regions_to_ds9(d,regfile)

		return d


	def __add_text_to_ds9(self,d,x,y,text,color='red',width=2, physical_image_offset = True):
		'''
		add 'text' to ds9 frame d at position (x,y)

		Notes:
			Use self.physical_image_offset_x and self.physical_image_offset_y to deal with the difference between image coordinate and physical coordinate
		'''
		if physical_image_offset:
			x_diff = self.physical_image_offset_x
			y_diff = self.physical_image_offset_y
		else:
			x_diff = 0
			y_diff = 0

		x_physical = x + x_diff
		y_physical = y + y_diff
		d.set('regions command {text %s %s #text="%s" color="%s" width=%s}'%(x_physical,y_physical,text,color,width))


	def load_fwhm_stars_to_ds9(self, image_key, which_dir='raw_image', newds9=True):
		'''
		load the star list from which FWHM value derived to ds9
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		fwhmstarsfile= os.path.join(self.stars_dir, image_key.split('.')[0]+'_fwhm.stars')
		source_regfile = os.path.join(self.stars_dir, image_key.split('.')[0] + '_fwhmstars.reg')
		fwhmstars = np.loadtxt(fwhmstarsfile)
		create_ds9_region_file(fwhmstars, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		self.__load_source_regions_to_ds9(d,source_regfile)

		return d

	def load_source_regions_to_ds9(self,image_key,regtype = 'prephot',whichone='staronly', which_dir='raw_image', newds9=True):

		'''
		load source list to ds9

		Input:
			regtype: 'prephot' or 'postphot'
			whichone:'staronly' or 'all', only work for psfphotmetry
				 'staronly' mean identification restricted to stars
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		d = self.__display_image_with_ds9(input_image, newds9=newds9)

		if regtype == 'prephot':
			source_regfile = os.path.join(self.stars_dir, image_key.split('.')[0] + self.regfile_suffix)
		elif regtype == 'postphot':

			if self.photometry_method is None:
				self.photometry_method = raw_input('Please input the photometry method used(apphot/psfphot):')

			if self.photometry_method == 'apphot':
				source_regfile = os.path.join(self.aperture_photometry_dir, image_key.split('.')[0] + self.regfile_suffix)
				if not os.path.isfile(source_regfile):
					raise IOError("region file %s doesn't exist"%source_regfile)

			elif self.photometry_method == 'psfphot':
				if whichone == 'staronly':
					source_regfile = os.path.join(self.psf_photometry_dir, image_key.split('.')[0] + self.regfile_suffix)
				else:
					source_regfile = os.path.join(self.psf_photometry_dir, image_key.split('.')[0] + '_all'+ self.regfile_suffix)

				if not os.path.isfile(source_regfile):
					raise IOError("region file %s doesn't exist"%source_regfile)
			else:
				raise IOError('Invalid input for photometry_method')
		else:
			raise IOError('Invalid input for regtype')

		self.__load_source_regions_to_ds9(d,source_regfile)

		x = self.photometry_info[image_key]['x']
		y = self.photometry_info[image_key]['y']
		xy = [x,y]
		if x != 0.0 and y != 0.0:
			self.__load_source_regions_to_ds9(d,xy,radius =20, color='red', width=2)

		return d



	def display_multi_images_and_highlight_sn(self,flts = None,jd_start=None,jd_end = None, display_drop_with_flag = 0, updatetable=1):
		'''

		display images with sn position labeled as well obstime and flt information shown
		Default is all registered images will be displayed unless flts/jd_start/jd_end is given
		'''
		if updatetable:
			self.__dict2table()
		flts_valid = np.unique(self.result_table_newsaved['flt'])

		if flts is None:
			img_table_flts = self.result_table_newsaved
		else:
			img_table_flts = Table()

			if isinstance(flts,str):
				flts = [flts]
			for flt in flts:
				if flt not in flts_valid:
					raise KeyError("no image taken with %s filter"%flt)
				table_toadd = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt)
				img_table_flts = vstack([img_table_flts, table_toadd])

		img_table_flts_jd = img_table_flts
		if jd_start is not None:
			img_table_flts_jd = self.__select_rows_from_table(img_table_flts_jd, 'obstime',jd_start,mode='gt' )
 		if jd_end is not None:
			img_table_flts_jd = self.__select_rows_from_table(img_table_flts_jd, 'obstime',jd_end,mode='lt' )

		img_table_flts_jd_dropflag = self.__select_rows_from_table(img_table_flts_jd, 'drop', display_drop_with_flag)
		img_keys = img_table_flts_jd_dropflag['name']

		self.__display_multi_images_and_highlight_sn(img_keys)


	def __display_multi_images_and_highlight_sn(self,img_keys):
		'''
		display multiple images in ds9
		'''
		pioneer_img = img_keys[0]

		d = self.display_image_with_ds9(pioneer_img)
		x = self.photometry_info[pioneer_img]['x']
		y = self.photometry_info[pioneer_img]['y']
		xy_reg = [x,y]
		self._photometry__load_source_regions_to_ds9(d, xy_reg, radius=25, color='red', width=2)

		x_text = x + 100
		y_text = y + 100
		jd = self.photometry_info[pioneer_img]['obstime']
		t = Time(jd,scale='utc',format='jd')
		t_isot = t.isot
		fwhm_this = self.photometry_info[pioneer_img]['fwhm']
		flt_this = self.photometry_info[pioneer_img]['flt']
		text  = flt_this + '@' + t_isot + 'with fwhm=%s'%fwhm_this
		self._photometry__add_text_to_ds9(d,x_text,y_text,text,color='green',width=2,)

		for img in img_keys:
		        if img == pioneer_img:
		                continue
		        d.set('frame new')
		        current_img = self.images[img]
		        d.set('file %s'%current_img)
		        d.set('zoom to fit')

		        x = self.photometry_info[img]['x']
		        y = self.photometry_info[img]['y']
		        xy_reg = [x,y]
		        self._photometry__load_source_regions_to_ds9(d, xy_reg, radius=25, color='red', width=2)

		        x_text = x + 100
		        y_text = y + 100
		        jd = self.photometry_info[img]['obstime']
		        t = Time(jd,scale='utc',format='jd')
		        t_isot = t.isot
			fwhm_this = self.photometry_info[img]['fwhm']
		        flt_this = self.photometry_info[img]['flt']
		        text  = flt_this + '@' + t_isot + ' with fwhm=%s'%fwhm_this
		        self._photometry__add_text_to_ds9(d,x_text,y_text,text,color='green',width=2,)


	def __display_image_with_ds9(self,input_image, newds9=True):
		'''
		display input_image in ds9
		'''
		if newds9:
             		d = pyds9.DS9()
               		d.set('fits %s'%input_image)
               		d.set('zoom to fit')
                	d.set('scale zscale')
			self.ds9D = d
		else:
			if not self.ds9D:
				raise ValueError("no available ds9 instance in self.ds9D")
			d = self.ds9D

        	return d


	def display_image_with_ds9(self,image_key, which_dir = 'raw_image', newds9=True):
        	'''
		display image with pyds9
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		return d


	def __get_internal_image(self, image_key, which_dir='raw_image'):
		'''
		GET the corresponding filename of different kind defined by which_dir
		INPUTS:
			which_dir, 'raw_iamge':
				   'modified_image':
				   'interp_image':
				   'conv_image':
				   'subtraced_image':
				   'astrometry_image':
				   'psfout_image': PSF photometry output residual image
				   'psf_subras': emperical psf subtraster from dophot_C
		'''
		if which_dir == 'raw_image':
			input_image = self.images[image_key]
		elif which_dir == 'modified_image':
			input_image = os.path.join(self.modified_image_dir, image_key)
		elif which_dir == 'subtracted_image':
			input_image = os.path.join(self.subtraction_dir, 'sub_'+image_key)
		elif which_dir == 'interp_image':
			input_image = os.path.join(self.subtraction_dir, 'interp_'+image_key)
		elif which_dir == 'conv_image':
			input_image = os.path.join(self.subtraction_dir, 'conv_'+image_key)
		elif which_dir == 'astrometry_image':
			input_image = os.path.join(self.template_dir,'cal_'+image_key)
		elif which_dir == 'psfout_image':
			input_image = os.path.join(self.psf_photometry_dir, 'out_'+image_key)
		elif which_dir == 'psf_subras':
			input_image = os.path.join(self.psf_photometry_dir, image_key.split('.')[0] + '_psf.fits')
		else:
			raise IOError('Non valiad input for which_dir...')

		if not os.path.exists(input_image):
			print "warning: the request image does not exist"

		return input_image


	def label_single_xy_region(self,image,xy):
		d = self.__display_image_with_ds9(image)
		d = self.__load_source_regions_to_ds9(d,xy,radius=15,color='red',width=2)
		return d


	def reln_image_to_datadir(self, verbose=0):
		'''
		re-link the registered images in photometry result table to the workplace folder
		'''
		photret_legacy = self.result_table
		realimgs = photret_legacy['realimg'].data
		dstimgs_retfile = photret_legacy['name'].data

		from_dir = self.repository_dir_current_sn
		to_dir = self.raw_image_dir
		dstimgs_old_all = os.listdir(to_dir)

		for img in dstimgs_old_all:
			dstimg_old_abs = os.path.join(to_dir,img)
                        os.remove(dstimg_old_abs)

		for srcimg,dstimg in zip(realimgs,dstimgs_retfile):
			img_to = os.path.join(to_dir,dstimg)
			img_from = os.path.join(from_dir,srcimg)
			command = "ln -s %s %s"%(img_from,img_to)
			if verbose:
				print command
			os.system(command)


	def ln_image_to_datadir(self,from_dir,to_dir):
		'''
		spaces in filenames can be annoying (cause unexpected errors in bash shell) and remove them first
		'''
		self.__filename_modification_remove_spaces(from_dir)

		ori_images = os.listdir(from_dir)
		lned_images = os.listdir(to_dir)
		num = len(lned_images)

		if num > 999:
			numdig = 4
		else:
			numdig = 3

		lned_images_realpath = [os.path.realpath(os.path.join(to_dir,img)) for img in lned_images]
		for image in ori_images:
			image_abs = os.path.realpath(os.path.join(from_dir,image))
			if image_abs not in lned_images_realpath:
				obs_indicator = '0'*(numdig-len(str(num)))+str(num)
				new_lned_image = obs_indicator + '.fits'
				from_image = os.path.join(from_dir,image)
				to_image = os.path.join(to_dir,new_lned_image)

				print "ln -s %s %s"%(from_image,to_image)
				command = "ln -s %s %s"%(from_image,to_image)
				os.system(command)
				num = num+1


	def __filename_modification_remove_spaces(self,file_dir):
		'''
		remove space in filename
		'''
		filenames = os.listdir(file_dir)
		filenames_new = [filename.replace(" ","") for filename in filenames]
		for fromfile,tofile in zip(filenames,filenames_new):
			if fromfile != tofile:
				src = os.path.join(file_dir,fromfile)
				dst = os.path.join(file_dir,tofile)
				os.rename(src,dst)


	def get_images_dict(self,data_dir):
		'''
		return the images info dict containing the relative image names as keys and absolute names as values
		'''
		images_dict = OrderedDict()
		filelist = os.listdir(data_dir)
		for img in filelist:
			images_dict[img] = os.path.join(data_dir,img)
		return images_dict


	def get_std_obs_match_filenames(self):
		'''
		prepare the filenames for calibration
		'''
		for img_key in self.images.keys():
			flt = self.photometry_info[img_key]['flt']
			stdobs_matchfile = os.path.join(self.std_ref_dir, 'std_%s_%s.txt'%(flt, img_key.split('.')[0]))
			stdobs_trans_coeffile = os.path.join(self.std_ref_dir, 'std_%s_%s_trans_coef.txt'%(flt, img_key.split('.')[0]))
			self.stdobs_match_files[img_key] = stdobs_matchfile
			self.stdobs_transcoef_files[img_key] = stdobs_trans_coeffile


	def get_secondary_stds(self, flt, secondary_source, secondary_photmethod='apphot', radec_format = False, xy_format=True, cut_radius=None):
		'''
		When primary standard stars are not enough for calibration purpose, we can use images from other resources to link the primary standard stars with stars in input image

		INPUTS:
			flt:
			secondary_source: the telescope/instrument name from which the secondard standard stars are from
			photmethod: photometry method
			radec_format: the coordiante in ra, dec format
			xy_format: the coordinate in image x,y format
			cut_radius: if not None, only return stars with distance to the supernova on the secondary-source image within 'cut_radius'
		'''
		sn_current = self.current_sn
		tel_current = self.current_telescope
		apphot_dir = self.aperture_photometry_dir
		psfphot_dir = self.psf_photometry_dir
		photmethod_current = self.photometry_method

		second_phot_infofile = self.result_table_file.replace(tel_current, secondary_source)
		if photmethod_current != secondary_photmethod:
			if secondary_photmethod == 'psfphot':
				second_phot_infofile = second_phot_infofile[:-4] + '_PSF' + second_phot_infofile[-4:]
			else:
				second_phot_infofile = second_phot_infofile.replace('_PSF', '')

		if not os.path.exists(second_phot_infofile):
			print "Sorry, you request secondary standard photometry %s not exist"%second_phot_infofile
			return

		sdphotinfo = Table.read(second_phot_infofile, format='ascii.fixed_width')
		sdtpls = sdphotinfo[sdphotinfo['template'] == 1]

		sdtpl_flt = sdtpls[sdtpls['flt'] == flt]
		if len(sdtpl_flt) != 1:
			print "template image not found in secondary source, or more than one found..."
			return

		sdtplflt = sdtpl_flt[0]
		sdimg = sdtplflt['name']
		if secondary_photmethod == 'psfphot':
			sdtplphot_dir = psfphot_dir.replace(tel_current, secondary_source)
			sdtplphot_file = os.path.join(sdtplphot_dir, sdimg[:-5]+self.psfphot_ret_file_suffix)
		elif secondary_photmethod == 'apphot':
			sdtplphot_dir = apphot_dir.replace(tel_current, secondary_source)
			sdtplphot_file = os.path.join(sdtplphot_dir, sdimg[:-5]+self.apphot_ret_file_suffix)
		else:
			raise ValueError("invalid input for secondary_photomethod")

		sdx = sdtplflt['x']
		sdy = sdtplflt['y']
		sdoffset = sdtplflt['calmag'] - sdtplflt['instmag']
		sdtplphot = np.loadtxt(sdtplphot_file)
		sdtplphot[:,2] = sdtplphot[:,2] + sdoffset
		if cut_radius is not None:
			mask = np.sqrt( (sdtplphot[:,0] - sdx)**2 + (sdtplphot[:,1] -sdy)**2 ) < cut_radius
			sdtplphot = sdtplphot[mask]
		if radec_format:
			print "this feature still under construction"
		if xy_format:
			sdstdfile = "%s_%s_stds.txt"%(secondary_source, flt)
			outfile = os.path.join( self.std_ref_dir, sdstdfile )
			print "secondary stds are saved here:%s"%outfile
			np.savetxt(outfile, sdtplphot)


	def match_sources_on_two_images(self, in_imgkey, ref_imgkey, option=None,mode='file', ref_col1=1, ref_col2=2, input_col1=1, input_col2=2, input_order_col=-3, ref_order_col=-3, max_distance=1, matched_xyxy_outfile=None, matched_trans_output=None,  match_result_display=False, show_match_stat=False):
		'''
		match two set of sources detected on input image and reference image

		INPUTS:
			matched_xyxy_outfile:  this is the simplified version of the matched source list with only x,y from two source list


		MIDDLE PRODUCTS:
			match_output: raw matched source list with whole information from input source list and reference source list
			match_trans_output: the matched and transformed source list
		'''
		in_imgkey_s = in_imgkey.split('.')[0]
		ref_imgkey_s = ref_imgkey.split('.')[0]
		ref_list = os.path.join(self.stars_dir, ref_imgkey_s + self.starfile_suffix)
		input_list = os.path.join(self.stars_dir, in_imgkey_s + self.starfile_suffix)
		match_output = os.path.join(self.stars_dir, "match_%s_%s.txt"%(in_imgkey_s, ref_imgkey_s))
		transcoeff_output = os.path.join(self.stars_dir, "match_%s_%s.coeff"%(in_imgkey_s, ref_imgkey_s))

		reflist_data = np.loadtxt(ref_list)
		Nr_ref, Nc_ref = reflist_data.shape

		#check here
		matched_data = self.__fitsh_grmatch(ref_list,input_list, match_output, transcoeff_output, mode =mode)
		wanted_columns = [ref_col1-1, ref_col2-1, input_col1+Nc_ref-1, input_col2+Nc_ref-1]
		matched_xyxy = matched_data[:,wanted_columns]

		if matched_xyxy_outfile is None:
			matched_xyxy_outfile = os.path.join(self.stars_dir, "match_%s_%s_xyxy.txt"%(in_imgkey_s, ref_imgkey_s))
		np.savetxt(matched_xyxy_outfile, matched_xyxy, fmt="%8.2f %8.2f %8.2f %8.2f")

		dxfit_values, dyfit_values = self.__extract_transformation_fitting(transcoeff_output)
		ref_xys = matched_xyxy[:,[0,1]]
		input_xys = matched_xyxy[:,[2,3]]
		transformed_xys = self.__transform_xys(dxfit_values, dyfit_values,ref_xys)

		matched_xyxy_trans = np.hstack((transformed_xys, input_xys))
		if matched_trans_output is None:
			matched_trans_output = os.path.join(self.stars_dir, "match_%s_%s_xyxy_trans.txt"%(in_imgkey_s, ref_imgkey_s))
		np.savetxt(matched_trans_output, matched_xyxy_trans, fmt="%8.2f %8.2f %8.2f %8.2f")

		if match_result_display:
			plt.plot(input_xys[:,0], input_xys[:,1], 'ro', label='input')
			plt.plot(transformed_xys[:,0], transformed_xys[:,1], 'bo', label='reference transformed')
			plt.xlabel('X')
			plt.ylabel('Y')
			plt.legend()
			plt.show()

		if show_match_stat:
			from matplotlib import gridspec
			fig = plt.figure(figsize=(12, 6))
			gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
			ax1 = plt.subplot(gs[0])
			ax1.set_xlabel('$\Delta X$')
			ax1.set_ylabel('$\Delta Y$')
			ax2 = plt.subplot(gs[1])
			ax2.set_xlabel('$\Delta X or \Delta Y$')
			ax2.set_ylabel('N')
			xdiff = input_xys[:,0]-transformed_xys[:,0]
			ydiff = input_xys[:,1]-transformed_xys[:,1]
			ax1.plot(xdiff, ydiff, 'ro')
			ax2.hist(xdiff, bins=20,  alpha=0.5, color='r',  label='X')
			ax2.hist(ydiff, bins=20, alpha = 0.5, color='b', label='Y')
			ax2.legend()
			plt.show()


	def flip_image_data(self, imgkey, invertwhich, which_dir='raw_image', output_image=None):
		'''
		flip array in the left/right direction or up/down direction
		'''
		input_image = self.__get_internal_image(imgkey, which_dir=which_dir)
		hdu_input = fits.open(input_image)
		data_input = hdu_input[0].data

		if invertwhich == 'x':
			data_new = np.fliplr(data_input)
		elif invertwhich == 'y':
			data_new = np.flipud(data_input)
		else:
			raise ValueError('invalid input for invertwhich')

		hdu_output = fits.PrimaryHDU(data=data_new)
		if output_image is None:
			output_image = os.path.join(self.modified_image_dir, imgkey)

		if os.path.exists(output_image):
			os.remove(output_image)

		hdu_output.writeto(output_image)


	def rotate_image_data(self, imgkey, rotdeg=90, which_dir = 'raw_image', output_image = None):
		'''

		INPUTS:
			imgkey:
			rotdeg: 90, 180, 270 clockwise direction

		'''
		input_image = self.__get_internal_image(imgkey, which_dir=which_dir)
		hdu_input = fits.open(input_image)
		data_input = hdu_input[0].data
		rotnum = rotdeg/90
		data_new = np.rot90(data_input, k=rotnum)

		hdu_output = fits.PrimaryHDU(data=data_new)
		if output_image is None:
			output_image = os.path.join(self.modified_image_dir, imgkey)

		if os.path.exists(output_image):
			os.remove(output_image)

		hdu_output.writeto(output_image)


	def interpolate_resample_to_the_template_image(self, imgkey,  matched_source_listfile=None, tpl_imgkey=None, which_dir_input = 'modified_image', verbose=True):
		'''
		Interpolate and resample the input image to the frame of reference image. The matched souce list matched_source_listfile or the reference image tpl_imgkey is required

		INPUTS:
			imgkey:
			matched_source_listfile:
			tpl_imgkey:
			which_dir_input:
		'''
		imgkey_s = imgkey.split('.')[0]

		if matched_source_listfile is None:
			if tpl_imgkey is None:
				raise ValueError("the matched source list file between input image and reference image is not provided...")
			else:
				ref_imgkey_s = tpl_imgkey.split('.')[0]
				matched_source_listfile = os.path.join(self.stars_dir, "match_%s_%s_xyxy.txt"%(imgkey_s, ref_imgkey_s))

		xygrid_parafile = os.path.join(self.parafile_dir, 'xygrid.par')
		xygrid_outfile = os.path.join(self.stars_dir, "%s_xygrid.out"%imgkey_s)
		xygrid_command = "xygrid %s %s %s"%(xygrid_parafile, matched_source_listfile, xygrid_outfile)
		if os.path.exists(xygrid_outfile):
			os.remove(xygrid_outfile)

		if verbose:
			print xygrid_command
		os.system(xygrid_command)

		if not os.path.exists(xygrid_outfile):
			raise IOError("xygrid didn't get the expected output...")
		else:
			resample2_infile = xygrid_outfile

		resample2_parafile = os.path.join(self.parafile_dir, 'resample2.par')
		instrument_parafile = os.path.join(self.parafile_dir, 'instrument.par')
		input_image = self.__get_internal_image(imgkey, which_dir=which_dir)
		interp_image = os.path.join(self.subtraction_dir, 'interp_'+imgkey)
		if os.path.exists(interp_image):
			os.remove(interp_image)

		resample2 = os.path.join(self.base.diapl_dir, 'resample2')
		resample2_command = "%s %s %s %s %s %s"%(resample2, resample2_parafile, instrument_parafile, resample2_infile, input_image, interp_image)
		if verbose:
			print resample2_command
		os.system(resample2_command)


	def get_phot_zpt_single_image(self, image_key, photmethod='apphot', aperture_size_fixed=False, renew_existing_photret=True, img_input_which_dir='raw_image', apphot_centering=True, link_stdref_obs_method='surrounding_search', std_catalog = 'apass', panstarrs_photmethod = 'PSF', ):
		'''
		get photometry zero-point
		'''
		if photmethod == 'apphot':
			photret_file = image_key.split('.')[0]  + self.apphot_ret_file_suffix
                        photret_file_abs = os.path.join(self.aperture_photometry_dir,photret_file)

                        if os.path.exists(photret_file_abs) and (not renew_existing_photret):
                                        print "existing photometry file will be used"
			else:
                    		print image_key
                       		xys = self.__get_xys_on_image(image_key)

				image = self.__get_internal_image(image_key, which_dir=which_dir)
				self.get_apphot_iraf_parameters(image_key, saveparfile=False)
                        	options = self.apphot_iraf_options.copy()
                	        print options

                	        photret = self.__aperture_photometry_apphot_iraf_single_image(image,xys,options, centering=apphot_centering)
                	        #xpos, ypos, flux, bkg,bkg_stdev, sky_counts, target_area,mag,mag_err

                	        photret_file_big = image_key.split('.')[0] +'_whole' + self.apphot_ret_file_suffix
                	        photret_file_big_abs = os.path.join(self.aperture_photometry_dir,photret_file_big)
                	        self.__delete_file_if_exist(photret_file_big_abs)
                	        np.savetxt(photret_file_big_abs,photret,fmt="%6.2f %6.2f %10.3f %8.3f %8.3f %8.3f %8.3f %6.3f %6.3f")
                	        photret_small = photret[:,[0,1,7,8]]

                	        self.__delete_file_if_exist(photret_file_abs)
                	        np.savetxt(photret_file_abs,photret_small,fmt="%6.2f %6.2f %6.3f %6.3f")

			flt = self.photometry_info[image_key]['flt']
			zpt, zpt_err = self.stdcal_single_image(image_key,flt, instmag=0, instmagerr=0, link_stdref_obs_method= link_stdref_obs_method, std_catalog = std_catalog, panstarrs_photmethod =  panstarrs_photmethod,  std_world2img_refimage=None,  stds_shift_x=0, stds_shift_y=0, std_obs_matched = False, verbose = 0, update_stdmag=False)

		elif photmethod == 'psfphot':
			print "not available yet..."
			zpt = 99.99
			zpt_err = 99.99
		else:
			raise ValueError("not recognised input photometry method: %s"%photmethod)


		return zpt, zpt_err


	def image_subtraction(self, which_dir='modified_image', xstamp=8, ystamp=8):
		self.__dict2table()
		self.__find_template_imagekey()

		faildict = {}
		for flt,tplimg in self.templates.items():
			failimgs = self.image_subtraction_flt(flt, tplimg, which_dir=which_dir, xstamp=xstamp, ystamp=ystamp)
			for img in failimgs:
				self.photometry_info[img]['drop'] = 13
			faildict[flt]= failimgs

		return faildict


	def image_subtraction_flt(self, flt, tplimg, updatetable=0, which_dir='modified_image', xstamp=8, ystamp=8):
		'''
		image subtraction process
		'''

		if updatetable:
               		self.__dict2table()

		fltimgs = self.result_table_newsaved[self.result_table_newsaved['flt']==flt]
		faillist = []
		for img in fltimgs['name']:
   			if self.photometry_info[img]['drop'] == 0:
				try:
        				self.match_sources_on_two_images(img, tplimg, match_result_display=0, show_match_stat=0, option=1)
				except:
					faillist.append(img)
					continue
        			if self.photometry_info[img]['bkg']<100:
                			self.add_constant_to_image(img, 100)
        			self.align_and_resample_image_given_matched_source_list(img, tplimg, which_dir=which_dir)
        			self.hotpants_image_subtraction_single(img, tpl_imgkey=tplimg, xstamp=xstamp,ystamp=ystamp, ilthresh=1, tlthresh=1, iuthresh=60000,tuthresh=60000)

		return faillist


	def aperture_photometry_on_subtracted_images_flt(self, flt, xc, yc, centering=False, apsizefixed=False):
		'''
		See self.__aperture_photometry_apphot_iraf_target
		'''
		self.__aperture_photometry_apphot_iraf_target(flt, which_dir = 'subtracted_image', aperture_size_fixed = apsizefixed, imagefwhm=None, centering=centering, x=xc, y=yc)



	def hotpants_image_subtraction_single(self, input_imgkey, tpl_imgkey=None, tplimage=None, fwhm_tplimg = None,  tuthresh=45000, iuthresh=45000, tlthresh=0, ilthresh=0, tgain=None, igain=None, trdnoise=None, irdnoise=None, rkernel=10, fitthresh=20, xregion=1, yregion=1, xstamp=10, ystamp=10, toconvolve=None, normalize='t', ngauss=3, polydegrees = [6,4,2], gaussiansigma=None):
		'''
		INPUTS:
			See self.__hotpants_image_subtraction for details
		'''

		in_imgkey_s = input_imgkey.split('.')[0]
		input_image = os.path.join(self.subtraction_dir, 'interp_%s.fits'%in_imgkey_s)

		if tpl_imgkey is not None:
			tplimage = os.path.join(self.template_dir, 'tplsub_'+tpl_imgkey)
			fwhm_tplimg = self.photometry_info[tpl_imgkey]['fwhm']
		elif tplimage is not None:
			if fwhm_tplimg is None:
				raise ValueError("fwhm of the template is required")
		else:
			raise ValueError("please provide valid template image")

		output_image = os.path.join(self.subtraction_dir, 'sub_'+input_imgkey)
		convolved_image = os.path.join(self.subtraction_dir, 'conv_'+input_imgkey)

		fwhm_inimg  = self.photometry_info[input_imgkey]['fwhm']
		if fwhm_inimg <= 0 or fwhm_inimg ==99.99:
			raise ValueError("non-realistic FWHM value for the input image...")

		telescope_gain = self.base.telescopes_info[self.current_telescope]['gain']
		telescope_rdnoise = self.base.telescopes_info[self.current_telescope]['readnoise']
		if tgain is None:
			tgain = telescope_gain
		if igain is None:
			igain = telescope_gain
		if trdnoise is None:
			trdnoise = telescope_rdnoise
		if irdnoise is None:
			irdnoise = telescope_rdnoise

		if gaussiansigma is None:
			sigdiff = np.sqrt(np.abs(fwhm_tplimg**2 - fwhm_inimg**2))
			sigma1 = np.round(0.5*sigdiff, 1)
			sigma2 = np.round(1*sigdiff, 1)
			sigma3 = np.round(2.0*sigdiff, 1)
			gaussiansigma = [sigma1, sigma2, sigma3]

		#ATTENTIONS!!!
		#hotpants will very likely go wrong if the smallest gaussian width < 0.3
		#it's embarrassing to say that I DON'T know why...
		#check the sigma1, and use np.max([sigma1, 0.3]) as the first gaaussian width
		sigma1 = gaussiansigma[0]
		if sigma1 < 0.3:
			print "!!!! The calculated gaussian width from the difference of the FWHM of two images are too small... "
			sigma1 = 0.3
			sigma2 = 0.6
			sigma3 = 1.2
			gaussiansigma = [sigma1, sigma2, sigma3]

		self.__hotpants_image_subtraction(input_image, tplimage, output_image, convolved_image, tuthresh=tuthresh, iuthresh=iuthresh, tlthresh=tlthresh, ilthresh=ilthresh, tgain=tgain, igain=igain, trdnoise=trdnoise, irdnoise=irdnoise, rkernel=rkernel, fitthresh=fitthresh, xregion=xregion, yregion=yregion, xstamp=xstamp, ystamp=ystamp, toconvolve=toconvolve, normalize=normalize, ngauss=ngauss, polydegrees = polydegrees, gaussiansigma=gaussiansigma)


	def __hotpants_image_subtraction(self, inimage, tplimage, outimage, outconvimg, tuthresh=45000, iuthresh=45000, tlthresh=0, ilthresh=0, tgain=1, igain=1, trdnoise=0, irdnoise=0, rkernel=10, fitthresh=20, xregion=1, yregion=1, xstamp=10, ystamp=10, toconvolve=None, normalize='t', ngauss=3, polydegrees = [6,4,2], gaussiansigma=[0.7, 1.5, 3.0], verbose=1):
		'''
		INPUTS:
			inimage: comparison image to be differenced
			tplimage: template image
			outimage: output difference image
			tuthresh: upper valid data count for template image
			iuthresh: upper valid data count for input image
			tlthresh: lower valid data count for template image
			ilthresh: lower valid data count for input image
			tgain: gain in template
			igain: gain in input image
			trdnoise: readnoise in template
			irdnoise: readnoise in input image
			rkernel: convolution kernel half width
			fitthresh: RMS threshold for good centroid in kernel fit
			xregion: number of image regions in x dimension
			yregion: number of image regions in y dimension
			xstamp: number of each region's stamps in x dimension
			ystamp: number of each region's stamps in y dimension
			toconvolve: force convolution on template or image
			normalize: normalize to template or image or unconvolved, 't', 'i' or 'u'
			oci: output convolved image
			ngauss: number of gaussians which compose kernel
			polydegrees: list of degrees of polynomial associated with gaussian
			gaussiansigma: list of widths of gaussian
		'''

		if len(polydegrees) != len(gaussiansigma):
			raise ValueError("the number of polynomial degrees associated with gaussian should be the same with the number of width of gaussian")

		if len(polydegrees) != ngauss:
			raise ValueError("input polynomial degree list not consistent with ngauss")

		nginput = str(ngauss)
		for polydeg, gsigma in zip(polydegrees, gaussiansigma):
			nginput = nginput + " %s %s"%(polydeg, gsigma)


		if toconvolve is None:
			toconv_input = ''
		elif toconvolve in ['t', 'i']:
			toconv_input = "-c %s"%toconvolve
		else:
			raise ValueError("input for toconvolve is invalid")

		hotpants = os.path.join(self.base.hotpants_dir, 'hotpants')
		hotpants_command = "%s -inim %s -tmplim %s -outim %s -tu %s -iu %s -tl %s -il %s -tg %s -ig %s -tr %s -ir %s -r %s -ft %s -nrx %s -nry %s -nsx %s -nsy %s %s -n %s -oci %s -ng %s"%(hotpants, inimage, tplimage, outimage, tuthresh, iuthresh, tlthresh, ilthresh, tgain, igain,  trdnoise, irdnoise, rkernel, fitthresh, xregion, yregion, xstamp, ystamp, toconv_input, normalize, outconvimg, nginput)

		if verbose:
			print hotpants_command

		os.system(hotpants_command)



	def prepare_template_image(self, tplmethod=1, tplflt=None, tplimgkey=None, tplimg_bricks=None, reference_imgkey=None,  which_dir='raw_image'):
		'''
		This function is used to prepare template image for image subtraction. Different procedures are required in different situations.
			(1) use one single image for template image
			(2) co-adding multiple images, which is ok when the input images are well aligned
			(3) align and resample multiple images and co-add them to build the template image, which is ok when the images are taken under stable conditions consecutively
			(4) convolution with kernel is needed before co-adding in 2)


		INPUTS:
			tplmethod: 1, 2, 3, 4 corresponding to procedures required above
			tplimgkey: the selected template image when tplmethod is 1
			tplimg_bricks: the raw materials to build the template image, list of names of key for the input images for example ['001.fits','002.fits','003.fits']
			reference_img: the reference image for alignment, resampling and even convolution
			which_dir: 'raw_image', 'modified_image', 'subtraction'
		'''

		if tplmethod == 1:
			self.__assign_single_image_as_template(tplimgkey, tplflt, which_dir=which_dir)
		elif tplmethod == 2:
			if which_dir == 'raw_image':
				imgdir = self.raw_image_dir
			elif which_dir == 'modified_image':
				imgdir = self.modified_image_dir
			elif which_dir == 'subtraction':
				imgdir = self.subtraction_dir
			else:
				raise ValueError("input for which_dir is not recognised...")

			imgs_input = open('tplimg_build_bricks_temp.txt', 'awt')
			for img in tplimg_bricks:
				imginput = os.path.join(imgdir, img)
				imgs_input.write(imginput+'\n')

			imgs_input.close()
			if tplflt is None:
				tplimg_output = os.path.join(self.template_dir, 'tplsub_temp.fits')
			else:
				tplimg_output = os.path.join(self.template_dir, 'tplsub_%s.fits'%tplflt)

			imcombine_iraf("@%s"%imgs_input, tplimg_output)
			self.__delete_file_if_exist(imgs_input)
		elif tplmethod == 3:
			if tplflt is None:
				tplimg_output = os.path.join(self.template_dir, 'tplsub_temp.fits')
			else:
				tplimg_output = os.path.join(self.template_dir, 'tplsub_%s.fits'%tplflt)


			self.align_and_stack_multiple_images(tplimg_bricks, reference_imgkey, tplimg_output, input_which_dir=which_dir)

		elif tplmethod == 4:
			print "under construction"


	def __assign_single_image_as_template(self, imgkey, flt, which_dir= 'raw_image'):
		'''
		just copy the given image to the destination directory
		'''
		self.__dict2table()
		info_sf = self.__select_rows_from_table(self.result_table_newsaved,'flt',flt,)
		names = info_sf['name']
		if imgkey not in names:
			raise ValueError("%s not in %s band"%(imgkey,flt))
		tplimage = self.__get_internal_image(imgkey, which_dir=which_dir)
		dstimg = os.path.join(self.template_dir, "tplsub_"+imgkey)
		shutil.copy(tplimage, dstimg)

		for img in names:
			self.photometry_info[img]['template'] = 0
		self.photometry_info[imgkey]['template'] = 1



	def align_and_stack_multiple_images(self, imgkeylist, refimgkey, output_image, input_which_dir='raw_image', renew_matched_sourcelist = False, renew_resampled_img=False):
		'''
		This function is to align and co-add multiple images.

		How?
			(1) align and resample the input images to the frame of reference image, see self.align_and_resample_image_given_matched_source_list
			(2) combine reference image and resampled input images, see imcombine_iraf

		INPUTS:
			imgkeylist: input images list without the reference image
			refimgkey: the reference image which will be one of the images to be combined
		'''

		tplimg_listfile = 'tplimg_build_bricks_temp.txt'
		self.__delete_file_if_exist(tplimg_listfile)

		imgs_input = open(tplimg_listfile, 'awt')

		refimage = self.__get_internal_image(refimgkey, which_dir=input_which_dir)
		imgs_input.write(refimage + '\n')

		ref_imgkey_s  = refimgkey.split('.')[0]

		for imgkey in imgkeylist:
			print imgkey
			if imgkey == refimgkey:
				continue

			in_imgkey_s = imgkey.split('.')[0]

			matched_xyxy_outfile = os.path.join(self.stars_dir, "match_%s_%s_xyxy.txt"%(in_imgkey_s, ref_imgkey_s))
			if (not os.path.exists(matched_xyxy_outfile)) or renew_matched_sourcelist:
				self.match_sources_on_two_images(imgkey, refimgkey, option=None, mode='file', ref_col1=1, ref_col2=2, input_col1=1, input_col2=2, input_order_col=-3, ref_order_col=-3, max_distance=1, matched_xyxy_outfile=None, match_result_display=False, show_match_stat=False)

			resample_outfile = os.path.join(self.subtraction_dir, 'interp_%s.fits'%in_imgkey_s)
			if (not os.path.exists(resample_outfile)) or renew_resampled_img:
				outimg = self.align_and_resample_image_given_matched_source_list(imgkey, refimgkey, which_dir=input_which_dir, verbose=0)
			imgs_input.write(resample_outfile+'\n')

		imgs_input.close()
		input_imglist_file = "@%s"%tplimg_listfile
		imcombine_iraf(input_imglist_file, output_image)
		self.__delete_file_if_exist(tplimg_listfile)



	def align_and_resample_image_given_matched_source_list(self, input_imgkey, ref_imgkey, matched_source_file = None, which_dir='raw_image',  verbose=1):
		'''
		This function is to align the input image with the reference image and resample the input image to the same image coordinate of the reference image. The matched source list should be built in advance and given as input
		INPUTS:
			input_imgkey:
			ref_imgkey:
			matched_source_file: the file containing matched sources with the format of x,y,x,y
		'''

		in_imgkey_s = input_imgkey.split('.')[0]
		ref_imgkey_s = ref_imgkey.split('.')[0]

		if matched_source_file is None:
			matched_source_file = os.path.join(self.stars_dir, "match_%s_%s_xyxy.txt"%(in_imgkey_s, ref_imgkey_s))

		if not os.path.exists(matched_source_file):
			raise IOError("matched source list file %s not exists"%matched_source_file)


		xygrid_para = os.path.join(self.parafile_dir, 'xygrid.par')
		xygrid_outfile = os.path.join(self.stars_dir, 'xygrid_%s_%s.txt'%(in_imgkey_s, ref_imgkey_s))

		xygrid = os.path.join(self.base.diapl_dir, 'xygrid')
		xygrid_command = "%s %s %s %s"%(xygrid, xygrid_para, matched_source_file, xygrid_outfile)

		try:
			if verbose:
				print xygrid_command

			os.system(xygrid_command)
		except:
			print "xygrid failure for image %s ..."%input_imgkey

		resample_para = os.path.join(self.parafile_dir, 'resample2.par')
		instrument_para = os.path.join(self.parafile_dir, 'instrument.par')
		resample_outfile = os.path.join(self.subtraction_dir, 'interp_%s.fits'%in_imgkey_s)

		input_image = self.__get_internal_image(input_imgkey, which_dir=which_dir)
		resample2 = os.path.join(self.base.diapl_dir, 'resample2')
		resample_command  = "%s %s %s %s %s %s"%(resample2, resample_para, instrument_para, xygrid_outfile, input_image, resample_outfile)

		try:
			if verbose:
				print resample_command

			os.system(resample_command)
		except:
			print "resample failure for image %s ..."%input_imgkey

		return resample_outfile



	def add_constant_to_image(self, imgkey, addvalue, input_which_dir='modified_image', output_image=None):
		'''
		add constant value addvalue to the inputimage
		'''
		input_image = self.__get_internal_image(imgkey, which_dir=input_which_dir)
		if output_image is None:
			output_image = input_image

		imarith_iraf(input_image, '+', addvalue, output_image)


	def smarts_nir_image_reduction_sky_subtraction(self, ngroup=2, tstart=None, tend=None, updatetable=1):
		'''
		SMARTS NIR images data reduction: remove the sky by subtracting two different images.
		The latest image to the target image will be used as the sky template


		'''
		if updatetable:
			self.__dict2table()

		for flt in np.unique(self.result_table_newsaved['flt']):
			flttable = self.result_table_newsaved[self.result_table_newsaved['flt']==flt]

			flttable.sort('obstime')
			self.smarts_nir_flt_image_reduction_sky_substraction(flttable, tstart=tstart, tend=tend, ngroup=ngroup)


	def smarts_nir_flt_image_reduction_sky_substraction(self, flttable, ngroup=2, tstart=None, tend=None, verbose=1):
		'''
		INPUTS:
			flttable:
			ngroup:
			tstart:
			tend:
			verbose:
		'''
		if tstart is not None:
			flttable = flttable[flttable['obstime']>tstart]

		if tend is not None:
			flttable = flttable[flttable['obstime']<tend]


		for i,obs in enumerate(flttable):
			img_sci = obs['name']
			img_sci_realname = obs['realimg']

			imod = np.mod(i,ngroup)
			if imod == (ngroup-1):
				isky = i - (ngroup-1)
			else:
				if i == (len(flttable)-1):
					isky = i -1
				else:
					isky = i + 1

			img_sky = flttable[isky]['name']
			img_sky_realname = flttable[isky]['realimg']

			if verbose:
				print "%s - %s"%(img_sci, img_sky)
				print "%s - %s"%(img_sci_realname, img_sky_realname)

			self.__subtract_images(img_sci, img_sky)



	def __subtract_images(self, img1, img2):
		'''
		subtract img1 by img2 and save the output image in self.modified_image_dir

		'''
		img_sci_abs = os.path.join(self.raw_image_dir, img1)
		img_sky_abs = os.path.join(self.raw_image_dir, img2)
		img_ret_abs = os.path.join(self.modified_image_dir, img1)

		imarith_iraf(img_sci_abs, '-', img_sky_abs, img_ret_abs)


	def demonext_coadd_images(self, input_which_dir='raw_image', output_dir=None, updatetable=1):
		'''
		try to co-add multiple images taken on the same night to get high SNR of the combined image
		'''
		if updatetable:
			self.__dict2table()
		dates_all  = [img[1:9]  for img in self.result_table_newsaved['realimg']]
		dates_unique = np.unique(dates_all)

		flts_unique = np.unique(self.result_table_newsaved['flt'])

		for date in dates_unique:
			for flt in flts_unique:
				self.__demonext_coadd_image_night_flt(date, flt, input_which_dir=input_which_dir, output_dir=output_dir)


	def __demonext_coadd_image_night_flt(self, date, flt, input_which_dir='raw_image', output_dir=None, clobber=False):


		if input_which_dir == 'raw_image':
			imgdir = self.raw_image_dir
		elif input_which_dir == 'modified_image':
			imgdir = self.modified_image_dir
		else:
			raise ValueError("the input image dir not recognised")

		input_imgkeylist = [img for img in self.images.keys() if date in self.photometry_info[img]['realimg'] and self.photometry_info[img]['flt']==flt]
		print input_imgkeylist
		refimgkey = input_imgkeylist[0]


                combined_image_file = '%s_%s.fits'%(flt, date)

		if output_dir is None:
			output_dir = os.path.join(self.repository_dir_current_telescope, "coadd_temp")
			print "the co-added images will be stored here: %s"%output_dir
			if os.path.exists(output_dir):
				os.remove(output_dir)

		if not os.path.exists(output_dir):
			os.mkdir(output_dir)

                combined_image = os.path.join(output_dir, combined_image_file)

                if os.path.exists(combined_image) and clobber:
                        os.remove(combined_image)


		self.align_and_stack_multiple_images(input_imgkeylist, refimgkey, combined_image, input_which_dir=input_which_dir, renew_matched_sourcelist = False, renew_resampled_img=False)


#### ===========================

	def one_by_one_check_and_photometry(self, flt, drop_signal=0, psfphot_threshold_default=50):
		'''
		one by one checking the image in filter 'flt' with given drop_signal and decide whether the photometry needed
		'''
		self.__dict2table()
		fltdata = self.result_table_newsaved[self.result_table_newsaved['flt']==flt]
		for img in fltdata['name']:
			if self.photometry_info[img]['drop'] == drop_signal:
				self.display_image_with_ds9(img, which_dir='raw_image')

				drop_default = 0
				drop = raw_input("Enter the drop signal for image %s (default: %s):" %(img,drop_default)) or drop_default
				if drop:
					self.photometry_info[img]['drop'] = int(drop)
				else:
        				bkg = float(raw_input("bkg:"))
        				fwhm = float(raw_input("fwhm:"))
        				self.photometry_info[img]['bkg'] = bkg
        				self.photometry_info[img]['fwhm'] = fwhm
					psfphot_threshold = raw_input('the threshold for peak for sources in psf photometry [default: %s]'%psfphot_threshold_default) or psfphot_threshold_default
					self.psfphot_thresmin = float(psfphot_threshold)
        				self.__psf_photometry_dophot_single(img, output_residuals=1, which_dir='raw_image')
        				self.__psf_photometry_dophot_single_target_ds9_pick(img)
