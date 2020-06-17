#! /usr/bin/python

import numpy as np
from pyraf import iraf
from iraf import digiphot,apphot,daophot
import os
from astropy.io import fits
from astropy.table import Table


def daophot_seepsf_iraf(psfimage, outputimg):
	'''
	INPUTS:
		psfimage:
		outputimg:
	'''

	daophot.seepsf(psfimage, outputimg)



if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()


	def_psfimage = ''
	parser.add_option('-i','--psfimage', dest = 'psfimage', type= 'string', default = def_psfimage, help='input psf image')

	def_outputimg = ''
	parser.add_option('-o','--outputimg', dest = 'outputimg', type= 'string', default = def_outputimg, help='output psf image')
	
	
	options, remainder = parser.parse_args()
	psfimage = options.psfimage
	outputimg = options.outputimg

	daophot_seepsf_iraf(psfimage, outputimg)
