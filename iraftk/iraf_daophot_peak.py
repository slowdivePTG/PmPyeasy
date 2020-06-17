#! /usr/bin/python

import numpy as np
from pyraf import iraf
from iraf import digiphot,apphot,daophot
import os
from astropy.io import fits
from astropy.table import Table


def daophot_peak_iraf(image, photfile, psfimage, peakfile, rejfile, fitrad=3):
	'''
	INPUTS:
		image:
		photfile:
		psfimage:
		peakfile:
		rejfile:
	'''

	daophot.datapars.unlearn()
	daophot.daopars.unlearn()
	daophot.daopars.fitrad = fitrad
	daophot.peak.unlearn()

	daophot.peak(image=image, photfile=photfile,  psfimage=psfimage, peakfile=peakfile, rejfile=rejfile, verify='no', update='no')	



if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()


	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_photfile = 'default'
	parser.add_option('--photfile', dest = 'photfile', type= 'string', default = def_photfile, help='input photometry file')
	
	def_psfimage = 'default'
	parser.add_option('--psfimage', dest = 'psfimage', type= 'string', default = def_psfimage, help='the input psf model image')
	
	def_peakfile = 'default'
	parser.add_option('--peakfile', dest = 'peakfile', type= 'string', default = def_peakfile, help='output photometry file')

	def_rejfile = 'default'
	parser.add_option('--rejfile', dest = 'rejfile', type= 'string', default = def_rejfile, help='output rejection file')


	def_fitrad = 5.0
	parser.add_option('-r','--fitrad', dest='fitrad', type=float, default=def_fitrad, help='fitting radius in scale units; default: %s'%def_fitrad)

	options, remainder = parser.parse_args()
	input_image = options.input_image
	photfile = options.photfile
	psfimage = options.psfimage
	peakfile = options.peakfile
	rejfile = options.rejfile
	fitrad = options.fitrad

	daophot_peak_iraf(input_image, photfile, psfimage, peakfile, rejfile, fitrad=fitrad)
