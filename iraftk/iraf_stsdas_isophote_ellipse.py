#! /anaconda/bin/python

import numpy as np
import pyraf
from pyraf import iraf
from iraf import stsdas
from iraf import analysis
from iraf import isophote
from iraf import ellipse
from iraf import bmodel

def stsdas_analysis_isophote_ellipse(input_image, output, x0, y0, minsma, maxsma, ellip0=0.2, pa0=45, sma0=10,step=0.1, linear=False,  recenter=False, xylearn=False, physical=True, conver=0.05, minit=10, maxit=50, hcenter=False, hellip=False,hpa=False, wander='INDEF', maxgerr=1.0,  olthresh=0.0, integrmode='bi-linear', usclip=3.0, lsclip=3.0, nclip=0, fflag=0.5, mag0=0.0, refer=1.0, zerolevel=0.0, interactive=False):
	'''
	Fit elliptical isophote to galaxy image
	    
	INPUTS:
		x0, y0, ellip0, pa0, sma0: initial ellipse parameters
		x0,y0:center coordinate X,Y; [>=1.0]
		ellip0: ellipticity (defined as 1-b/a), [0.05 to 1.0]
		pa0: position angle [-90 to 90]
		sma0: semi-major axis length, [>=5.0]

		minsma: minimum semi-major axis length to be measured, if set to zero, if central pixel intensity will be measured
		maxsma: maximum semi-major axis length to be measured
	'''
	
	ellipse.unlearn()
	ellipse.x0 = x0
	ellipse.y0 = y0
	ellipse.ellip0 = ellip0
	ellipse.pa0 = pa0  #-90 < PA <= 90
	ellipse.sma0 = sma0
	ellipse.minsma = minsma
	ellipse.maxsma = maxsma
	ellipse.step = step
	ellipse.linear = linear
	ellipse.recenter = recenter
	ellipse.xylearn = xylearn
	ellipse.physical = physical
	ellipse.unlearn() #set the algorithm control parameters for the ellipse task
	ellipse.conver = conver
	ellipse.minit = minit
	ellipse.maxit = maxit
	ellipse.hcenter = hcenter
	ellipse.hellip = hellip
	ellipse.hpa = hpa
	ellipse.wander = wander
	ellipse.maxgerr = maxgerr
	ellipse.olthresh = olthresh #if set to zero, the x0,y0 value found in geompar are used without questioning
	ellipse.unlearn() #set image sampling parameters for the ellipse task
	ellipse.integrmode = integrmode
	ellipse.usclip = usclip
	ellipse.lsclip = lsclip
	ellipse.nclip = nclip
	ellipse.fflag = fflag
	ellipse.unlearn() #pset with paramters that define the magnitude scale
	ellipse.mag0 = mag0
	ellipse.refer = refer
	ellipse.zerolevel = zerolevel
	#ellipse.geompar.unlearn() #set the geometri parameters for the ellipse task
	#ellipse.geompar.x0 = x0
	#ellipse.geompar.y0 = y0
	#ellipse.geompar.ellip0 = ellip0
	#ellipse.geompar.pa0 = pa0  #-90 < PA <= 90
	#ellipse.geompar.sma0 = sma0
	#ellipse.geompar.minsma = minsma
	#ellipse.geompar.maxsma = maxsma
	#ellipse.geompar.step = step
	#ellipse.geompar.linear = linear
	#ellipse.geompar.recenter = recenter
	#ellipse.geompar.xylearn = xylearn
	#ellipse.geompar.physical = physical
	#ellipse.controlpar.unlearn() #set the algorithm control parameters for the ellipse task
	#ellipse.controlpar.conver = conver
	#ellipse.controlpar.minit = minit
	#ellipse.controlpar.maxit = maxit
	#ellipse.controlpar.hcenter = hcenter
	#ellipse.controlpar.hellip = hellip
	#ellipse.controlpar.hpa = hpa
	#ellipse.controlpar.wander = wander
	#ellipse.controlpar.maxgerr = maxgerr
	#ellipse.controlpar.olthresh = olthresh #if set to zero, the x0,y0 value found in geompar are used without questioning
	#ellipse.samplepar.unlearn() #set image sampling parameters for the ellipse task
	#ellipse.samplepar.integrmode = integrmode
	#ellipse.samplepar.usclip = usclip
	#ellipse.samplepar.lsclip = lsclip
	#ellipse.samplepar.nclip = nclip
	#ellipse.samplepar.fflag = fflag
	#ellipse.magpar.unlearn() #pset with paramters that define the magnitude scale
	#ellipse.magpar.mag0 = mag0
	#ellipse.magpar.refer = refer
	#ellipse.magpar.zerolevel = zerolevel
	
	ellipse(input_image, output, x0=x0, y=y0, olthresh=olthresh, maxsma=maxsma, interactive=interactive)


def stsdas_analysis_isophote_bmodel(table, output, parent_image=None, backgr=0):
	'''
	build a model image from the results of isophotal analysis
	'''
	bmodel.unlearn()
	if parent_image is not None:
		bmodel.parent = parent_image
	
	bmodel(table=table, output=output, backgr=backgr)
