#! /usr/bin/python

#Here I put some catalog search function for convenient catalog search
#First created by Ping Chen on 2016-11-11

#supported methods:
#1. local APASS
#2. IRSA catalogs (API and wget)
#3. General Catalog Access from MAST (API and wget)
#3. astroquery.vizier
#4. astroquery.mast


#asrtoquery.mast.Catalogs The Catalogs class provides access to a subset of the astronomical catalogs stored at MAST. The catalogs currently available through this interface are:
#The Hubble Source Catalog (HSC)
#The GALEX Catalog (V2 and V3)
#The Gaia (DR1 and DR2) and TGAS Catalogs
#The TESS Input Catalog (TIC)
#The TESS Candidate Target List (CTL)
#The Disk Detective Catalog
#PanSTARRS (DR1, DR2)


import os
import numpy as np
import wget

from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
from astroquery.mast import Catalogs

def query_local_APASS(ra,dec,distance,output_file='/home/asassn/ASASSN/photometry_pipeline/txtfiles/apass_example_out.tbl', local_database_dir='/home/asassn/ASASSN/standard_star'):
	'''
	Simle search of local APASS database; this will return entries within square region defined by ra,dec and distance

	INPUTS:
	    	ra: the right ascension of the center of search region in degree unit
	    	dec: the declination of the center of search region in degree unit
	    	distance: the half width of the square in degree unit
		output_file: the filename of output catalog
		local_database_dir: the path to the local APASS data

	PRODUCTS:
		The search results will be store in the output file
	'''

	#define the search squrare
	ramin = ra - distance
	ramax = ra + distance
	decmin = dec - distance
	decmax = dec + distance

	ra_range = np.array([ramin,ramax])
	dec_range = np.array([decmin,decmax])
	print ra_range
	print dec_range

	#first get the search region index
	fileindex = __get_region_index(ra_range,dec_range)
	print "fileindex involved in the search:", fileindex

	#then get the related filename of APASS data
	filenames = __get_region_file(local_database_dir,fileindex)
	print "filenames involved in the search:", filenames

	if os.path.isfile(output_file):
		os.remove(output_file)

	resfid = open(output_file,'a')
	resfid.write('filedname      RAJ2000   e_RAJ2000  DEJ2000   e_DEJ2000     nobs  mobs   Vmag    B-V     Bmag    g_mag    r_mag    i_mag    Vmag   e_B-V   e_B_mag  e_g_mag    e_r_mag  e_i_mag \n')

	#print filenames.keys()
	for i,key in enumerate(filenames.keys()):
	    std_ref =  os.path.join(local_database_dir,filenames[key])

	    print std_ref
	    fid = open(std_ref)
	    line= fid.readline()
	    line= fid.readline()

	    ralow = np.min(ra_range)
	    rahigh = np.max(ra_range)
	    declow = np.min(dec_range)
	    dechigh = np.max(dec_range)

	    print ralow,rahigh
	    print declow,dechigh

	    while line:
	        records = line.split()
	        try:
	            ra = float(records[1])
	            dec = float(records[3])
	        except:
	    	    line = fid.readline()
	    	    continue

	        if ra>ralow and ra<rahigh:
	            if dec>declow and dec<dechigh:
	    	        print ra,dec
	                resfid.write(line)
	        line  = fid.readline()
	    fid.close()
	resfid.close()

	return 0


def __get_region_index(ra,dec):
	'''
	get the corresponding file given ra dec range of the target image
	ra and dec are in units of degree
	'''
	ral = ra[0]
	rah = ra[1]
	decl = dec[0]
	dech = dec[1]
	segments = np.arange(-90,90,5)

	diff1 = decl-segments
	if np.all(diff1):
	    loc = np.where(diff1==diff1[diff1>0][-1])
	    blg1 = segments[loc]
	else:
	    loc = np.where(diff1==0)
	    blg1 = segments[loc]

	diff2 = dech-segments
	if np.all(diff2):
	    loc = np.where(diff2==diff2[diff2>0][-1])
	    blg2 = segments[loc]
	else:
	    loc = np.where(diff2==0)
	    blg2 = segments[loc]

	index = np.array([blg1,blg2])
	print index
	return index


def __get_region_file(stdref_database_dir,fileindex,):
    '''
    prepare
    '''
    std_ref_files = os.listdir(stdref_database_dir)
    filenames = {}
    fileindex = fileindex.ravel()

    for blg in fileindex:
        if blg>0:
	    blg_str = str(blg)
	    print blg_str
	    if blg_str == '5':
		blg_str = '05'
            prefix = 'zp'+ blg_str
	elif blg == 0:
	    prefix = 'zp00'
        else:
	    blg_str = str(blg)
	    if str(blg) == '-5':
            	blg_str = '-05'

            prefix = 'z'+blg_str.replace('-','m')

	print prefix

        filename = [ref for ref in std_ref_files if prefix in ref]
        filenames[prefix] = filename[0]

    return filenames



def query_VO_SCS(RA,Dec,SR,table='fp_psc',out_format='csv', output_file = '/home/asassn/ASASSN/photometry_pipeline/txtfiles/voscs_example_out.tbl'):
	'''
	IRSA catalogs can be searched via the VO Simple Cone Search.
	Please refer to the following link for details:
	https://irsa.ipac.caltech.edu/docs/vo_scs.html

	The SCS API will return all of the available columns.

	The base URL for the Simple Cone Search is:
	http://irsa.ipac.caltech.edu/SCS?

	INPUTS:
		paramters	values					default description
		table:		fp_psc, iraspsc, etc.			fp_psc	IRSA catalog string identifier "catname" (see below).
		RA:		0.0-360.0				NA	Right Ascension in deg (ICRS, but see Note below).
		Dec:		-90.0-90.0				NA	Declination in deg (ICRS, but see Note below).
		SR:		0.0 < SR < 0.75				NA	Cone search radius in deg.
		format:		votable, ipac_table, csv, tsv, fits	csv	Format of the output table (optional).

		Some Popular Catalogs
		catname for "table" parameter	Description
		wise_allwise_p3as_psd		AllWISE Source Catalog
		fp_psc				2MASS Point Source Catalog
		glimpse_s07			GLIMPSE I Spring 07 Catalog (Spitzer)
		cosmos_phot			COSMOS Photometry Catalog
		iraspsc				IRAS Point Source Catalog

	UUTPUT:


	'''

	url_template = "http://irsa.ipac.caltech.edu/SCS?table=which_table&RA=RA_value_degree&DEC=DEC_value_degree&SR=SR_value_degree&format=out_format"
	url = url_template.replace('which_table',table)
	url = url.replace('RA_value_degree', str(RA))
	url = url.replace('DEC_value_degree', str(Dec))
	url = url.replace('SR_value_degree', str(SR))
	url = url.replace('out_format',out_format)

	if os.path.isfile(output_file):
		os.remove(output_file)

	wget.download(url, out=output_file)



def panstarrs_query_Vizier(ra_deg, dec_deg, rad_deg, outfile=None):
	"""
	Query PanSTARRS @ VizieR using astroquery.vizier
	:param ra_deg: RA in degrees
	:param dec_deg: Declination in degrees
	:param rad_deg: field radius in degrees
	:return: astropy.table object
	"""
	Vizier.ROW_LIMIT = -1
	field = coord.SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg), frame='fk5')
	catalog =  Vizier.query_region(field, radius=rad_deg*u.deg, catalog="II/349/ps1")[0]

	if outfile is None:
		outfile = 'PanSTARRS_RA%s_Dec%s_grizy.csv'%(str(ra_deg, dec_deg))

	catalog.write(outfile, format='ascii.csv')


def apass_query_Vizier(ra_deg, dec_deg, rad_deg, outfile=None):
	"""
	Query APASS @ VizieR using astroquery.vizier
	:param ra_deg: RA in degrees
	:param dec_deg: Declination in degrees
	:param rad_deg: field radius in degrees
	:return: astropy.table object
	"""
	Vizier.ROW_LIMIT = -1
	field = coord.SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg), frame='fk5')
	catalog =  Vizier.query_region(field, radius=rad_deg*u.deg, catalog="II/336/apass9")[0]
	if outfile is None:
		outfile = 'APASS_RA%s_Dec%s_griBV.csv'%(str(ra_deg, dec_deg))

	catalog.write(outfile, format='ascii.csv')


def twomass_query_Vizier(ra_deg, dec_deg, rad_deg, outfile=None):
	"""
	Query 2MASS @ VizieR using astroquery.vizier
	:param ra_deg: RA in degrees
	:param dec_deg: Declination in degrees
	:param rad_deg: field radius in degrees
	:return: astropy.table object
	"""
	Vizier.ROW_LIMIT = -1
	field = coord.SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg), frame='fk5')
	catalog =  Vizier.query_region(field, radius=rad_deg*u.deg, catalog="II/246/out")[0]
	if outfile is None:
		outfile = '2MASS_RA%s_Dec%s_JHK.csv'%(str(ra_deg, dec_deg))

	catalog.write(outfile, format='ascii.csv')


def gaia2_query_Vizier(ra_deg, dec_deg, rad_deg, outfile=None):
	"""
	Query Gaia DR2 @ VizieR using astroquery.vizier
	:param ra_deg: RA in degrees
	:param dec_deg: Declination in degrees
	:param rad_deg: field radius in degrees
	:return: astropy.table object
	"""
	Vizier.ROW_LIMIT = -1
	field = coord.SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg), frame='fk5')
	catalog =  Vizier.query_region(field, radius=rad_deg*u.deg, catalog="I/345/gaia2")[0]
	if outfile is None:
		outfile = 'Gaia2_RA%s_Dec%s_table_Vizier.csv'%(str(ra_deg, dec_deg))

	catalog.write(outfile, format='ascii.csv')

def gaia_query_mast_Catalogs(ra_deg, dec_deg, rad_deg, outfile=None, version=2):
	'''
	An optional version parameter allows you to select which version you want, the default is the highest version.
	version=2 equals DR2
	'''
	catalog_data = Catalogs.query_region("%s %s"%(ra_deg, dec_deg), radius=rad_deg, catalog="Gaia", version=version)
	if outfile is None:
		outfile = 'Gaia%s_RA%s_Dec%s_table_MAST.csv'%(version,str(ra_deg), str(dec_deg))

	catalog.write(outfile, format='ascii.csv')





def ps1_query_mast_Catalogs(ra_deg, dec_deg, rad_deg, outfile=None, data_release='dr2', table='mean'):
	'''
	The PanSTARRS Catalog has multiple data releases as well as multiple queryable tables.
	An optional data release parameter allows you to select which data release is desired,
	with the default being the latest version (dr2). The table to query is a required parameter.

	DR 1 Mean, and Stack: 'mean', 'stack'
	DR 2 Mean, Stack, and Detection: 'mean', 'stack', 'detection', 'forced_mean'

	See details below:
	https://catalogs.mast.stsci.edu/docs/panstarrs.html
	'''
	catalog_data = Catalogs.query_region("%s %s"%(ra_deg, dec_deg), radius=rad_deg, catalog="Panstarrs", data_release=data_release, table=table)
	if outfile is None:
		outfile = 'PS1_%s_RA%s_Dec%s_table_MAST.csv'%(data_release,str(ra_deg), str(dec_deg))

	catalog.write(outfile, format='ascii.csv')


def query_General_MAST(RA, Dec, SR, FORMAT=None, catalog=None, filename=None, maxobj=None, magFaintEnd = None, magBrightEnd = None, mindet = None):
	'''
	General Catalog Access : http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?Parameters...

	Required Parameter List
	1 of the following 3 queries - VO ConeSearch, BoxSearch, IDsearch

	RA=ra(deg) &DEC=dec(deg) &SR=search radius(deg)
	BBOX=raMin(deg),decMin(deg),raMax(deg),decMax(deg)
	ID=catID
	Optional Parameters

	FORMAT= VOTABLE(default) | HTML | KML | CSV | TSV | JSON | TEXT(limited set of catalogs)
	CATALOG=GSC23(default) | GSC11 | GSC12 | USNOB | SDSS | FIRST | 2MASS | IRAS | GALEX | GAIA | TGAS | WISE
	| CAOM_OBSCORE | CAOM_OBSPOINTING | PS1V3OBJECTS | PS1V3DETECTIONS
	FILENAME=outputname (directs output to file)
	MAXOBJ=n (limits number of entries returned by brightest magnitude)
	MAGRANGE=bright,faint (limits number of entries returned by limits)
	MINDET=n (minimum numbr of detections PanSTARRS only)
	e.g.

	http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?RA=83.633083&DEC=22.0145&SR=.01&FORMAT=html&CAT=PS1V3OBJECTS&MINDET=25&MAXOBJ=5
	http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?CAT=PS1V3OBJECTS&ID=134410836341049171&FORMAT=HTML
	http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?ra=10&dec=20&sr=0.05&format=html
	http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?bbox=10,20,10.5,21&format=csv&catalog=gsc11
	http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?id=NBQI004317

	More details:
	http://gsss.stsci.edu/Software/WebServices.htm


	NOTE: Currently only Cone search implemented!

	INPUTS:
		RA:  right ascension, degree
		Dec: declination, degree
		SR:  search radius, degree
		FORMAT: default VOTABLE
		catalog: default: GSC23
		filename: output filename
		maxobj: limits number of entries returned by brightest magnitude
		magFaintEnd:
		magBrightEnd:
		mindet: minimum numbr of detections PanSTARRS only
	'''

	#check RA, Dec and SR input
	RA_deg_str = str(RA)
	Dec_deg_str = str(Dec)
	SR_deg_str = str(SR)

	RA_deg_yn  = RA_deg_str.replace('.','',1).isdigit()

	Dec_temp = Dec_deg_str.replace('.','',1)
	if Dec_temp[0] == '-':
		Dec_temp = Dec_temp[1:]

	Dec_deg_yn = Dec_temp.isdigit()
	SR_deg_yn  = SR_deg_str.replace('.','',1).isdigit()

	if (not RA_deg_yn) or (not Dec_deg_yn) or (not SR_deg_yn):
		raise ValueError('RA, Dec and SR must be in the format of degree')


	URL_origin = 'http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?RA=RA_input&DEC=Dec_input&SR=search_radius'

	url = URL_origin.replace('RA_input', RA_deg_str)
	url = url.replace('Dec_input', Dec_deg_str)
	url = url.replace('search_radius', SR_deg_str)

	if FORMAT is not None:
		url = url + '&FORMAT=%s'%FORMAT

	if catalog is not None:
		url = url + '&CATALOG=%s'%catalog

#	if filename is not None:
#		if os.path.isfile(filename):
#			os.remove(filename)
#		url = url + '&FILENAME=%s'%filename

	if maxobj is not None:
		url = url + '&MAXOBJ=%s'%maxobj

	if magFaintEnd is not None and magBrightEnd is not None:
		url = url + '&MAGRANGE=%s,%s'%(magFaintEnd, magBrightEnd)

	if catalog == 'PS1V3OBJECTS' and mindet is not None:
		url = url + '&MINDET=%s'%mindet


	print "The catalog will be downloaded from the address: %s"%url

	if filename is not None:
		if os.path.isfile(filename):
			os.remove(filename)

	wget.download(url, out=filename)


	outlines = open(filename).readlines()[1:]
	#print outlines

	fid = open(filename,'wt')
	for line in outlines:
		fid.write(line)

	fid.close()



if __name__ == '__main__':
	'''
	use the Simple Cone Search (SCS) API to search and download 2MASS catalog
	'''
	import optparse
	parser = optparse.OptionParser()

	def_RA = ''
        parser.add_option('-a', '--ra','--RA', dest="RA", type="string", default=def_RA,   help = "Target Right Ascension in degree [%s]" % def_RA)

        def_DEC = ''
        parser.add_option('-d', '--dec','--DEC', dest="DEC", type="string", default=def_DEC, help = "Target Declination in degrees [%s]" % def_DEC)

        def_SR = 0.15
        parser.add_option('-r', '--radius', dest="search_radius", type="float", default=def_SR,   help = "cone search radius in degree [%s]" % def_SR)

	def_catalog = 'fp_psc'
	parser.add_option('-c', '--catalog', dest='search_catalog', type='string', default=def_catalog, help="the catalog in interest, avaiable ones: apass, fp_psc, glimpse_s07, cosmos_phot, iraspsc [%s]"%def_catalog)

	def_degree = False
        parser.add_option('-D','--degree', dest="radec_degree",action="store_true", default=def_degree, help="whether the RA and Dec inputs are in degree[%s]"%def_degree)

	options, remainder = parser.parse_args()

	RA_str = options.RA
	Dec_str = options.DEC
	SR = options.search_radius
	catalog = options.search_catalog
	radec_degree = options.radec_degree

	if RA_str == '' or Dec_str == '':
		raise IOError( "RA and Dec are required" )

	if not radec_degree:
		from common import radec_format_transformation
		RA, Dec = radec_format_transformation(RA_str, Dec_str)
	else:
		RA = float(RA_str)
		Dec = float(Dec_str)



	valid_catalogs = ['apass', 'fp_psc', 'glimpse_s07', 'cosmos_phot', 'iraspsc', 'GSC23','GSC11','GSC12','USNOB','SDSS','FIRST','2MASS','IRAS','GALEX','GAIA','TGAS','WISE','CAOM_OBSCORE','CAOM_OBSPOINTING','PS1V3OBJECTS', 'PS1V3DETECTIONS']

	IRSA_catalogs = ['fp_psc', 'glimpse_s07', 'cosmos_phot', 'iraspsc']
	MAST_catalogs = ['GSC23','GSC11','GSC12','USNOB','SDSS','FIRST','2MASS','IRAS','GALEX','GAIA','TGAS','WISE','CAOM_OBSCORE','CAOM_OBSPOINTING','PS1V3OBJECTS', 'PS1V3DETECTIONS']

	if catalog not in valid_catalogs:
		raise ValueError("the required catalog not supported yet")

	if catalog == 'apass':
		query_local_APASS(RA, Dec, SR, )
	elif catalog in IRSA_catalogs:
		query_VO_SCS(RA,Dec, SR, table=catalog)
	else:
		FORMAT = 'csv'
		filename = 'VO_MAST_' + catalog + '_example.csv'
		query_General_MAST(RA, Dec, SR, FORMAT=FORMAT, catalog=catalog, filename=filename)
