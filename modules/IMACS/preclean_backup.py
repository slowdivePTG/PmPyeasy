#PmPyeasy.iraf.ccdproc(PmPyeasy.WorkDir+'*_SQ???.fits', output="", ccdtype="",  fixpix='no', oversca='no', trim='yes', zerocor='no', darkcor='no', flatcor='no', trimsec='[100:200,100:200]', zero='bias', flat='nflat')
PmPyeasy.iraf.ccdproc.unlearn()
PmPyeasy.iraf.ccdproc(PmPyeasy.WorkDir+'*_SQ???.fits', output="", ccdtype="",  fixpix='no', oversca='yes', trim='no', zerocor='no', darkcor='no', flatcor='no', biassec="image",  readaxis="line", trimsec='[300:800,*]', zero='bias', flat='nflat')

PmPyeasy.iraf.hedit.unlearn()
#PmPyeasy.iraf.hedit(PmPyeasy.WorkDir+'*_SQ???.fits', 'EPOCH', "(EQUINOX)", add='yes', addonly='yes', verify='no')
PmPyeasy.iraf.hedit(PmPyeasy.WorkDir+'*_SQ???.fits', 'ccdsec,datasec', delete='yes', verify='no', update='yes')

#PmPyeasy.iraf.hedit(PmPyeasy.WorkDir+'*_SQ???.fits', "DATE-OBS", '(@"DATE-OBS"//"T"//@"TIME-OBS")', add='yes', addonly='yes', verify='no',show='no', update='yes')

PmPyeasy.iraf.hedit.unlearn()
PmPyeasy.iraf.hedit(PmPyeasy.WorkDir+'*_SQ???.fits', 'OBSERVAT', '(@"SITENAME")', add='yes', addonly='yes', verify='no')

PmPyeasy.iraf.setjd.unlearn()
PmPyeasy.iraf.setjd(PmPyeasy.WorkDir+'*_SQ???.fits', time="TIME-OBS", epoch="EQUINOX")