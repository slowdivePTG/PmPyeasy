import iraf task 

$from pyraf import iraf

#some tasks can import directly from iraf, for example

$from iraf import imcombine
$from iraf import apphot

#some task can't be imported before its parent package is imported 

$from iraf import imred 
$from imred import bias
$from bias import colbias
