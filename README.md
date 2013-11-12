planetc
=======

This module contains code for implementing the analytic Mandel & Agol (2002) equations
for planetary transit lightcurves. Uses a C backend with Python/Cython wrappers.

The C backend of planetc was originally based on the one written by Neale Gibson, as part
of his "MyFunctions" package. There has been a fairly substantial reformatting of those 
original C backends for the ones contained in planetc, but the content remains essentially 
the same. On the other hand, the Python routines that perform calls to the C backend routines
have been modified significantly, and the Python/C API wrappers have been completely rewritten 
with the nicer Cython syntax. Tests confirm that the output from planetc exactly matches the
output from the original MyFunctions routines, with identical computation time.

See INSTALL for instructions on how to install the package. 
