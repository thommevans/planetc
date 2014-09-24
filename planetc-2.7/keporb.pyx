# The following two lines give more efficient C-like handling of
# array indices, and must come before any lines of code (but #'s
# and whitespace are OK):

#cython: boundscheck=False
#cython: wraparound=False

#####################################################################
#
#  COMMENTS:
#
#  2012nov TME 
#
#  This Cython script takes the Python-esque routines below, which
#  themselves plug into the external C module keporb_src.c, and
#  creates a fully C-translated script keporb.c, which can then be
#  compiled into a Python-readable *.so binary using a setup.py
#  script.
#
#  The keporb_src.c has been adapted directly from Neale Gibson's
#  PlanetOrbit_functions.c module, with some small changes (mostly
#  cosmetic). Time trials comparing the speeds of routines in Neale's
#  PlanetOrbit.so module and the current keporb.so module confirm
#  that the performances are near identical. The Cython syntax 
#  required for this module, however, is **far** simpler than the
#  Python/C API syntax that Neale's PlanetOrbit.c wrapper uses,
#  certainly for a user more familiar with Python than C.
#
#####################################################################

import numpy as np
cimport numpy as np

# Read in the external C functions from the keporb_src.c file,
# using the header information in keporb_src.h. Within this
# module only, we will prefix these C functions with c_, so that
# we can define Python functions below with the same names as
# the originals:
cdef extern from "keporb_src.h":
    
    double c_NormSep "NormSep" ( double MeanAnom, double aRs, double ecc, double omega, double incl )
    double c_Xcoord "Xcoord" ( double MeanAnom, double aRs, double ecc, double omega )
    double c_Ycoord "Ycoord" ( double MeanAnom, double aRs, double ecc, double omega, double incl )
    double c_Zcoord "Zcoord" ( double MeanAnom, double aRs, double ecc, double omega, double incl )
    double c_calc_EccAnom "calc_EccAnom" ( double MeanAnom, double ecc )
    double c_calc_TrueAnom "calc_TrueAnom" ( double EccAnom, double ecc )
    double c_calc_cos_TrueAnom "calc_cos_TrueAnom" ( double EccAnom, double ecc )
    double c_calc_sin_TrueAnom "calc_sin_TrueAnom" ( double EccAnom, double ecc )    


# Define ordinary Python functions that take type-specified
# Python variables as input, the call the C functions, and
# return Python variables as output:

def NormSep( np.ndarray[ np.double_t, ndim=1 ] MeanAnom, double aRs, double ecc, double omega, double incl ):
    """
    Calculates the on-sky separation between the centres 
    of the primary and secondary discs, in units of the
    primary radius.
    """

    cdef unsigned int i, n
    n = MeanAnom.shape[0]
    cdef np.ndarray[ np.double_t, ndim=1 ] output = np.zeros( n, dtype=np.double )

    for i in xrange( n ):
        output[i] = c_NormSep( MeanAnom[i], aRs, ecc, omega, incl )
        
    return output


def Xcoord( np.ndarray[ np.double_t, ndim=1 ] MeanAnom, double aRs, double ecc, double omega ):
    """
    Calculates the X coordinate of the secondary on the
    sky plane, in the XYZ coordinate system where:
      - The origin is located at the system barycentre.
      - The X axis is the 'horizontal' on-sky axis, running
        along the line of nodes.
      - The Y axis is the other on-sky axis.
      - The Z axis is perpendicular to the sky plane,
        increasing in the direction away from the observer.
    """

    cdef unsigned int i, n
    n = MeanAnom.shape[0]
    cdef np.ndarray[ np.double_t, ndim=1 ] output = np.zeros( n, dtype=np.double )

    for i in xrange( n ):
        output[i] = c_Xcoord( MeanAnom[i], aRs, ecc, omega )
        
    return output


def Ycoord( np.ndarray[ np.double_t, ndim=1 ] MeanAnom, double aRs, double ecc, double omega, double incl ):
    """
    Calculates the Y coordinate of the secondary on the
    sky plane, in the XYZ coordinate system where:
      - The origin is located at the system barycentre.
      - The X axis is the 'horizontal' on-sky axis, running
        along the line of nodes.
      - The Y axis is the other on-sky axis.
      - The Z axis is perpendicular to the sky plane,
        increasing in the direction away from the observer.
    """

    cdef unsigned int i, n
    n = MeanAnom.shape[0]
    cdef np.ndarray[ np.double_t, ndim=1 ] output = np.zeros( n, dtype=np.double )

    for i in xrange( n ):
        output[i] = c_Ycoord( MeanAnom[i], aRs, ecc, omega, incl )
        
    return output


def Zcoord( np.ndarray[ np.double_t, ndim=1 ] MeanAnom, double aRs, double ecc, double omega, double incl ):
    """
    Calculates the Z coordinate of the secondary on the
    sky plane, in the XYZ coordinate system where:
      - The origin is located at the system barycentre.
      - The X axis is the 'horizontal' on-sky axis, running
        along the line of nodes.
      - The Y axis is the other on-sky axis.
      - The Z axis is perpendicular to the sky plane,
        increasing in the direction away from the observer.
    """

    cdef unsigned int i, n
    n = MeanAnom.shape[0]
    cdef np.ndarray[ np.double_t, ndim=1 ] output = np.zeros( n, dtype=np.double )

    for i in xrange( n ):
        output[i] = c_Zcoord( MeanAnom[i], aRs, ecc, omega, incl )
        
    return output


def calc_EccAnom( np.ndarray[ np.double_t, ndim=1 ] MeanAnom, double ecc ):
    """
    Calculates the eccentric anomaly by solving the Kepler
    equation numerically with a Newton-Raphson iteration
    method, as outlined in Solar System Dynamics by Murray
    & Dermott (1999).
    """

    cdef unsigned int i, n
    n = MeanAnom.shape[0]
    cdef np.ndarray[ np.double_t, ndim=1 ] output = np.zeros( n, dtype=np.double )

    for i in xrange( n ):
        output[i] = c_calc_EccAnom( MeanAnom[i], ecc )
        
    return output


def calc_TrueAnom( np.ndarray[ np.double_t, ndim=1 ] EccAnom, double ecc ):
    """
    Calculates true anomaly.
    """

    cdef unsigned int i, n
    n = EccAnom.shape[0]
    cdef np.ndarray[ np.double_t, ndim=1 ] output = np.zeros( n, dtype=np.double )

    for i in xrange( n ):
        output[i] = c_calc_TrueAnom( EccAnom[i], ecc )
        
    return output


def calc_cos_TrueAnom( np.ndarray[ np.double_t, ndim=1 ] EccAnom, double ecc ):
    """
    Calculates the cosine of the true anomaly.
    This is just a wrapper for a C backend that
    calculates the cosine of the true anomaly
    as part of the NormSep backend, but is
    included here for completeness.
    """

    cdef unsigned int i, n
    n = EccAnom.shape[0]
    cdef np.ndarray[ np.double_t, ndim=1 ] output = np.zeros( n, dtype=np.double )

    for i in xrange( n ):
        output[i] = c_calc_cos_TrueAnom( EccAnom[i], ecc )
        
    return output


def calc_sin_TrueAnom( np.ndarray[ np.double_t, ndim=1 ] EccAnom, double ecc ):
    """
    Calculates the sine of the true anomaly.
    This is just a wrapper for a C backend that
    calculates the sine of the true anomaly as
    part of the NormSep backend, but is included
    here for completeness.
    """

    cdef unsigned int i, n
    n = EccAnom.shape[0]
    cdef np.ndarray[ np.double_t, ndim=1 ] output = np.zeros( n, dtype=np.double )

    for i in xrange( n ):
        output[i] = c_calc_sin_TrueAnom( EccAnom[i], ecc )
        
    return output

