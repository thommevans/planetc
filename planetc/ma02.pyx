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
#
#####################################################################

import numpy as np
cimport numpy as np

# Read in the external C functions from the ma02_quad_src.c,
# ma02_nonlin_src.c and ma02_utils_src.c source files,
# using the header information in ma02_src.h. Within this
# module only, we will prefix these C functions with c_, so that
# we can define Python functions below with the same names as
# the originals:
cdef extern from "ma02_src.h":
    double c_F_quad "F_quad" ( double NormSep, double RpRs, double gam1, double gam2 )
    double c_F_nonlin "F_nonlin" ( double NormSep, double RpRs, double c1, double c2, double c3, double c4 )

# Define ordinary Python functions that take type-specified
# Python variables as input, the call the C functions, and
# return Python variables as output:

def F_quad( np.ndarray[ np.double_t, ndim=1 ] NormSep, double RpRs, double gam1, double gam2 ):
    """
    """

    cdef unsigned int i, n
    n = NormSep.shape[0]
    cdef np.ndarray[ np.double_t, ndim=1 ] output = np.zeros( n, dtype=np.double )

    for i in xrange( n ):
        output[i] = c_F_quad( NormSep[i], RpRs, gam1, gam2 )
        
    return output

def F_nonlin( np.ndarray[ np.double_t, ndim=1 ] NormSep, double RpRs, double c1, double c2, double c3, double c4 ):
    """
    """

    cdef unsigned int i, n
    n = NormSep.shape[0]
    cdef np.ndarray[ np.double_t, ndim=1 ] output = np.zeros( n, dtype=np.double )

    for i in xrange( n ):
        output[i] = c_F_nonlin( NormSep[i], RpRs, c1, c2, c3, c4 )
        
    return output
