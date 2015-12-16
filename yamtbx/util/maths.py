"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import numpy

vectors_angle = lambda x, y: abs(numpy.arccos(numpy.dot(x,y)/numpy.linalg.norm(x)/numpy.linalg.norm(y)))
    
def kabsch_superpose(P, Q): # P,Q: vstack'ed matrix
    """
    Usage:
    P = numpy.vstack([a2, b2, c2])
    Q = numpy.vstack([a1, b1, c1])
    m = kabsch_superpose(P, Q)
    newP = numpy.dot(m, P)
    """

    A = numpy.dot(numpy.transpose(P), Q)
    U, s, V = numpy.linalg.svd(A)
    tmp = numpy.identity(3)
    tmp[2,2] = numpy.sign(numpy.linalg.det(A))
    R = numpy.dot(numpy.dot(numpy.transpose(V), tmp), numpy.transpose(U))
    return R
# kabsch_superpose

def rotmat_to_axis_angle(R):
    """
    Reference:
    http://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis-angle
    
    For rotation axis u,
    Ru = u, thus, u is an eigen vector of R.
    for angle, tr(R) = 1 + 2*cos(angle)
    """
    assert R.shape == (3, 3)
    assert numpy.abs(numpy.linalg.det(R) - 1) < 1e-6

    w, v = numpy.linalg.eig(R)
    #print "w,v"
    #print w
    #print v
    rotaxis = None
    for i in xrange(3):
        if abs(numpy.imag(w[i])) < 1e-6 and abs(numpy.real(w[i])-1) < 1e-6:
            #print "i==", i
            rotaxis = map(numpy.real, v[:,i])
            break
    if rotaxis is None: raise "Rotation axis not found"
    
    rotaxis /= numpy.linalg.norm(rotaxis)
    angle = numpy.arccos((R.trace()-1.)/2.)
    
    #print "axis=", rotaxis, "angle= %.3f deg" % numpy.rad2deg(angle)
    return rotaxis, angle
