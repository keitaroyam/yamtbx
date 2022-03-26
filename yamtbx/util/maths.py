"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import unicode_literals
import numpy
import math
from cctbx.array_family import flex

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
    for i in range(3):
        if abs(numpy.imag(w[i])) < 1e-6 and abs(numpy.real(w[i])-1) < 1e-6:
            #print "i==", i
            rotaxis = list(map(numpy.real, v[:,i]))
            break
    if rotaxis is None: raise "Rotation axis not found"
    
    rotaxis /= numpy.linalg.norm(rotaxis)
    angle = numpy.arccos((R.trace()-1.)/2.)
    
    #print "axis=", rotaxis, "angle= %.3f deg" % numpy.rad2deg(angle)
    return rotaxis, angle

def rodrigues(u, theta):
    # http://mathworld.wolfram.com/RodriguesRotationFormula.html

    w = u / numpy.linalg.norm(u)
    st = numpy.sin(theta)
    ct = numpy.cos(theta)

    rot = numpy.array([[ct + w[0]**2*(1.-ct),       w[0]*w[1]*(1.-ct)-w[2]*st, w[1]*st+w[0]*w[2]*(1.-ct)],
                       [w[2]*st+w[0]*w[1]*(1.-ct),  ct+w[1]**2*(1.-ct),        -w[0]*st+w[1]*w[2]*(1.-ct)],
                       [-w[1]*st+w[0]*w[2]*(1.-ct), w[0]*st+w[1]*w[2]*(1.-ct), ct+w[2]**2*(1.-ct)]])
    return rot
# rodrigues()

def weighted_correlation_coefficient(x, y, w):
    # may be computationally unstable?
    if isinstance(x, numpy.ndarray):
        assert isinstance(y, numpy.ndarray)
        assert isinstance(w, numpy.ndarray)

        m_x = numpy.average(x, weights=w)
        m_y = numpy.average(y, weights=w)
        # 1/sum_w is omitted
        cov = numpy.sum(w*(x-m_x)*(y-m_y))
        var_x = numpy.sum(w*(x-m_x)**2)
        var_y = numpy.sum(w*(y-m_y)**2)

        return cov/numpy.sqrt(var_x)/numpy.sqrt(var_y)
    elif isinstance(x, flex.double):
        assert isinstance(y, flex.double)
        assert isinstance(w, flex.double)
        
        sum_w = flex.sum(w)
        m_x = flex.sum(w*x)/sum_w
        m_y = flex.sum(w*y)/sum_w

        cov = flex.sum(w*(x-m_x)*(y-m_y))
        var_x = flex.sum(w*flex.pow2(x-m_x))
        var_y = flex.sum(w*flex.pow2(y-m_y))

        return cov/math.sqrt(var_x)/math.sqrt(var_y)
    else:
        return weighted_correlation_coefficient(numpy.array(x), numpy.array(y), numpy.array(w))

# weighted_correlation_coefficient()
