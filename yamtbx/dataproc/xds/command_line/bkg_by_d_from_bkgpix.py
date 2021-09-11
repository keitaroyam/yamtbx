#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import sys
import os
import numpy
from yamtbx.dataproc.xds.xparm import XPARM
from yamtbx.dataproc import cbf

def radius_to_d(r_in_px, xparm):
    assert xparm.qx == xparm.qy
    theta = numpy.arctan(r_in_px*xparm.qx/xparm.distance) / 2.
    d = xparm.wavelength / 2. / numpy.sin(theta)
    return d
# radius_to_d()

def xy_to_d(x_in_px, y_in_px, xparm):
    orgx, orgy = xparm.origin
    x, y = x_in_px - orgx, y_in_px - orgy
    r = numpy.sqrt(x**2+y**2)
    return radius_to_d(r, xparm)
# radius_to_d()

def calc_edge_resolution(xparm, nx, ny):
    orgx, orgy = xparm.origin
    min_len = min((orgx, nx-orgx, orgy, ny-orgy))
    return radius_to_d(min_len, xparm)

def run(bkgpix_in, xparm_in, nbins):
    bkgpix, nx, ny = cbf.load_cbf_as_flex(bkgpix_in)
    xparm = XPARM(xparm_in)

    d_min = calc_edge_resolution(xparm, nx, ny)
    print("# edge resolution=", d_min)
    s2_step = (1./d_min**2) / nbins

    bins = [[] for i in range(nbins)]

    for i in range(0, bkgpix.size(), 2):
        val = bkgpix[i]
        if val < 0:
            continue

        x,y = i%nx, int(i/nx)
        d = xy_to_d(x, y, xparm)
        if d < d_min:
            continue

        s2 = 1./d**2
        idx = int(s2/s2_step)
        if idx >= nbins: idx = nbins - 1
        bins[idx].append(val/100.)

    for i in range(nbins):
        dmax, dmin =  1./numpy.sqrt(i*s2_step), 1./numpy.sqrt((i+1)*s2_step)
        print("%7.2f %7.2f %.4f" % (dmax, dmin, sum(bins[i])/len(bins[i])))



if __name__ == "__main__":
    imgin = sys.argv[1]
    xparm = os.path.join(os.path.dirname(imgin), "XPARM.XDS")
    run(imgin, xparm, 100)
