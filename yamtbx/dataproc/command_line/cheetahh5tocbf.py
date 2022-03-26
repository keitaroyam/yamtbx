#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

from libtbx import easy_mp

from yamtbx.dataproc import cbf

import h5py
import numpy
import os
import sys

nproc = 11

def save_cbf(wavelen, data, cbfout):
    header = """\
# Detector: NOT PILATUS but MPCCD
# Wavelength %f A
""" % wavelen

    height, width = data.shape
    #print "Saving", cbfout, height, width, data
    cbf.save_numpy_data_as_cbf(data.reshape(width*height), width, height, "hdf5_converted", str(cbfout),
                               pilatus_header=header)
# save_cbf()

def convert_single(h5in, root, cbfout):
    f = h5py.File(h5in, "r")

    wavelen = f["%s/photon_wavelength_A"%root][()]
    data = numpy.array(f["%s/data"%root])
    save_cbf(wavelen, data, cbfout)
    f.close()
    print("Processed: %s %s" % (os.path.basename(h5in), root))
# convert_single()

def convert(h5in, par=False):
    f = h5py.File(h5in, "r")
    if "/LCLS/data" in f:
        convert_single(h5in, root="/LCLS", cbfout=os.path.basename(h5in)+".cbf")
    else:
        for tag in f:
            convert_single(h5in, root="/%s"%tag,
                           cbfout="%s_%s.cbf" % (os.path.basename(h5in), tag))
# convert()

def run(h5_files):
    if len(h5_files) == 0: return

    if len(h5_files) > 1:
        easy_mp.pool_map(fixed_func=convert,
                         args=h5_files,
                         processes=nproc)
    else:
        h5in = h5_files[0]
        tags = list(h5py.File(h5in, "r").keys())
        fun = lambda x: convert_single(h5in, root="/%s"%x,
                                       cbfout="%s_%s.cbf" % (os.path.basename(h5in), x))
        for tag in tags: fun(tag)
        return # parallel reading of single file seems buggy..
        easy_mp.pool_map(fixed_func=fun,
                         args=tags,
                         processes=nproc)
# run()

if __name__ == "__main__":
    import sys

    h5_files = [x for x in sys.argv[1:] if x.endswith(".h5")]

    if len(h5_files) == 0:
        print("Usage: %s h5_files" % sys.argv[0])
        quit()

    run(h5_files)
