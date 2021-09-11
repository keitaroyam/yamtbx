"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
import h5py
import numpy
from iotbx.detectors import ImageFactory
import os
CBF_BYTE_OFFSET= 0x0070

def run(img_in):
    im = ImageFactory(img_in)
    im.read()
    print(dir(im))
    print(im.size2, im.size1)
    data = numpy.array(im.linearintdata, dtype=numpy.uint16).reshape(im.size2, im.size1)
    print(data, data.dtype)
    
    prefix = os.path.basename(img_in)
    of = h5py.File("%s_byteoffset.h5"%prefix, "w")
    grp = of.create_group("data")
    dset = grp.create_dataset(prefix, data.shape, dtype=data.dtype, compression=CBF_BYTE_OFFSET)
    dset[...] = data
    of.close()

if __name__ == "__main__":
    import sys
    img_in = sys.argv[1]
    run(img_in)
