#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

import h5py
import numpy

def run(h5in):
    h5 = h5py.File(h5in, "r")
    prefix = os.path.splitext(os.path.basename(h5in))[0]

    for k in sorted(h5["entry/data"].keys()):
        data = h5["entry/data"][k] # .shape[0]
        for i in range(data.shape[0]):
            binout = "%s_%s_%.6d.bin"%(prefix, k, i)
            img = data[i,]
            img = img.astype(numpy.int32)
            img[img==2**(data.dtype.itemsize*8)-1] = -1

            img.tofile(binout)
            print("Saved:", binout)
            print("Instruction for fit2d")
            print("  Width:", data.shape[2])
            print(" Height:", data.shape[1])
            print("   Type: Integer (%d byte)" % img.dtype.itemsize)
            print(" Signed: %s" % ("No" if img.dtype.kind=="u" else "Yes")) 
            print() 

    return 

    print("""Instruction for R

R
to.read<-file("%(filename)s","rb")
d <- readBin(to.read, integer(), size=%(size)d, signed=%(signed)s, n=%(width)d*%(height)d, endian = "little")
hist(d)
""" % dict(filename=binout, size=arr.dtype.itemsize, signed="TRUE", width=ndimfast, height=ndimmid))
# run()

if __name__ == "__main__":
    import sys
    import os

    h5in = sys.argv[1]
    run(h5in)
