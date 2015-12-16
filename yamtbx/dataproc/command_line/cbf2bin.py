#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc import cbf

def run(cbfin, binout):
    # This function only returns signed int.
    arr, ndimfast, ndimmid = cbf.load_minicbf_as_numpy(cbfin, quiet=False)
    arr.tofile(binout)
    print "Saved:", binout
    print
    print "Instruction for fit2d"
    print "  Width:", ndimfast
    print " Height:", ndimmid
    print "   Type: Integer (%d byte)" % arr.dtype.itemsize
    print " Signed: Yes" 
    print 
    print """Instruction for R

R
to.read<-file("%(filename)s","rb")
d <- readBin(to.read, integer(), size=%(size)d, signed=%(signed)s, n=%(width)d*%(height)d, endian = "little")
hist(d)
""" % dict(filename=binout, size=arr.dtype.itemsize, signed="TRUE", width=ndimfast, height=ndimmid)
# run()

if __name__ == "__main__":
    import sys
    import os

    cbfin = sys.argv[1]
    if len(sys.argv) > 2:
        binout = sys.argv[2]
    else:
        binout = os.path.basename(cbfin) + ".bin"

    run(cbfin, binout)
