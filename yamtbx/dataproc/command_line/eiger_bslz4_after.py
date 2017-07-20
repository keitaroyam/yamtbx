"""
(c) RIKEN 2016. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import h5py
from yamtbx.dataproc import eiger
import numpy
import time
import shutil
import os
import traceback
import tempfile

def is_bslz4_applied(h5obj, dspath):
    import bitshuffle.h5
    ds = h5obj[dspath]
    p = h5py.h5d.DatasetID(ds.id.id).get_create_plist()
    n = p.get_nfilters()
    if n < 1: return False
    return p.get_filter(0)[0] == bitshuffle.h5.H5FILTER
# is_bslz4_applied()

def run_safe(infile, check_data=True):
    startt = time.time()
    h5in = h5py.File(infile, "r")

    if is_bslz4_applied(h5in, "/entry/data/data"):
        print "SKIPPING. Already bslz4'd: %s" % infile
        return

    tmpfd, outfile = tempfile.mkstemp(prefix=os.path.basename(infile), dir="/dev/shm")
    os.close(tmpfd)

    # copy and compress data
    data = h5in["/entry/data/data"]
    eiger.create_data_file(outfile, data, data.chunks,
                           h5in["/entry/data/data"].attrs["image_nr_low"],
                           h5in["/entry/data/data"].attrs["image_nr_high"])

    h5in.close()

    eltime = time.time() - startt
    size1 = os.path.getsize(infile) / 1024**2
    size2 = os.path.getsize(outfile) / 1024**2

    if check_data:
        # Check data and overwrite if ok
        h5in = h5py.File(infile, "r")
        h5out = h5py.File(outfile, "r")
        if (h5in["/entry/data/data"][:] == h5out["/entry/data/data"][:]).all():
            print "OK. overwriting with compressed file: %s # %.3f sec %.2f MB -> %.2f MB (%.1f %%)" % (infile, eltime, size1, size2, size2/size1*100.)
            shutil.move(outfile, infile)
        else:
            print "Error! data not match: %s # %.3f sec" % (infile, eltime)
    else:
        print "Overwriting with compressed file: %s # %.3f sec %.2f MB -> %.2f MB (%.1f %%)" % (infile, eltime, size1, size2, size2/size1*100.)
        shutil.move(outfile, infile)

    if os.path.isfile(outfile):
        print " temporary file removed: %s" % outfile
        os.remove(outfile)
    
# run()

def run_from_args(argv):
    for f in argv:
        try:
            run_safe(f)
        except:
            print "Exception with %s" % f
            print traceback.format_exc()


if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
