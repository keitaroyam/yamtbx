"""
(c) RIKEN 2016. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
import h5py
import numpy
import bitshuffle.h5
import time
import shutil
import os
import traceback
import tempfile
import glob
import math
from yamtbx.dataproc import eiger

def run(infile, nframes, tmpdir="/dev/shm"):
    wdir = tempfile.mkdtemp(prefix="h5split", dir=tmpdir)
    orgdir = os.path.normpath(os.path.dirname(infile))
    print("Workdir: %s" % wdir)
    
    # Copy original master file to wdir
    infile_tmp = os.path.join(wdir, os.path.basename(infile))
    shutil.copyfile(infile, infile_tmp)

    h5in = h5py.File(infile_tmp, "a")
    h5org = h5py.File(infile, "r")

    datasets = []
    n_all = max([h5org["/entry/data"][k].attrs["image_nr_high"] for k in list(h5org["/entry/data"].keys())])
    lookup = [0 for x in range(n_all)]
    org_files = [infile]

    print("Reading original data")
    for i, k in enumerate(sorted(h5org["/entry/data"].keys())):
        print(" %s %s" % (k, h5org["/entry/data"][k].shape))
        datasets.append(h5org["/entry/data"][k])
        low, high = h5org["/entry/data"][k].attrs["image_nr_low"], h5org["/entry/data"][k].attrs["image_nr_high"]
        for j in range(low, high+1): lookup[j-1] = i

        del h5in["/entry/data"][k]
        org_files.append(os.path.join(orgdir, h5org["/entry/data"].get(k, getlink=True).filename))

    # Write data
    cur_idx = 0
    for i in range(int(math.ceil(n_all/float(nframes)))):
        outname = "data_%.6d" % (i+1)
        print("preparing", outname)
        newlow, newhigh = i*nframes+1, min((i+1)*nframes, n_all)
        if lookup[newlow-1] == lookup[newhigh-1]:
            lidx = len([x for x in lookup[:newlow-1] if x==lookup[newlow-1]])
            ridx = lidx + (newhigh - newlow + 1)
            data = datasets[lookup[newlow-1]][lidx:ridx]
            print(" data_%.6d [%6d, %6d)" % (lookup[newlow-1]+1, lidx, ridx))
        else:
            data = None
            for j in range(lookup[newlow-1], lookup[newhigh-1]+1):
                if j == lookup[newlow-1]:
                    lidx = len([x for x in lookup[:newlow] if x==j]) - 1
                    ridx = None # till end
                elif j == lookup[newhigh-1]:
                    lidx = 0
                    ridx = len([x for x in lookup[:newhigh] if x==j])
                else:
                    lidx = 0
                    ridx = None # till end

                print(" data_%.6d [%6s, %6s)" % (j+1, lidx, ridx))
                if data is None: data = datasets[j][lidx:ridx]
                else: data = numpy.concatenate((data, datasets[j][lidx:ridx]))
        
        eiger.create_data_file(os.path.join(wdir, outname+".h5"), data, datasets[0].chunks, newlow, newhigh)
        h5in["/entry/data/%s"%outname] = h5py.ExternalLink(outname+".h5", "/entry/data/data")
        print(" wrote %s %s" % (outname+".h5", data.shape))

    h5in.close()

    bdir = os.path.join(orgdir, "split_org_%s"%time.strftime("%y%m%d-%H%M%S"))
    os.mkdir(bdir)
    print("Moving old files")
    for f in org_files:
        print(" %s to %s" % (f, bdir))
        shutil.move(f, bdir)

    print("Moving new files")
    for f in glob.glob(os.path.join(wdir, "*")):
        print(" %s to %s" % (f, orgdir))
        shutil.move(f, orgdir)

    os.rmdir(wdir)

# run()

if __name__ == "__main__":
    import sys
    run_safe(sys.argv[1], int(sys.argv[2]))
