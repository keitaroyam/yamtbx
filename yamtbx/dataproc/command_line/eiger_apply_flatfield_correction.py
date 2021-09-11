from __future__ import print_function
from __future__ import unicode_literals
from builtins import range
import h5py
import bitshuffle.h5
from yamtbx.dataproc import cbf
from yamtbx.dataproc import eiger
import numpy
import os

def run(infile):
    h5 = h5py.File(infile, "r")

    if h5["/entry/instrument/detector/flatfield_correction_applied"].value:
        print("Correction was already applied.")
        return 

    ff = h5["/entry/instrument/detector/detectorSpecific/flatfield"][:]

    outdir = "flatfield_applied_%s" % os.path.splitext(os.path.basename(infile))[0]
    if not os.path.exists(outdir): os.mkdir(outdir)

    cbf.save_numpy_data_as_cbf((ff*1000).astype(numpy.uint32).flatten(), size1=ff.shape[1], size2=ff.shape[0], title="flatfield_x1000",
                               cbfout=os.path.join(outdir, "flatfield.cbf"))

    for k in sorted(h5["/entry/data"].keys()):
        outfile = os.path.join(outdir, h5["/entry/data"].get(k, getlink=True).filename)
        ds = h5["/entry/data"][k]
        
        h5o = h5py.File(outfile, "w")
        h5o.create_group("/entry")
        h5o["/entry"].attrs["NX_class"] = "NXentry"
        h5o.create_group("/entry/data")
        h5o["/entry/data"].attrs["NX_class"] = "NXdata"
        dataset = h5o.create_dataset("/entry/data/data", ds.shape,
                                       compression=bitshuffle.h5.H5FILTER,
                                       compression_opts=(0, bitshuffle.h5.H5_COMPRESS_LZ4),
                                       chunks=ds.chunks, dtype=ds.dtype)
        dataset.attrs["image_nr_low"] = ds.attrs["image_nr_low"]
        dataset.attrs["image_nr_high"] = ds.attrs["image_nr_high"]

        # Don't want to load all data to memory at the same time
        for i in range(ds.shape[0]):
            print("processing %6d/%d in %s" % (i+1, ds.shape[0], k))
            dataset[i,] = ds[i,] * ff + .5

# run()
    

if __name__ == "__main__":
    import sys
    masterh5 = sys.argv[1]
    run(masterh5)
