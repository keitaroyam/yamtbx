import iotbx.phil
from libtbx import easy_mp
from yamtbx.dataproc import eiger
import h5py
import numpy
import os
master_params_str = """\
bin = 4
 .type = int(value_min=1)
 .help = binning number
dead_area_treatment = 0
 .type = int
 .help = "0: treat as zero, 1: treat as -inf"
nproc = 1
 .type = int(value_min=1)
"""

def software_binning(data, binning, dead_area_treatment):
    assert dead_area_treatment in (0, 1)

    u = l = 0
    b = r = None
    nframes, dimy, dimx = data.shape
    if dimy%binning > 0:
        res = dimy%binning
        u = res//2 # upper
        b = -(res - u) # bottom
    if dimx%binning > 0:
        res = dimx%binning
        l = res//2
        r = -(res - l)

    newdata = data[:, u:b, l:r]
    
    # https://gist.github.com/derricw/95eab740e1b08b78c03f
    newshape = (nframes, dimy//binning, dimx//binning)
    compression_pairs = [(nframes, 1), (dimy//binning, binning), (dimx//binning, binning)]
    flattened = reduce(lambda x,y:x+y, compression_pairs)

    mask = numpy.zeros(newshape[1]*newshape[2]*binning*binning, dtype=numpy.uint16).reshape(newshape[1]*binning,newshape[2]*binning)
    mask[newdata[0,]==2**(newdata.dtype.itemsize*8)-1] = 1
    mask = mask.reshape(newshape[1], binning, newshape[2], -1).sum(3).sum(1).astype(numpy.uint16)

    if dead_area_treatment == 0:
        newdata[newdata==2**(newdata.dtype.itemsize*8)-1] = 0

    newdata = newdata.reshape(flattened)
    for i in xrange(len(newshape)): newdata = newdata.sum(-1*(i+1), dtype=numpy.uint32)

    if dead_area_treatment == 0:
        newdata[:, mask==numpy.max(mask)] = 2**(newdata.dtype.itemsize*8)-1
    elif dead_area_treatment == 1:
        newdata[:, mask>0] = 2**(newdata.dtype.itemsize*8)-1
        
    return newdata, u, l
# software_binning()

def run(params, h5file):
    outfile = os.path.splitext(h5file)[0] + "_bin%d.h5" % params.bin
    h5 = h5py.File(h5file, "r")
    h5out = h5py.File(outfile, "w")

    h5out.create_group("/entry")
    h5out["/entry"].attrs["NX_class"] = "NXentry"

    h5out.create_group("/entry/data")
    h5out["/entry/data"].attrs["NX_class"] = "NXdata"

    def worker(k):
        print "Converting %s" % k
        data = h5["/entry/data"][k][:]
        data, u, l = software_binning(data, params.bin, params.dead_area_treatment)
        dfile = os.path.splitext(h5["/entry/data"].get(k, getlink=True).filename)[0]+"_bin%d.h5"%params.bin
        eiger.create_data_file(os.path.join(os.path.dirname(outfile), dfile),
                               data, (1, data.shape[1], data.shape[2]),
                               h5["/entry/data"][k].attrs["image_nr_low"],
                               h5["/entry/data"][k].attrs["image_nr_high"])
        return u, l
    # worker()

    map_res = easy_mp.pool_map(fixed_func=worker,
                               args=h5["/entry/data"].keys(),
                               processes=params.nproc)

    f_xyconv = lambda x,y: ((x-map_res[0][0])/params.bin, (y-map_res[0][1])/params.bin)

    for k in h5["/entry/data"]:
        dfile = os.path.splitext(h5["/entry/data"].get(k, getlink=True).filename)[0]+"_bin%d.h5"%params.bin
        h5out["/entry/data"][k] = h5py.ExternalLink(dfile, "/entry/data/data")


    print "Converting master.h5"
    if "/entry/sample" in h5: h5out["/entry"].copy(h5["/entry/sample"], "sample")

    h5out.create_group("/entry/instrument")
    h5out["/entry/instrument"].attrs["NX_class"] = "NXinstrument"
    h5out["/entry/instrument"].copy(h5["/entry/instrument/beam"], "beam")

    h5out.create_group("/entry/instrument/detector")
    h5out["/entry/instrument/detector"].attrs["NX_class"] = "NXdetector"
    for k in h5["/entry/instrument/detector"]:
        if k not in ("detectorSpecific", "beam_center_x", "beam_center_y", "x_pixel_size", "y_pixel_size"):
            h5out["/entry/instrument/detector"].copy(h5["/entry/instrument/detector"][k], k)

    beamx, beamy = f_xyconv(h5["/entry/instrument/detector/beam_center_x"].value, h5["/entry/instrument/detector/beam_center_y"].value)
    h5out["/entry/instrument/detector"].create_dataset("beam_center_x", (), dtype=numpy.float32, data=beamx)
    h5out["/entry/instrument/detector/beam_center_x"].attrs["units"] = "pixel"
    h5out["/entry/instrument/detector"].create_dataset("beam_center_y", (), dtype=numpy.float32, data=beamy)
    h5out["/entry/instrument/detector/beam_center_y"].attrs["units"] = "pixel"
    h5out["/entry/instrument/detector"].create_dataset("x_pixel_size", (), dtype=numpy.float32,
                                                       data=h5["/entry/instrument/detector/x_pixel_size"].value*params.bin)
    h5out["/entry/instrument/detector/x_pixel_size"].attrs["units"] = "m"
    h5out["/entry/instrument/detector"].create_dataset("y_pixel_size", (), dtype=numpy.float32, 
                                                       data=h5["/entry/instrument/detector/y_pixel_size"].value*params.bin)
    h5out["/entry/instrument/detector/y_pixel_size"].attrs["units"] = "m"

    h5out.create_group("/entry/instrument/detector/detectorSpecific")
    h5out["/entry/instrument/detector/detectorSpecific"].attrs["NX_class"] = "NXcollection"
    for k in h5["/entry/instrument/detector/detectorSpecific"]:
        if not k.startswith("detectorModule_") and k not in ("flatfield", "pixel_mask", "countrate_correction_count_cutoff", "number_of_excluded_pixels", "x_pixels_in_detector", "y_pixels_in_detector"):
            h5out["/entry/instrument/detector/detectorSpecific"].copy(h5["/entry/instrument/detector/detectorSpecific"][k], k)

    h5out["/entry/instrument/detector/detectorSpecific"].create_dataset("countrate_correction_count_cutoff", (), dtype=numpy.uint32, 
                                                                        data=h5["/entry/instrument/detector/detectorSpecific/countrate_correction_count_cutoff"].value*params.bin)
    #h5out["/entry/instrument/detector/detectorSpecific"].create_dataset("number_of_excluded_pixels", (), dtype=numpy.uint32, 
    #                                                                    data=)
    h5out["/entry/instrument/detector/detectorSpecific"].create_dataset("x_pixels_in_detector", (), dtype=numpy.uint32, 
                                                                        data=h5["/entry/instrument/detector/detectorSpecific/x_pixels_in_detector"].value//params.bin)
    h5out["/entry/instrument/detector/detectorSpecific"].create_dataset("y_pixels_in_detector", (), dtype=numpy.uint32, 
                                                                        data=h5["/entry/instrument/detector/detectorSpecific/y_pixels_in_detector"].value//params.bin)
# run()

def run_from_args(argv):
    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    h5files = cmdline.remaining_args

    for f in h5files:
        run(params, f)
# run_from_args()

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
