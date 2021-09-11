#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
from yamtbx.dataproc import cbf
import h5py
import numpy
import os

def make_dummy_pilatus_header(h5):
    safe_str = lambda x,y,d: x[y].value if y in x else d
    safe_val = lambda x,y,d: x[y][0] if y in x else d

    try:
        det = h5["entry"]["instrument"]["detector"]
    except KeyError:
        return ""

    return """\
# Detector: NOT PILATUS but %(detector)s
# %(date)s
# Pixel_size %(pixel_size_x)e m x %(pixel_size_y)e m
# %(sensor_material)s sensor, thickness %(sensor_thickness)e m^M
# Exposure_time %(exp_time)f s
# Wavelength %(wavelength)f A
# Detector_distance %(distance)f m
# Beam_xy (%(orgx).2f, %(orgy).2f) pixels
# Angle_increment %(osc_width)f deg.
""" % dict(detector= safe_str(det, "description", ""),
           date= safe_str(det["detectorSpecific"], "data_collection_date", ""),
           pixel_size_x= safe_val(det, "x_pixel_size", 0),
           pixel_size_y= safe_val(det, "y_pixel_size", 0),
           sensor_material= safe_str(det, "sensor_material", ""),
           sensor_thickness= safe_val(det, "sensor_thickness", 0),
           exp_time= safe_val(det, "count_time", 0),
           wavelength= safe_val(h5["entry"]["instrument"]["beam"], "wavelength", 0) if "beam" in h5["entry"]["instrument"] else 0,
           distance= safe_val(det, "detector_distance", 0),
           orgx= safe_val(det, "beam_center_x", 0),
           orgy= safe_val(det, "beam_center_y", 0),
           osc_width=-1, # Should be found in f["entry"]["sample"]["rotation_angle_step"]
           )

# make_dummy_pilatus_header()

def get_mask_info(val):
    """
https://www.dectris.com/nexus.html#main_head_navigation

Contains a bit field for each pixel to signal dead, blind or high or otherwise unwanted or undesirable pixels.
The bits have the following meaning:

0 gap (pixel with no sensor)
1 dead
2 under responding
3 over responding
4 noisy
5-31 -undefined-	
"""
    ret = []
    if val&1:
        ret.append("gap (0)")
    if val&2:
        ret.append("dead (1)")
    if val&4:
        ret.append("under responding (2)")
    if val&8:
        ret.append("over responding (3)")
    if val&16:
        ret.append("noisy (4)")

    return ",".join(ret)
# get_mask_info()

def run(h5in, cbf_prefix):
    f = h5py.File(h5in, "r")
    if "instrument" not in f["entry"]:
        print("Error: This is not master h5 file.")
        return

    dname = os.path.dirname(cbf_prefix)
    if dname != "" and not os.path.exists(dname):
        os.makedirs(dname)
        print("dir created:", dname)

    # Analyze pixel_mask
    pixel_mask = numpy.array(f["entry"]["instrument"]["detector"]["detectorSpecific"]["pixel_mask"])
    print("No. of unuseful pixels:")
    for val in set(pixel_mask.flatten()):
        if val==0: continue
        print("", get_mask_info(val), (pixel_mask==val).sum())
    print()

    # Extract and write data
    data = [x for x in f["entry"] if x.startswith("data_")]
    count = 0
    for key in sorted(data):
        print("Extracting", key)
        im = f["entry"][key]
        print(" shape=", im.shape)
        print(" dtype=", im.dtype)
        nframes, height, width = im.shape
        for i in range(nframes):
            count += 1
            cbfout = "%s_%.6d.cbf" % (cbf_prefix, count)
            data = im[i,].astype(numpy.int32)
            data[pixel_mask>0] = -1
            cbf.save_numpy_data_as_cbf(data.reshape(width*height), width, height, "hdf5_converted", cbfout, 
                                       pilatus_header=make_dummy_pilatus_header(f))
           
            print(" ", cbfout)
        print()
# run()

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: %s master-h5 prefix" % os.path.basename(sys.argv[0]))
        print()
        print("for example: %s series_10_master.h5 /tmp/series10" % os.path.basename(sys.argv[0]))
        print(" then writes /tmp/series10_000001.cbf, ....")
        quit()

    h5in = sys.argv[1]

    if len(sys.argv) > 2:
        cbf_prefix = sys.argv[2]
    else:
        p = os.path.basename(h5in).replace("_master.h5","")
        cbf_prefix = "cbf_%s/%s" % (p,p) 

    run(h5in, cbf_prefix)
