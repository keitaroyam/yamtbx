"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import h5py
import numpy
import os
from libtbx import easy_mp
from yamtbx.dataproc.XIO import XIO
from yamtbx.util import read_path_list

def make_geom(f, geom_out):
    h, junk = read_image(f, read_data=False)
    
    s = """\
photon_energy = /LCLS/photon_energy_eV

rigid_group_q0 = q0
rigid_group_collection_connected = q0
rigid_group_collection_independent = q0

q0/res = %(res).4f
q0/min_fs = 0
q0/max_fs = %(fsmax)d
q0/min_ss = 0
q0/max_ss = %(ssmax)d
q0/corner_x = %(cornerx).2f
q0/corner_y = %(cornery).2f
q0/fs = -x
q0/ss = -y
q0/clen = %(clen)f
q0/adu_per_eV = /LCLS/adu_per_eV
q0/max_adu = %(max_adu)d
""" % dict(fsmax=h["size1"]-1, ssmax=h["size2"]-1,
           cornerx=h["orgx"], cornery=h["orgy"],
           clen=h["distance"]/1000., res=1./(h["pixel_size"]*1.e-3), max_adu=65535)

    open(geom_out, "w").write(s)
# make_geom()

def read_image(f, read_data=True):
    h = {}
    data = None

    im = XIO.Image(f)
    h["size1"] = im.header["Width"]
    h["size2"] = im.header["Height"]
    h["pixel_size"] = im.header["PixelX"]
    h["distance"] = im.header["Distance"]
    h["wavelength"] = im.header["Wavelength"]
    h["beamx"] = im.header["BeamX"]
    h["beamy"] = im.header["BeamY"]
    h["orgx"], h["orgy"] = h["beamx"]/h["pixel_size"], h["beamy"]/h["pixel_size"]

    if read_data:
        data = numpy.array(im.getData(), dtype=numpy.uint16).reshape(h["size2"],h["size1"])

    return h, data
# read_cmos_image()

def run_each(f):
    h5out = os.path.basename(f) + ".h5"
    h, data = read_image(f)

    of = h5py.File(h5out, "w")
    grp = of.create_group("LCLS")
    dset = grp.create_dataset("photon_energy_eV", (1,), dtype=numpy.float)
    dset[...] = 12398.4 / h["wavelength"]
    dset = grp.create_dataset("photon_wavelength_A", (1,), dtype=numpy.float)
    dset[...] = h["wavelength"]
    dset = grp.create_dataset("adu_per_eV", (1,), dtype=numpy.float)
    dset[...] = 0.294 * h["wavelength"]/12398.4 ### according to XDS, GAIN = 0.294 @ 1 A
    dset = grp.create_dataset("detector_distance_m", (1,), dtype=numpy.float)
    dset[...] = h["distance"] / 1000.
    dset = grp.create_dataset("beam_xy_px", (2,), dtype=numpy.float)
    dset[...] = (h["orgx"], h["orgy"])
    dset = grp.create_dataset("original_file", (1,), "S%d"%len(f))
    dset[...] = f

    grp = of.create_group("data")
    dset = grp.create_dataset("data", data.shape, dtype=data.dtype, compression="gzip")
    dset[...] = data

    of.close()
    print "Processed:", f
# run_each()

def run(opts, files):
    if len(files) == 1 and files[0].endswith(".lst"):
        files = read_path_list(files[0])

    make_geom(files[0], os.path.basename(files[0])+".geom")
    
    easy_mp.pool_map(fixed_func=run_each,
                     args=files,
                     processes=opts.nproc)
# run()

if __name__ == "__main__":
    import sys
    import optparse

    parser = optparse.OptionParser()
    parser.add_option("-j", "--nproc", dest="nproc", type=int,
                      help="number of parallel runs")

    opts, args = parser.parse_args()

    
    run(opts, args)
