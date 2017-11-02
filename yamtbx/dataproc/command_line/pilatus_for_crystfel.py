#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import iotbx.phil
from libtbx import easy_mp
from cctbx.array_family import flex
from cctbx import miller

from yamtbx.dataproc import cbf
from yamtbx.dataproc import XIO

import h5py
import numpy
import os
import sys

master_params_str="""\
geom_out = pilatus.geom
 .type = path
beam_out = pilatus.beam
 .type = path
byteoffset = True
 .type = bool
decompose = True
 .type = bool
nproc = 1
 .type = int
"""

def make_geom(header, geom_out):
    s = """\
0/res = %(res).1f ; 1m /x micron
0/min_fs = 0
0/max_fs = %(fsmax)d
0/min_ss = 0
0/max_ss = %(ssmax)d
0/corner_x = %(cornerx).2f ; units are pixels
0/corner_y = %(cornery).2f
0/fs = -x
0/ss = -y
0/clen = %(clen)f
0/adu_per_eV = /LCLS/adu_per_eV
0/max_adu = %(max_adu)d
""" % dict(fsmax=header["Width"]-1, ssmax=header["Height"]-1,
           cornerx=header["BeamX"]/header["PixelX"], cornery=header["BeamY"]/header["PixelY"],
           clen=header["Distance"]/1000., res=1/(header["PixelX"]*1.e-3), max_adu=1048576)
    open(geom_out, "w").write(s)
# make_geom()

def make_geom_decomposed(header, shape, geom_out):
    corner = header["BeamX"]/header["PixelX"], header["BeamY"]/header["PixelY"]

    s = """\
clen = %(clen)f
res = %(res).1f ; 1m /x micron
adu_per_eV = /LCLS/adu_per_eV
max_adu = %(max_adu)d
;corner_x = %(cornerx).2f ; units are pixels
;corner_y = %(cornery).2f

""" % dict(clen=header["Distance"]/1000., res=1/(header["PixelX"]*1.e-3), max_adu=1048576,
           cornerx=corner[0], cornery=corner[1])

    panel = 0
    # PILATSU 6M
    if shape == (195*12, 5*487):
        for j in xrange(12):
            for i in xrange(5):
                panel += 1
                s += """\
%(panel)d/min_fs = %(min_fs)d
%(panel)d/max_fs = %(max_fs)d
%(panel)d/min_ss = %(min_ss)d
%(panel)d/max_ss = %(max_ss)d
%(panel)d/fs = x
%(panel)d/ss = y
%(panel)d/corner_x = %(cornerx).2f ; units are pixels
%(panel)d/corner_y = %(cornery).2f

""" % dict(panel=panel,
           min_ss=j*195, max_ss=(j+1)*195-1,
           min_fs=i*487, max_fs=(i+1)*487-1, # why need to swap?
           cornerx=-corner[0]+494*i, cornery=-corner[1]+212*j
           )
    else:
        raise "Never reaches here."
    open(geom_out, "w").write(s)
# make_geom_decomposed()

def make_geom_decomposed_for_raw(header, shape, geom_out):
    # For raw pilatus image (gaps with -1s).
    # This function is of no use, because CrystFEL complains if there is gap.. blugh.
    s = """\
clen = %(clen)f
res = %(res).1f ; 1m /x micron
adu_per_eV = /LCLS/adu_per_eV
max_adu = %(max_adu)d

""" % dict(clen=header["Distance"]/1000., res=1/(header["PixelX"]*1.e-3), max_adu=1048576)

    panel = 0
    corner = header["BeamX"]/header["PixelX"], header["BeamY"]/header["PixelY"]
    # PILATSU 6M
    if shape == (2527, 2463):
        for j in xrange(12):
            for i in xrange(5):
                panel += 1
                s += """\
%(panel).2d/min_fs = %(min_fs)d
%(panel).2d/max_fs = %(max_fs)d
%(panel).2d/min_ss = %(min_ss)d
%(panel).2d/max_ss = %(max_ss)d
%(panel).2d/corner_x = %(cornerx).2f ; units are pixels
%(panel).2d/corner_y = %(cornery).2f
%(panel).2d/fs = -x
%(panel).2d/ss = -y

""" % dict(panel=panel, cornerx=corner[0], cornery=corner[1],
           min_fs=j*212, max_fs=j*212+195-1,
           min_ss=i*494, max_ss=i*494+487-1
           )

    open(geom_out, "w").write(s)
# make_geom_decomposed_for_raw()

def make_beam(beam_out):
    s = """\
beam/fluence = 7.0e10
beam/radius = 7.7e-7
beam/photon_energy = /LCLS/photon_energy_eV
beam/bandwidth = 5.2e-3
beam/divergence = 4.0e-4
profile_radius = 0.001e9
"""
    open(beam_out, "w").write(s)
# make_beam()

def decompose_panels(data):
    new_data = None
    # PILATSU 6M
    if data.shape == (2527, 2463):
        new_data = numpy.zeros(shape=(195*12, 487*5), dtype=data.dtype)
        for j in xrange(12):
            for i in xrange(5):
                new_data[195*j:195*(j+1), 487*i:487*(i+1)] = data[j*212:j*212+195, i*494:i*494+487]

    return new_data
# decompose_panels()

def convert(cbfin, params):
    if params.byteoffset:
        import yamtbx_byteoffset_h5_ext
        import pyublas

    h5out = os.path.basename(cbfin) + ".h5"

    data, ndimfast, ndimmid = cbf.load_minicbf_as_numpy(cbfin)
    data = data.reshape((ndimmid, ndimfast))
    if params.decompose:
        data = decompose_panels(data)
    header = XIO.Image(cbfin).header

    of = h5py.File(h5out, "w")

    grp = of.create_group("LCLS")
    dset = grp.create_dataset("photon_energy_eV", (1,), dtype=numpy.float)
    dset[...] = 12398.4/header["Wavelength"]
    dset = grp.create_dataset("photon_wavelength_A", (1,), dtype=numpy.float)
    dset[...] = header["Wavelength"]
    dset = grp.create_dataset("adu_per_eV", (1,), dtype=numpy.float)
    dset[...] = header["Wavelength"] / 12398.4
    dset = grp.create_dataset("detector_distance_m", (1,), dtype=numpy.float)
    dset[...] = header["Distance"]/1000.
    dset = grp.create_dataset("beam_xy_px", (2,), dtype=numpy.float)
    dset[...] = (header["BeamX"]/header["PixelX"], header["BeamY"]/header["PixelY"])
    dset = grp.create_dataset("osc_step_deg", (1,), dtype=numpy.float)
    dset[...] = header["PhiWidth"]

    cbfin = os.path.abspath(cbfin)
    dset = grp.create_dataset("original_file", (1,), "S%d"%len(cbfin))
    dset[...] = cbfin

    grp = of.create_group("data")
    if params.byteoffset:
        assert data.dtype == numpy.int32
        # hid_t group_id, const std::string &name, const pyublas::numpy_vector<int> &data, int width, int height)
        yamtbx_byteoffset_h5_ext.write_byteoffset_data(grp.id.id, "data", data.ravel(), data.shape[1], data.shape[0])
    else:
        dset = grp.create_dataset("data", data.shape, dtype=data.dtype)#, compression="CBF")
        dset[...] = data

    of.close()

    print "Processed: %s" % os.path.basename(cbfin)

    return data.shape
# convert()

def run(cbf_files, params):
    print "Attention - assuming cbf files given belong to a single dataset"
    print
    print "%d cbf files were given." % len(cbf_files)
    print

    if params.byteoffset:
        import yamtbx_byteoffset_h5_ext
        import pyublas

    last_shape = easy_mp.pool_map(fixed_func=lambda x: convert(x, params),
                                  args=cbf_files,
                                  processes=params.nproc)[-1]

    if params.decompose:
        make_geom_decomposed(XIO.Image(cbf_files[0]).header, last_shape, params.geom_out)
    else:
        make_geom(XIO.Image(cbf_files[0]).header, params.geom_out)

    make_beam(params.beam_out)

    print "Done."
    print
    print "Check %s and %s!" % (params.geom_out, params.beam_out)

# run()

if __name__ == "__main__":
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    if len(args) == 1 and args[0].endswith(".lst"):
        cbf_files = []
        for l in open(args[0]):
            if "#" in l: l = l[:l.index("#")]
            cbf_files.append(l.strip())
    else:
        cbf_files = filter(lambda x: x.endswith(".cbf"), args)

    if len(cbf_files) == 0:
        print "Usage: %s cbf_files" % sys.argv[0]
        quit()

    run(cbf_files, params)
