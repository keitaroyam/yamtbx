#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

"""
Convert XPARM.XDS information for labelit programs.

Use PHENIX_TRUST_OTHER_ENV=1 to start this script.
"""

import re, sys, os, pickle, math
import numpy
from cctbx import crystal_orientation
import cctbx.sgtbx.bravais_types
from iotbx.detectors.context.endstation import EndStation
import iotbx.detectors
from labelit.symmetry.subgroup import MetricSubgroup
from iotbx.detectors.context.config_detector import beam_center_convention_from_image_object
from scitbx import matrix

from yamtbx.dataproc.xds import *
from yamtbx.dataproc.xds.xparm import XPARM

ftype = numpy.float64

"""
import pickle
H = pickle.load(open("LABELIT_possible"))
G = pickle.load(open("LABELIT_pickle"))
"""

def get_spot_convention(imagefile):
    # Maybe 2 if RAXIS else 0..
    return 2
    from spotfinder.command_line.signal_strength import master_params
    phil_params = master_params.extract()

    imageobject = iotbx.detectors.ImageFactory(G["file"][start_frame])

    beam_center_convention_from_image_object(imageobject, phil_params)

    print "SPOT CONVENTION=", phil_params.spot_convention

    return phil_params.spot_convention
# get_spot_convention()

class XPARM_to_labelit(XPARM):
    def __init__(self, xparm_file):
        XPARM.__init__(self, xparm_file)
        self.bin = 1
        self.set_bin()
    # __init__()

    def set_bin(self):
        # By default, labelit calls setBin(2) if size1 > 4000.
        # See: $PHENIX/cctbx_project/spotfinder/diffraction/imagefiles.py +245

        if self.nx > 4000 or self.ny > 4000:
            self.bin = 2

        # We must *NOT* modify qx,qy, nx,ny since it would be used later. (e.g. beam center conversion)
        print "BINNING=", self.bin, "QX, QY, NX, NY=", self.qx*self.bin, self.qy*self.bin, self.nx//self.bin, self.ny//self.bin

    # set_bin()

    def get_bravais(self):
        return str(cctbx.sgtbx.bravais_types.bravais_lattice(1))
        return str(cctbx.sgtbx.bravais_types.bravais_lattice(self.spacegroup))
    # get_bravais()

    def get_system(self):
        return cctbx.sgtbx.space_group_info(1).group().crystal_system().lower()
        return cctbx.sgtbx.space_group_info(self.spacegroup).group().crystal_system().lower()
    # get_system()

    def get_labelit_orient(self):
        UNI_I = matrix.sqr((0,1,0,1,0,0,0,0,-1)).inverse()
        a, b, c = tuple(self.a_axis), tuple(self.b_axis), tuple(self.c_axis)
        matXDS = matrix.sqr(a+b+c)
        print "matXDS=", matXDS[4]
        matRossmann = (UNI_I * matXDS.transpose()).transpose()
        orient = crystal_orientation.crystal_orientation(matRossmann, False) # reciprocal flag

        rotation_ax = self.get_endstation().rot_axi
        orient = orient.rotate_thru(rotation_ax, -self.starting_angle*math.pi/180.)
        return orient
    # get_labelit_orient()

    def get_labelit_xy_beam(self):
        # Definition of beam center is different between XDS and MOSFLM!
        n = numpy.array((0.,0.,1.)) # FIXME. Not always true. needs cross_prod(DIRECTION_OF_DETECTOR_X-AXIS=, DIRECTION_OF_DETECTOR_Y-AXIS=)
        b = self.incident_beam
        offset = abs(self.distance) * (1./numpy.dot(b,n) * b - n)
        print "BEAM CENTER OFFSET=", offset[0]/self.qx, offset[1]/self.qy
        return self.origin[1]*self.qy + offset[1], self.origin[0]*self.qx + offset[0]
    # get_labelit_xy_beam()

    def get_endstation(self):
        UNI_I = matrix.sqr((0,1,0,1,0,0,0,0,-1)).inverse()
        e = EndStation()
        #e.set_camera_convention() always 1?
        e.set_rotation_axis(UNI_I*matrix.col(self.rotation_axis))
        print "endstation.rot_axi=", e.rot_axi
        print "endstation.cam_con", tuple(e.cam_con)
        return e
    # get_endstation()

    def get_distance(self):
        return abs(self.distance)

    def get_pixel_size(self):
        return self.qx * self.bin

    def get_size1_size2(self):
        return self.nx // self.bin, self.ny // self.bin  # FIXME. reversed?? OK?

# class XPARM_to_labelit

def template_to_filename(img_template, iframe):
    re_var = re.compile("\?+")

    # like: "????"
    var = re_var.search(img_template).group()
    repl = "%%.%dd" % len(var)

    # replace e.g. ???? => 0001
    return img_template.replace(var, repl % iframe)
# template_to_filenames()


if __name__ == "__main__":

    #wdir = "/Users/yam/crystal/Lysozyme/Lyz080513B/xds_process_Lyz"
    wdir = sys.argv[1]

    inp = dict(get_xdsinp_keyword(os.path.join(wdir, "XDS.INP")))
    xparm = XPARM_to_labelit(os.path.join(wdir, "GXPARM.XDS"))

    G = {}
    G["pixel_size"] = xparm.get_pixel_size()
    start_frame = int(inp["DATA_RANGE"].split()[0])
    G["file"] = {start_frame: template_to_filename(inp["NAME_TEMPLATE_OF_DATA_FRAMES"], start_frame)}
    G["size1"], G["size2"] = xparm.get_size1_size2()

    #img = iotbx.detectors.ImageFactory(G["file"][start_frame])

    G["spot_convention"] = get_spot_convention(G["file"][start_frame])
    G["endstation"] = xparm.get_endstation()
    G["recommended_grid_sampling"] = 0.018128529847134645 #0.029 # Tekitou
    G["xbeam"], G["ybeam"] = xparm.get_labelit_xy_beam()
    G["distance"] = xparm.get_distance() #inp["DETECTOR_DISTANCE"]
    G["twotheta"] = 0.
    G['wavelength'] = inp["X-RAY_WAVELENGTH"]
    G['ref_maxcel'] = 362.259397 #114.373193 # FIXME? Maybe doesn't matter.
    G['deltaphi'] = inp["OSCILLATION_RANGE"]
    G["best_integration"] = {}
    G["best_integration"]["mosaicity"] = 0.15 # FIXME? Maybe doesn't matter.
    G["best_integration"]["orient"] = xparm.get_labelit_orient()
    print "orient:", G["best_integration"]["orient"]
    print G


    H = []
    h = MetricSubgroup()
    h["counter"] = 1
    h["orient"]  = xparm.get_labelit_orient()
    #h["mosaicity"] = 0.075 # FIXME
    h["bravais"] = xparm.get_bravais() # like oP
    # Needed?
    h["status"] = "ok"
    h["residual"] = 0.1 # FIXME
    h["system"] = xparm.get_system()
    h["refined x beam"], h["refined y beam"] = xparm.get_labelit_xy_beam()
    h["refined distance"] = xparm.get_distance()
    h["max_angular_difference"] = 0.3 # FIXME
    h["count_GOOD"] = 239 # FIXME
    H.append(h)
    print H

    pickle.dump(G, open(os.path.join(wdir, "LABELIT_pickle"), "w"))
    pickle.dump(H, open(os.path.join(wdir, "LABELIT_possible"), "w"))


    print
    print "Run:"
    print "labelit.precession_photo bravais_choice=1 image_range=%d,%d intensity_full_scale=512 plot_section=H,K,0 pdf_output.file=HK0_xds.pdf" % tuple(map(lambda x:int(x), inp["DATA_RANGE"].split()))
    print "labelit.precession_photo bravais_choice=1 image_range=%d,%d intensity_full_scale=512 plot_section=0,K,L pdf_output.file=0KL_xds.pdf" % tuple(map(lambda x:int(x), inp["DATA_RANGE"].split()))
    print "labelit.precession_photo bravais_choice=1 image_range=%d,%d intensity_full_scale=512 plot_section=H,0,L pdf_output.file=H0L_xds.pdf" % tuple(map(lambda x:int(x), inp["DATA_RANGE"].split()))
