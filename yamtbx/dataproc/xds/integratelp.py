"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
"""
NOTE on unit cell constraints determination:
  XDS doesn't handle "real" rhombohedral space group (right?).
  So, No need to support R3 or R32. They are handled as H3 or H32, maybe.
"""

import re
import numpy
import collections

from cctbx import sgtbx

def rotations_to_missetting_angles(vals):
    a = numpy.array(vals)
    t = numpy.deg2rad(numpy.linalg.norm(a))
    u = a / numpy.linalg.norm(a)
    ct, st = numpy.cos(t), numpy.sin(t)
    ux, uy, uz = u
    R = numpy.matrix([[ct+ux*ux*(1-ct),    ux*uy*(1-ct)-uz*st, ux*uz*(1-ct)+uy*st],
                      [uy*uz*(1-ct)+uz*st, ct+uy*uy*(1-ct),    uy*uz*(1-ct)-ux*st],
                      [uz*ux*(1-ct)-uy*st, uz*uy*(1-ct)+ux*st, ct+uz*uz*(1-ct)   ]
                      ])

    phi = numpy.zeros(3) # missetting angles

    if 1. - numpy.abs(R[2,0]) < 0.0000001:
        phi[0] = 0.
        phi[1] = numpy.sin(-R[2,0]) * numpy.pi/2.
        phi[2] = numpy.arctan2(-R[0,1], R[1,1])
    else:
        x = -R[2,0]
        if abs(x) > 1: x = numpy.sign(x)
        phi[1] = numpy.arcsin(x)
        phi[2] = numpy.arctan2(R[1,0], R[0,0])
        phi[0] = numpy.arctan2(R[2,1], R[2,2])

    #print "org:", vals
    #print "phi:", numpy.rad2deg(phi)


    return numpy.rad2deg(phi)
# rotations_to_missetting_angles()

class IntegrateLp:
    def __init__(self, lpin):
        if lpin is not None:
            self.parse(lpin)
    # __init__()

    def parse(self, int_lp):
        re_im = re.compile("^ (.....)   0 +([0-9\.]+) +([0-9]+) +([0-9]+) +([0-9]+) +([0-9]+) +([0-9]+) +([0-9\.]+) +([0-9\.]+)")
        re_cell = re.compile("^ UNIT CELL PARAMETERS *([0-9\.]+) *([0-9\.]+) *([0-9\.]+) *([0-9\.]+) *([0-9\.]+) *([0-9\.]+)")
        re_rotation = re.compile("^ CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION *([-0-9\.]+) *([-0-9\.]+) *([-0-9\.]+)") #
        re_mosaicity = re.compile("^ CRYSTAL MOSAICITY \(DEGREES\) *([0-9\.]+)") #
        re_axis = re.compile("^ LAB COORDINATES OF ROTATION AXIS *([-0-9\.]+) *([-0-9\.]+) *([-0-9\.]+)") #
        re_beam = re.compile("^ DIRECT BEAM COORDINATES \(REC\. ANGSTROEM\) *([-0-9\.]+) *([-0-9\.]+) *([-0-9\.]+)") #
        re_dist = re.compile("^ CRYSTAL TO DETECTOR DISTANCE \(mm\) *([-0-9\.]+)")
        re_dev_spot = re.compile("^ STANDARD DEVIATION OF SPOT    POSITION \(PIXELS\) *([0-9\.]+)")
        re_dev_spindle = re.compile("^ STANDARD DEVIATION OF SPINDLE POSITION \(DEGREES\) *([0-9\.]+)")
        re_orig = re.compile("^ DETECTOR ORIGIN \(PIXELS\) AT *([0-9\.]+) *([0-9\.]+)")

        images = [] # as key of params
        self.cell_changes = []
        self.blockparams = collections.OrderedDict()
        clear_flag = False

        self.frames = []
        self.scales, self.overloads, self.strongs, self.rejecteds, self.sigmads, self.sigmars = [], [], [], [], [], []

        self.space_group = None

        # Read INTEGRATE.LP file
        for l in open(int_lp):
            r_im = re_im.search(l)
            r_cell = re_cell.search(l)
            r_rotation = re_rotation.search(l)
            r_dist = re_dist.search(l)
            r_spot = re_dev_spot.search(l)
            r_spindle = re_dev_spindle.search(l)
            r_orig = re_orig.search(l)

            if l.startswith(" SPACE_GROUP_NUMBER="):
                sgnum = int(l.strip().split()[-1])
                if sgnum > 0:
                    self.space_group = sgtbx.space_group_info(sgnum).group()

            if r_im:
                if clear_flag:
                    images = []
                    clear_flag = False
                image, scale, nbkg, novl, newald, nstrong, nrej, sigmad, sigmar = r_im.groups()
                images.append(int(image))

                # for plot
                self.frames.append(int(image))
                self.scales.append(scale)
                self.overloads.append(int(novl))
                self.strongs.append(int(nstrong))
                self.rejecteds.append(int(nrej))
                self.sigmads.append(sigmad)
                self.sigmars.append(sigmar)

            if r_cell:
                #a, b, c, alpha, beta, gamma = r_cell.groups()
                self.blockparams.setdefault(tuple(images), {})["cell"] = r_cell.groups()
                self.cell_changes.append((images, r_cell.groups()))
                clear_flag = True

            if r_rotation:
                self.blockparams.setdefault(tuple(images), {})["rotation"] = r_rotation.groups()
                misset = rotations_to_missetting_angles(map(float, r_rotation.groups()))
                self.blockparams.setdefault(tuple(images), {})["misset"] = map(lambda x:"%.2f"%x, misset)
                clear_flag = True

            if r_dist:
                self.blockparams.setdefault(tuple(images), {})["dist"] = r_dist.group(1)
                clear_flag = True
            if r_spot:
                self.blockparams.setdefault(tuple(images), {})["spot"] = r_spot.group(1)
                clear_flag = True
            if r_spindle:
                self.blockparams.setdefault(tuple(images), {})["spindle"] = r_spindle.group(1)
                clear_flag = True
            if r_orig:
                self.blockparams.setdefault(tuple(images), {})["orig"] = r_orig.groups()
                clear_flag = True

            if l.startswith(" SIGMAB (degree)"):
                self.blockparams.setdefault(tuple(images), {})["sigmab9"] = l.strip().split()[-9:]
                clear_flag = True

            if l.startswith(" SIGMAR (degree)"):
                self.blockparams.setdefault(tuple(images), {})["sigmar9"] = l.strip().split()[-9:]
                clear_flag = True

    # parse_integrate_lp()

# class IntegrateLp
