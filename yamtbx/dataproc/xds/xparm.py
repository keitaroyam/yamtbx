"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
import os
import numpy
import itertools
import copy
from yamtbx.dataproc.xds import get_xdsinp_keyword
from yamtbx.util import safe_float

class Segment(object):
    def __init__(self):
        self.x1, self.x2, self.y1, self.y2 = 1, 1, 1, 1
        self.orgxs, self.orgys, self.fs = 0, 0, 0
        self.eds_x = numpy.array([1.,0,0])
        self.eds_y = numpy.array([0,1.,0])

class XPARM(object):
    def __init__(self, xparm_file=None):
        self.segments = []
        self.xds_ver = "March 30, 2013"
        if xparm_file is not None:
            self.parse_xparm_file(xparm_file)
        else:
            self.starting_frame = 1
            self.starting_angle = 0.
            self.osc_range = -1.
            self.rotation_axis = numpy.array((1.,0.,0.))
            self.wavelength = -1.
            self.incident_beam = numpy.array((0.,0.,1.))
            self.nx, self.ny = 0, 0
            self.qx, self.qy = 0., 0.
            self.distance = 0.
            self.origin = numpy.array((0., 0.))
            self.X_axis = numpy.array((1., 0., 0.))
            self.Y_axis = numpy.array((0., 1., 0.))
            self.Z_axis = numpy.array((0., 0., 1.))
            self.spacegroup = 1
            self.unit_cell = numpy.array((100., 100., 100., 90., 90., 90.))
            self.a_axis = numpy.array((100., 0., 0.))
            self.b_axis = numpy.array((0., 100., 0.))
            self.c_axis = numpy.array((0., 0., 100.))
    # __init__()

    def space_group_info(self):
        from cctbx import sgtbx
        g = sgtbx.space_group_info(self.spacegroup).group()
        return g.info()

    def n_segments(self):
        return len(self.segments) if self.segments else 1

    def parse_xparm_file(self, xparm_file):
        lines = open(xparm_file).readlines()

        is_new_format = "XPARM.XDS" in lines[0]

        if not is_new_format:
            starting_frame, starting_angle, osc_range, rotx, roty, rotz = lines[0].split()
            wavelength, ibeamx, ibeamy, ibeamz = lines[1].split()
            nx, ny, qx, qy = lines[2].split()
            distance, orgx, orgy = lines[3].split()
            Xx, Xy, Xz = lines[4].split()
            Yx, Yy, Yz = lines[5].split()
            Zx, Zy, Zz = lines[6].split()
            spacegroup, a, b, c, alpha, beta, gamma = lines[7].split()
            ax, ay, az = lines[8].split()
            bx, by, bz = lines[9].split()
            cx, cy, cz = lines[10].split()
        else:
            self.xds_ver = lines[0][22:].strip()
            starting_frame, starting_angle, osc_range, rotx, roty, rotz = lines[1].split()
            wavelength, ibeamx, ibeamy, ibeamz = lines[2].split()
            spacegroup, a, b, c, alpha, beta, gamma = lines[3].split()
            ax, ay, az = lines[4][:15], lines[4][15:30], lines[4][30:45]
            bx, by, bz = lines[5][:15], lines[5][15:30], lines[5][30:45]
            cx, cy, cz = lines[6][:15], lines[6][15:30], lines[6][30:45]
            nseg, nx, ny, qx, qy = lines[7].split()
            orgx, orgy, distance = lines[8].split()
            Xx, Xy, Xz = lines[9].split()
            Yx, Yy, Yz = lines[10].split()
            Zx, Zy, Zz = lines[11].split()

        self.starting_frame = int(starting_frame)
        self.starting_angle = float(starting_angle)
        self.osc_range = float(osc_range)
        self.rotation_axis = numpy.array((float(rotx), float(roty), float(rotz)))
        self.wavelength = float(wavelength)
        self.incident_beam = numpy.array((float(ibeamx), float(ibeamy), float(ibeamz)))
        self.nx = float(nx)
        self.ny = float(ny)
        self.qx = float(qx)
        self.qy = float(qy)
        self.distance = float(distance)
        self.origin = numpy.array((float(orgx), float(orgy)))
        self.X_axis = numpy.array((float(Xx), float(Xy), float(Xz)))
        self.Y_axis = numpy.array((float(Yx), float(Yy), float(Yz)))
        self.Z_axis = numpy.array((float(Zx), float(Zy), float(Zz)))
        self.spacegroup = int(spacegroup)
        self.unit_cell = numpy.array((float(a), float(b), float(c), float(alpha), float(beta), float(gamma)))
        self.a_axis = numpy.array(list(map(safe_float, (ax, ay, az))))
        self.b_axis = numpy.array(list(map(safe_float, (bx, by, bz))))
        self.c_axis = numpy.array(list(map(safe_float, (cx, cy, cz))))

        # Segments
        if is_new_format and len(lines) > 14:
            for li in range(12, len(lines), 2):
                self.segments.append(Segment())

                s = lines[li]
                iseg, x1, x2, y1, y2 = list(map(int, (s[:10], s[10:20], s[20:30], s[30:40], s[40:50])))
                self.segments[-1].x1 = x1
                self.segments[-1].x2 = x2
                self.segments[-1].y1 = y1
                self.segments[-1].y2 = y2
                
                s = lines[li+1]
                orgxs, orgys, fs, x0, x1, x2, y0, y1, y2 = list(map(float, (s[:8], s[8:16], s[16:24],
                                                                       s[24:33], s[33:42], s[42:51],
                                                                       s[51:60], s[60:69], s[69:78])))
                self.segments[-1].orgxs = orgxs
                self.segments[-1].orgys = orgys
                self.segments[-1].fs = fs
                self.segments[-1].eds_x = numpy.array([x0,x1,x2])
                self.segments[-1].eds_y = numpy.array([y0,y1,y2])
    # parse_xparm_file()

    def set_info_from_xdsinp_or_inpstr(self, xdsinp=None, inpstr=None):
        assert (xdsinp,inpstr).count(None) == 1
        t1 = lambda x: x.split()[0] # may have units that should be removed (if read from INTEGRATE.LP header)
        
        table = [("STARTING_FRAME", "starting_frame", lambda x: int(t1(x))),
                 ("STARTING_ANGLE", "starting_angle", lambda x: float(t1(x))),
                 ("OSCILLATION_RANGE", "osc_range", lambda x: float(t1(x))),
                 ("ROTATION_AXIS", "rotation_axis", lambda x: numpy.array([float(y) for y in x.split()[:3]])),
                 ("DIRECTION_OF_DETECTOR_X-AXIS", "X_axis", lambda x: numpy.array([float(y) for y in x.split()[:3]])),
                 ("DIRECTION_OF_DETECTOR_Y-AXIS", "Y_axis", lambda x: numpy.array([float(y) for y in x.split()[:3]])),
                 ("X-RAY_WAVELENGTH", "wavelength", lambda x: float(t1(x))),
                 ("INCIDENT_BEAM_DIRECTION", "incident_beam", lambda x: numpy.array([float(y) for y in x.split()[:3]])),
                 ("NX", "nx", lambda x: int(t1(x))),
                 ("NY", "ny", lambda x: int(t1(x))),
                 ("QX", "qx", lambda x: float(t1(x))),
                 ("QY", "qy", lambda x: float(t1(x))),
                 ("DETECTOR_DISTANCE", "distance", lambda x: float(t1(x))),
                 ("SPACE_GROUP_NUMBER", "spacegroup", lambda x: int(t1(x))),
                 ("UNIT_CELL_CONSTANTS", "unit_cell", lambda x: numpy.array([float(y) for y in x.split()[:6]])),
                 ("UNIT_CELL_A-AXIS", "a_axis", lambda x: numpy.array([float(y) for y in x.split()[:3]])),
                 ("UNIT_CELL_B-AXIS", "b_axis", lambda x: numpy.array([float(y) for y in x.split()[:3]])),
                 ("UNIT_CELL_C-AXIS", "c_axis", lambda x: numpy.array([float(y) for y in x.split()[:3]]))
                 ]
        inp_raw = get_xdsinp_keyword(xdsinp=xdsinp, inp_str=inpstr)
        inp = dict(inp_raw)# I believe dict() removes duplicated parameters and keeps last.

        for k, at, f in table:
            if k in inp and inp[k].strip() != "":
                setattr(self, at, f(inp[k]))
        if "ORGX" in inp:
            self.origin[0] = float(inp["ORGX"])
        if "ORGY" in inp:
            self.origin[1] = float(inp["ORGY"])

        if "DIRECTION_OF_DETECTOR_X-AXIS" in inp:
            self.Z_axis = numpy.cross(self.X_axis, self.Y_axis)
            self.Z_axis /= numpy.linalg.norm(self.Z_axis)

        # Segment
        for k, v in inp_raw:
            if k == "SEGMENT":
                sp = list(map(int, v.split()))
                self.segments.append(Segment())
                self.segments[-1].x1 = sp[0]
                self.segments[-1].x2 = sp[1]
                self.segments[-1].y1 = sp[2]
                self.segments[-1].y2 = sp[3]
            elif k == "SEGMENT_ORGX":
                self.segments[-1].orgxs = float(v)
            elif k == "SEGMENT_ORGY":
                self.segments[-1].orgys = float(v)
            elif k == "SEGMENT_DISTANCE":
                self.segments[-1].fs = float(v)
            elif k == "DIRECTION_OF_SEGMENT_X-AXIS":
                self.segments[-1].eds_x = numpy.array(list(map(float, v.split())))               
            elif k == "DIRECTION_OF_SEGMENT_Y-AXIS":
                self.segments[-1].eds_y = numpy.array(list(map(float, v.split())))               
    # set_info_from_xdsinp()


    def parse_integratelp_header(self, lpin):
        #get_xdsinp_keyword()
        inp_str = ""
        for l in open(lpin):
            inp_str += l
            if "PROCESSING OF IMAGES" in l: break

        self.set_info_from_xdsinp_or_inpstr(inpstr=inp_str)
    # parse_integratelp_header()

    def set_info_from_xdsinp(self, xdsinp):
        self.set_info_from_xdsinp_or_inpstr(xdsinp=xdsinp)
    # set_info_from_xdsinp()

    def xparm_str(self, old_format=False):
        assert not old_format # Currently, only new format is supported!
        xparm_str = """ XPARM.XDS    VERSION %s
%6d%14.4f%10.4f%10.6f%10.6f%10.6f
%15.6f%15.6f%15.6f%15.6f
%6d%12.4f%12.4f%12.4f%8.3f%8.3f%8.3f
%15.6f%15.6f%15.6f
%15.6f%15.6f%15.6f
%15.6f%15.6f%15.6f
%10d%10d%10d%12.6f%12.6f
%15.6f%15.6f%15.6f
%15.6f%15.6f%15.6f
%15.6f%15.6f%15.6f
%15.6f%15.6f%15.6f
""" % (self.xds_ver,
       self.starting_frame, self.starting_angle, self.osc_range, self.rotation_axis[0], self.rotation_axis[1], self.rotation_axis[2], 
       self.wavelength, self.incident_beam[0], self.incident_beam[1], self.incident_beam[2], 
       self.spacegroup, self.unit_cell[0], self.unit_cell[1], self.unit_cell[2], self.unit_cell[3], self.unit_cell[4], self.unit_cell[5],
       self.a_axis[0], self.a_axis[1], self.a_axis[2], 
       self.b_axis[0], self.b_axis[1], self.b_axis[2], 
       self.c_axis[0], self.c_axis[1], self.c_axis[2],
       self.n_segments(), self.nx, self.ny, self.qx, self.qy,
       self.origin[0], self.origin[1], self.distance, 
       self.X_axis[0], self.X_axis[1], self.X_axis[2], 
       self.Y_axis[0], self.Y_axis[1], self.Y_axis[2], 
       self.Z_axis[0], self.Z_axis[1], self.Z_axis[2], 
       )

        if not self.segments:
            xparm_str += """\
%10d%10d%10d%10d%10d
%8.2f%8.2f%8.2f%9.5f%9.5f%9.5f%9.5f%9.5f%9.5f
""" % (1, 1, self.nx, 1, self.ny,
       0., 0., 0., 1., 0., 0., 0., 1., 0.)
        else:
            for iseg, seg in enumerate(self.segments):
                xparm_str += """\
%10d%10d%10d%10d%10d
%8.2f%8.2f%8.2f%9.5f%9.5f%9.5f%9.5f%9.5f%9.5f
""" % (iseg+1, seg.x1, seg.x2, seg.y1, seg.y2,
       seg.orgxs, seg.orgys, seg.fs,
       seg.eds_x[0], seg.eds_x[1], seg.eds_x[2],
       seg.eds_y[0], seg.eds_y[1], seg.eds_y[2])
        
        return xparm_str
    # xparm_str()

    def crystal_symmetry(self):
        from cctbx import crystal
        return crystal.symmetry(tuple(self.unit_cell),
                                self.spacegroup)
    # crystal_symmetry()

    def update_cell_based_on_axes(self):
        from yamtbx.util.maths import vectors_angle

        a, b, c = [numpy.linalg.norm(x) for x in (self.a_axis, self.b_axis, self.c_axis)]
        al = vectors_angle(self.b_axis, self.c_axis)
        be = vectors_angle(self.c_axis, self.a_axis)
        ga = vectors_angle(self.a_axis, self.b_axis)
        al, be, ga = [numpy.rad2deg(x) for x in (al, be, ga)]
        
        self.unit_cell = numpy.array((a,b,c,al,be,ga))
        print() 
# class XPARM

def prep_xparm_objects_from_integrate_lp(lpfile, xparm_ref=None):
    if xparm_ref:
        xp_ref = XPARM(xparm_ref)
    else:
        xp_ref = XPARM()
        xp_ref.parse_integratelp_header(lpfile)
    
    keys = {"beam direction": "DIRECT BEAM COORDINATES (REC. ANGSTROEM)",
            "beam center": "DETECTOR ORIGIN (PIXELS) AT",
            "distance": "CRYSTAL TO DETECTOR DISTANCE (mm)",
            "rotation axis": "LAB COORDINATES OF ROTATION AXIS",
            "a axis": "COORDINATES OF UNIT CELL A-AXIS",
            "b axis": "COORDINATES OF UNIT CELL B-AXIS",
            "c axis": "COORDINATES OF UNIT CELL C-AXIS",
            "cell": "UNIT CELL PARAMETERS",
            "spacegroup": "SPACE GROUP NUMBER"
            }
    all_data = []

    flag_read = False
    for l in open(lpfile):
        if "PROCESSING OF IMAGES" in l:
            flag_read = False
            l = l.strip()
            range_current = [int(x.strip()) for x in l[l.index("PROCESSING OF IMAGES")+len("PROCESSING OF IMAGES"):].split("...")]
            
            flag_read = True
            all_data.append([tuple(range_current), {}])
        elif "STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED" in l:
            flag_read = False
        elif flag_read:
            for key, s in list(keys.items()):
                if s in l:
                    l = l.strip()
                    val = [float(x.strip()) for x in l[l.index(s)+len(s):].split()]
                    all_data[-1][1][key] = val

    ret = []

    for r, data in all_data:
        xp = copy.copy(xp_ref)
        
        if "rotation axis" in data: xp.rotation_axis = numpy.array(data["rotation axis"])
        if "beam direction" in data: xp.incident_beam = numpy.array(data["beam direction"])
        if "spacegroup" in data: xp.spacegroup = int(data["spacegroup"][0])
        if "cell" in data: xp.unit_cell = numpy.array(data["cell"])
        if "a axis" in data: xp.a_axis = numpy.array(data["a axis"])
        if "b axis" in data: xp.b_axis = numpy.array(data["b axis"])
        if "c axis" in data: xp.c_axis = numpy.array(data["c axis"])
        if "beam center" in data: xp.origin = numpy.array(data["beam center"])
        if "distance" in data: xp.distance = data["distance"][0]

        ret.append((r, xp))

    return ret
# prep_xparm_objects_from_integrate_lp()

def get_xparm_from_integrate_lp(lpfile, frame):
    assert 0 < frame

    xparm_objs = prep_xparm_objects_from_integrate_lp(lpfile)
    xp = [x for x in xparm_objs if x[0][0] <= frame <= x[0][1]][0][1]

    return xp.xparm_str()
# get_xparm_from_integrate_lp()


def find_best_ucaxes_transform(xparm1, xparm2):
    """
    Reference: matrix iteration method in whirligig in CrystFEL-0.6.0 src/whirligig.c: gatinator()

    calculates the matrix m, which gives the best result; axes1 and m * axes2
    returns m, axes1, axes2, rms
    """

    a1, b1, c1 = xparm1.a_axis, xparm1.b_axis, xparm1.c_axis
    a2, b2, c2 = xparm2.a_axis, xparm2.b_axis, xparm2.c_axis

    P = numpy.vstack([a2, b2, c2])
    Q = numpy.vstack([a1, b1, c1])

    min_score = numpy.inf
    best_mat = None

    for x in itertools.product(range(-1,2), repeat=9):
        #if x != (1,0,0,0,1,0,0,0,1): continue
        m = numpy.array(x).reshape(3,3)
        if numpy.linalg.det(m) != 1: continue
        #print "Trying"
        #print m

        score = numpy.sum((numpy.dot(m, Q) - P)**2)
        #if x == (1,0,0,0,1,0,0,0,1): print score
        if score < min_score:
            min_score = score
            best_mat = m

    return best_mat, Q, P, numpy.sqrt(min_score/3.)
# find_best_ucaxes_transform()
