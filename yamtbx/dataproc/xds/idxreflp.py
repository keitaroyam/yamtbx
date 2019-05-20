"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import re
import math
import numpy
import copy
from yamtbx.dataproc.xds.xparm import XPARM
from yamtbx.dataproc.xds import get_xdsinp_keyword
from cctbx import uctbx

re_outof = re.compile("^ *([0-9]+) OUT OF *([0-9]+) SPOTS INDEXED.")

class IdxrefLp:

    def __init__(self, lpin):
        self.subtree_population = []
        self.clusters = []
        self.reduced_cell = None

        if lpin is not None:
            self.parse(lpin)
    # __init__()

    def parse(self, lpin):
        reading = ""
        self.subtree_population = []
        self.clusters = []
        self.reduced_cell = None

        lines = open(lpin).readlines()
        for i, l in enumerate(lines):
            if "OUT OF" in l:
                r = re_outof.search(l)
                if r: outof = map(int, r.groups())

            elif "#  COORDINATES OF VECTOR CLUSTER   FREQUENCY       CLUSTER INDICES" in l:
                reading = "clusters"
            elif reading == "clusters":
                if l.strip()=="":
                    reading = ""
                    continue
                x, y, z, freq, h, k, l = map(float, (l[6:16], l[16:26], l[26:36],
                                                     l[37:46], l[47:56] , l[56:66], l[66:76]))
                self.clusters.append(((x,y,z), freq, (h,k,l)))
            elif "PARAMETERS OF THE REDUCED CELL (ANGSTROEM & DEGREES)" in l:
                reading = "reduced_cell"
            elif reading == "reduced_cell":
                cell = map(float, (l[:10], l[10:20], l[20:30], l[30:40], l[40:50], l[50:60]))
                self.reduced_cell = cell # uctbx.unit_cell complains if zero is included!
                reading = ""
            elif "SUBTREE    POPULATION" in l:
                reading = "subtree_population"
            elif "NUMBER OF ACCEPTED SPOTS FROM LARGEST SUBTREE" in l:
                reading = ""
            elif reading == "subtree_population" and l.strip() != "":
                sp = map(int, l.split())
                self.subtree_population.append(sp[1])
    # parse()

    def first_subtree_fraction(self):
        if len(self.subtree_population) == 0:
            return 0

        return self.subtree_population[0] / float(sum(self.subtree_population))
    # first_subtree_fraction()

    def cluster_integerness(self):
        # hist = [[100,0,0,0,0, 0,0,0,0,0], (for k), (for l)]
        hist = map(lambda x: numpy.zeros(6, dtype=numpy.int), xrange(3))
        
        for xyz, freq, hkl in self.clusters:
            fracs = map(lambda x: abs(math.modf(x)[0]), hkl) # fractional part
            fracs = map(lambda x: 1-x if x > 0.5 else x, fracs)
            fracs = map(lambda x: int(x*10+0.5), fracs) # take the first decimal place (after rounding)
            for i, h in enumerate(fracs):
                hist[i][h] += int(freq)
        return hist
    # cluster_integerness()

    def is_cell_maybe_half(self):
        if self.subtree_population == [] or self.clusters == []:
            return False

        intness = self.cluster_integerness()

        if any(map(lambda h: h[0]/sum(h) > 0.4 and h[5]/sum(h) > 0.4, intness)):
            sum_popu = sum(self.subtree_population)
            if all(map(lambda x:x/sum_popu > 0.4, self.subtree_population[:2])):
                return True

        return False
    # is_cell_maybe_too_short()

    def deduce_correct_cell_based_on_integerness(self):
        # currently only support doubled case..

        reduced_cell = copy.copy(self.reduced_cell)

        for i, h in enumerate(self.cluster_integerness()):
            sum_intness = sum(h)
            if h[0]/sum_intness > 0.4 and h[5]/sum_intness > 0.4:
                reduced_cell[i] *= 2
        
        return uctbx.unit_cell(reduced_cell)
    # deduce_correct_cell_based_on_integerness()

# class IdxrefLp

class SpotXds:
    def __init__(self, sptxdsin):
        self.items = []
        self.calc_d = None
        if sptxdsin is not None:
            self.parse(sptxdsin)
    # __init__()

    def parse(self, sptxdsin):
        self.items = []

        for l in open(sptxdsin):
            sp = l.split()
            x, y, z, intensity, iseg, h, k, l = (None,)*8

            if len(sp) == 4:
                x, y, z, intensity = sp
            elif len(sp) == 5:
                x, y, z, intensity, iseg = sp
            elif len(sp) == 7:
                x, y, z, intensity, h, k, l = sp
            elif len(sp) == 8:
                x, y, z, intensity, iseg, h, k, l = sp
                
            x, y, z, intensity = map(float, (x, y, z, intensity))
            if iseg is not None: iseg = int(iseg)
            if h is not None: h, k, l = map(int, (h, k, l))
                
            self.items.append(((x, y, z), intensity, iseg, (h, k, l)))
    # parse()
    
    def write(self, out, frame_selection=[]):
        for (i, (xyz, intensity, iseg, hkl)) in enumerate(self.items):
            if frame_selection and int(xyz[2])+1 not in frame_selection: continue
            out.write("% .2f % .2f % .2f" % xyz)
            out.write(" % .1f" % intensity)
            if iseg is not None: out.write(" %d" % iseg)
            if None not in hkl: out.write(" %d %d %d" % hkl)
            out.write("\n")
    # write()

    def collected_spots(self, with_resolution=True):
        if len(self.items) == 0: return []

        if with_resolution:
            assert self.calc_d is not None
            
        ret = []
        for xyz, intensity, iseg, hkl in self.items:
            if with_resolution:
                ret.append(xyz+(intensity, self.calc_d(*xyz[:2]),))
            else:
                ret.append(xyz+(intensity, -1,))
            
        return ret
    # collected_spots()

    def indexed_and_unindexed_by_frame(self):
        if len(self.items) == 0: return []
        data = {}
        for xyz, intensity, iseg, hkl in self.items:
            frame = int(xyz[2])+1
            if frame not in data:
                data[frame] = [0, 0]
            if hkl[0] is None or all(map(lambda h: h==0, hkl)):
                data[frame][1] += 1 # unindexed
            else:
                data[frame][0] += 1 # indexed

        ret = data.items()
        ret.sort(key=lambda x: x[0])
        return ret
    # indexed_and_unindexed_by_frame()

    def spots_by_frame(self):
        if len(self.items) == 0: return
        data = {}
        for xyz, intensity, iseg, hkl in self.items:
            frame = int(xyz[2])+1
            data[frame] = data.get(frame, 0) + 1

        return data
    # spots_by_frame()

    def indexed_and_unindexed_on_detector(self, with_resolution=True):
        if len(self.items) == 0: return

        data = {"indexed": [], "unindexed": []}

        if with_resolution:
            assert self.calc_d is not None
            

        for xyz, intensity, iseg, hkl in self.items:
            tmp = xyz[:2]
            if with_resolution:
                tmp += self.calc_d(*xyz[:2]),
            else:
                tmp += -1,

            if hkl[0] is None or all(map(lambda h: h==0, hkl)):
                data["unindexed"].append(tmp)
            else:
                data["indexed"].append(tmp)

        return data
    # indexed_and_unindexed_on_detector()

    def indexed_and_unindexed_by_frame_on_detector(self):
        if len(self.items) == 0: return

        data = {"indexed": {}, "unindexed": {}}

        for xyz, intensity, iseg, hkl in self.items:
            frame = int(xyz[2])+1
            tmp = xyz[:2]

            key = "unindexed" if hkl[0] is None or all(map(lambda h: h==0, hkl)) else "indexed"
            data[key].setdefault(frame, []).append(tmp)

        return data
    # indexed_and_unindexed_by_frame_on_detector()

    def set_xparm(self, xparm_in):
        xparm = XPARM(xparm_in)
        self.calc_d = lambda x, y: xparm.wavelength/2./math.sin(0.5*math.atan(math.sqrt((x-xparm.origin[0])**2+(y-xparm.origin[1])**2)*xparm.qx/abs(xparm.distance)))
    # set_xparm()

    def set_xdsinp(self, xdsinp):
        inp = dict(get_xdsinp_keyword(xdsinp))
        wavelength = float(inp["X-RAY_WAVELENGTH"])
        orgx, orgy = map(float, (inp["ORGX"], inp["ORGY"]))
        qx = float(inp["QX"])
        distance = abs(float(inp["DETECTOR_DISTANCE"]))

        # XXX no support for non-normal incident beam or multipanel detector
        self.calc_d = lambda x, y: wavelength/2./math.sin(0.5*math.atan(math.sqrt((x-orgx)**2+(y-orgy)**2)*qx/distance))
    # set_xdsinp()


# class SpotXds
