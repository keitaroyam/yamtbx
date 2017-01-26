"""
(c) RIKEN 2016. All rights reserved. 
Author: Keitaro Yamashita
This software is released under the new BSD License; see LICENSE.
"""
import re
import numpy
import sys

#float_or_1 = lambda x: 1 if x == "" else float(x)
float_or_1 = lambda x: (1,-1)[len(x)] if x == "" or x == "-" else float(x)

class Sensor:
    re_x = re.compile("([-0-9\.]*) *x")
    re_y = re.compile("([-+0-9\.]*) *y")

    def __init__(self, clen=None, res=None):
        for k in ("min_fs", "min_ss", "max_fs", "max_ss", "fs", "ss", "corner_x", "corner_y", "coffset"):
            setattr(self, k, None)
            
        self.clen = clen
        self.res = res
    # __init__()

    def set_info(self, key, valstr):
        if key == "min_fs":
            self.min_fs = int(float(valstr))
        elif key == "min_ss":
            self.min_ss = int(float(valstr))
        elif key == "max_fs":
            self.max_fs = int(float(valstr))
        elif key == "max_ss":
            self.max_ss = int(float(valstr))
        elif key == "fs":
            self.fs = [0., 0.]
            r_x = self.re_x.search(valstr)
            r_y = self.re_y.search(valstr)
            if r_x: self.fs[0] = float_or_1(r_x.group(1))
            if r_y: self.fs[1] = float_or_1(r_y.group(1))
            if not r_x and not r_y: print "Wrong parameter value for fs:", valstr
        elif key == "ss":
            self.ss = [0., 0.]
            r_x = self.re_x.search(valstr)
            r_y = self.re_y.search(valstr)
            if r_x: self.ss[0] = float_or_1(r_x.group(1))
            if r_y: self.ss[1] = float_or_1(r_y.group(1))
            if not r_x and not r_y: print >>sys.stderr, "Wrong parameter value for ss:", valstr
        elif key == "corner_x":
            self.corner_x = float(valstr)
        elif key == "corner_y":
            self.corner_y = float(valstr)
        elif key == "res":
            self.res = float(valstr)
        elif key == "clen":
            self.clen = float(valstr)
        elif key == "coffset":
            self.coffset = float(valstr)
        else:
            print >>sys.stderr, "Warning - Unknown paramter:", key
    # read_info()
# class Sensor

class Geomfile:
    def __init__(self, geom_in=None):
        self.clen = None
        self.res = None
        if geom_in is not None:
            self.read_file(geom_in)
    # __init__()

    def read_file(self, geom_in):
        self.sensors = {}

        for l in open(geom_in):
            if ";" in l: l = l[:l.find(";")]
            l = l.strip()

            if "=" in l:
                sp = map(lambda x:x.strip(), l.split("="))
                if len(sp) > 2:
                    print "Warning: more than one = in this line. Ignored:", l
                    continue

                lhs, rhs = sp
                if lhs == "clen":
                    self.clen = rhs #float
                elif lhs == "res":
                    self.res = float(rhs)
                elif lhs == "adu_per_eV":
                    self.adu_per_eV = rhs #float
                elif lhs == "max_adu":
                    self.max_adu = rhs # float
                elif "/" in lhs:
                    sp = lhs.split("/")
                    if len(sp) > 2:
                        print >>sys.stderr, "Warning: more than one / in left hand. Ignored:", l
                        continue

                    sid, lhs2 = sp
                    if sid not in self.sensors:
                        self.sensors[sid] = Sensor(clen=self.clen, res=self.res) # XXX res and clen must be defined first (before the sensor defs; if sensor defs don't have clen nor res)

                    self.sensors[sid].set_info(lhs2, rhs)
    # read_file()

    def get_extent(self):
        inf = float("inf")
        fmin, fmax, smin, smax = inf, -inf, inf, -inf
        for sid in self.sensors:
            s = self.sensors[sid]
            c = numpy.array((s.corner_x, s.corner_y))
            fs, ss = numpy.array(s.fs), numpy.array(s.ss)
            #for x in xrange(s.max_fs-s.min_fs+1):
            #    for y in xrange(s.max_ss-s.min_ss+1):
            for x in (0, s.max_fs-s.min_fs):
                for y in (0, s.max_ss-s.min_ss):
                    v = c + fs * x + ss * y
                    fmin = min(fmin, v[0])
                    fmax = max(fmax, v[0])
                    smin = min(smin, v[1])
                    smax = max(smax, v[1])

        return fmin, fmax, smin, smax
    # get_extent()

    def get_image_extent(self):
        max_ss, max_fs = 0, 0
        for sid in self.sensors:
            s = self.sensors[sid]
            max_fs = max(max_fs, s.max_fs)
            max_ss = max(max_ss, s.max_ss)
        return max_fs, max_ss
    # get_image_extent()

    def calc_xy_from_fs_ss(self, fs, ss, panel=None, shifts=None):
        if panel is None:
            for sid in self.sensors:
                s = self.sensors[sid]
                if s.min_fs <= fs <= s.max_fs and s.min_ss <= ss <= s.max_ss:
                    panel = sid
                    break
        
        if panel is None:
            return float("nan"), float("nan")

        s = self.sensors[panel]
        fs, ss = fs - s.min_fs, ss - s.min_ss
                
        x = s.corner_x + s.fs[0]*fs + s.ss[0]*ss
        y = s.corner_y + s.fs[1]*fs + s.ss[1]*ss

        if shifts is not None: # shifts in meter
            x += shifts[0]*s.res
            y += shifts[1]*s.res
            
        return x,y
    # calc_xy_from_fs_ss()
# class Geomfile
