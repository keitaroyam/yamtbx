"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import numpy
from cctbx.array_family import flex
from cctbx import miller
from cctbx import crystal

class DenzoXfile:
    """
    Reference:
    http://www.hkl-xray.com/sites/default/files/manual_online.pdf
    p. 44
    """
    def __init__(self, filein, *args, **kwargs):
        if filein is not None:
            self.read_file(filein, *args, **kwargs)
    # __init__()
    
    def read_file(self, filein, onlyfull=False, onlypartial=False, read_I_sum=False):
        assert (onlyfull,onlypartial).count(True) < 2
        
        ifs = open(filein)
        self.comment = ifs.readline()
        if self.comment.startswith(" 999 ") and len(self.comment.split()) == 3: raise RuntimeError("%s does not contains useful info" % filein)
        self.amatrix, self.umatrix = numpy.zeros((3,3)), numpy.zeros((3,3))
        for i in xrange(3):
            vals = map(float, ifs.readline().strip().split())
            self.amatrix[i,] = vals[0:3]
            self.umatrix[i,] = vals[3:6]

        #  osc. start osc. end dist. in pixel units wavelength cryst rotz roty rotx mosaicity
        l = ifs.readline()
        tmp = l[:12], l[12:24], l[24:36], l[36:48], l[48:58], l[58:68], l[68:78], l[78:88]
        self.osc_start, self.osc_end, self.distance, self.wavelength, self.rotz,self.roty,self.rotx, self.mosaicity = map(float, tmp)

        # these are the actual reflections now....
        # 1,2,3 h,k,l not reduced 4 flag 0 full, 1 partial 5 I by profile fitting 6 I by profile summation 7 chi2
        # of profile
        # fitting 8 sigma of I 9 cosine of incidence angle at detector 10 x coordinate of predicted centroid, in pixel units 11 y
        # coordinate of predicted centroid, in pixel units 12 Lorentz, polarization, obliquity factor 13 strength of
        # averaged profile, arbitrary units
        # In columns 5,6, and 8 the intensity of very strong reflections (i.e. overflows) will be written as an integer,
        # rather than as a floating point number - see example on line T above. Scalepack understands this Denzo
        # convention.

        self.miller_indices = flex.miller_index()
        self.data = flex.double()
        self.sigmas = flex.double()
        
        for l in ifs:
            if l.startswith(" 999 "):
                boxsize = map(int, l.strip().split())[1:]
                break
            
            #sp = l.strip().split()
            sp = l[:4], l[4:8], l[8:12], l[13], l[14:22], l[22:30], l[30:37], l[37:43], l[43:49], l[49:56], l[56:63], l[63:69], l[69:77]
            h, k, l, isfull = map(int, sp[:4])
            Ipf, Ips, chisq, sigma, cos, calx, caly, lpo, strength = map(float, sp[4:])

            if onlyfull and not isfull: continue
            if onlypartial and isfull: continue
            
            self.miller_indices.append((h,k,l))
            self.data.append(Ips if read_I_sum else Ipf)
            self.sigmas.append(sigma)

        self.intbox = [] # 1: background area, 0: guard area, 2: spot area
        for l in ifs:
            if l.startswith(" "):
                self.intbox.append(map(int, l.strip()))
            elif l.startswith("unit cell"):
                self.unit_cell = map(float, l[9:].strip().split())
            elif l.startswith("space group"):
                self.space_group = l[12:].strip()
    # read_file()
    
    def miller_set(self, anomalous_flag=True):
        return miller.set(crystal_symmetry=crystal.symmetry(unit_cell=self.unit_cell, space_group_symbol=self.space_group),
                          indices=self.miller_indices, anomalous_flag=anomalous_flag)
    # miller_set()
            
    def miller_array(self, anomalous_flag=True):
        return miller.array(miller_set=self.miller_set(anomalous_flag),
                            data=self.data, sigmas=self.sigmas)
    # miller_array()
# class DenzoXfile

if __name__ == "__main__":
    import sys
    x = DenzoXfile(sys.argv[1])
    a = x.miller_array()
    a.show_summary()

    print "A="
    print x.amatrix
    print "U="
    print x.umatrix
