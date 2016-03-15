"""
Ver. 2.0 format looks like:

CrystFEL reflection list version 2.0
Symmetry: 422
   h    k    l          I    phase   sigma(I)   nmeas
   0    0    1      -9.83        -      12.58       5
   0    0    2      -3.89        -       8.33     376
...
"""

from iotbx import crystal_symmetry_from_any
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from libtbx.utils import Sorry

class HKLfile:
    def __init__(self, symm_source=None, hklin=None):
        self.symm = None

        if symm_source is not None:
            self.symm = crystal_symmetry_from_any.extract_from(symm_source)
            #print self.symm.space_group().is_chiral()
            print "Symmetry from",symm_source
            print self.symm.show_summary(prefix=" ")

        if hklin is not None:
            self.read_file(hklin)
    # __init__()

    def read_file(self, hklin):
        assert self.symm is not None # XXX more careful check
        fin = open(hklin)
        line = fin.readline() # first line
        
        #TODO support other format.
        assert "CrystFEL reflection list version 2.0" in line

        line = fin.readline() # second line
        sym_str = line[line.find(":")+1:].strip()
        print "HKL pointgroup:",sym_str,

        anomalous_flag = None
        if sym_str == self.symm.space_group().laue_group_type():
            anomalous_flag = False
        elif sym_str == self.symm.space_group().point_group_type():
            anomalous_flag = True
        else:
            Sorry("Incompatible symmetry: %s and %s" % (sym_str, self.symm.space_group().info()))
        print "Anomalous:", anomalous_flag

        line = fin.readline() # third line
        assert "   h    k    l          I    phase   sigma(I)   nmeas" in line

        indices, Is, sigIs, ns = [], [], [], []
        for l in fin:
            if l.startswith("End of reflections"):
                break

            h,k,l,i,phase,sigi,nmeas = l.strip().split()
            h,k,l,nmeas = map(int, (h,k,l,nmeas))
            i,sigi = map(float, (i,sigi))
            indices.append((h,k,l))
            Is.append(i)
            sigIs.append(sigi)
            ns.append(nmeas)
        
        miller_set = miller.set(crystal_symmetry=self.symm,
                                indices=flex.miller_index(indices),
                                anomalous_flag=anomalous_flag)

        self.array = miller.array(miller_set=miller_set,
                                  data=flex.double(Is),
                                  sigmas=flex.double(sigIs)
                                  ).set_observation_type_xray_intensity()

        self.redundancies = miller.array(miller_set=miller_set,
                                  data=flex.int(ns))
    # read_file()

    def set_resolution(self, d_min=None, d_max=None):
        self.array = self.array.resolution_filter(d_min=d_min, d_max=d_max)
        self.redundancies = self.redundancies.resolution_filter(d_min=d_min, d_max=d_max)
    # set_resolution()
# class HKLfile

def write_file(filename, i_obs, nmeas):
    if i_obs.anomalous_flag():
        symm_str = i_obs.space_group().point_group_type()
    else:
        symm_str = i_obs.space_group().build_derived_laue_group().point_group_type()

    ofs = open(filename, "w")
    ofs.write("CrystFEL reflection list version 2.0\n")
    ofs.write("Symmetry: %s\n" % symm_str)
    ofs.write("   h    k    l          I    phase   sigma(I)   nmeas\n")
    
    assert all(i_obs.indices() == nmeas.indices())

    for hkl, i, sigma, n in zip(i_obs.indices(), i_obs.data(), i_obs.sigmas(), nmeas.data()):
        if n == 0: continue
        ofs.write("%4i %4i %4i %10.2f %s %10.2f %7i\n" % (hkl[0], hkl[1], hkl[2], i, "       -", sigma, n))
    
    ofs.write("End of reflections\n")
# write_file()

if __name__ == "__main__":
    import sys

    HKLfile(symm_source=sys.argv[1],
            hklin=sys.argv[2])
