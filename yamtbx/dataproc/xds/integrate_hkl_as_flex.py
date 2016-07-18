"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import re
import os
from cctbx.array_family import flex
from cctbx import uctbx
from cctbx import crystal
from cctbx import miller

def is_integrate_hkl(filein):
    if not os.path.isfile(filein): return False

    line = open(filein).readline()
    return "!OUTPUT_FILE=INTEGRATE.HKL" in line
# is_xds_ascii()

class reader:
    def __init__(self, filein, read_columns, read_data=True):
        assert len(set(["H","K","L"]).intersection(set(read_columns))) == 0

        read_data_flag = False
        re_column_info = re.compile("[A-Z],[A-Z]")
        column_names = []
        read_indices = []

        self.hkl = flex.miller_index()
        self.data = {}
        
        for l in open(filein):

            if l.startswith("!UNIT_CELL_CONSTANTS="):
                cell = map(lambda x:float(x), l[len("!UNIT_CELL_CONSTANTS="):].strip().split())
                self.unit_cell = uctbx.unit_cell(cell)
            elif l.startswith("!SPACE_GROUP_NUMBER="):
                self.space_group_number = int(l[len("!SPACE_GROUP_NUMBER="):].strip())
                assert 1 <= self.space_group_number <= 230
            elif l.startswith("!X-RAY_WAVELENGTH="):
                self.wavelength = float(l[len("!X-RAY_WAVELENGTH="):].strip())
            elif "DETECTOR_DISTANCE=" in l:
                self.distance = float(l[l.index("DETECTOR_DISTANCE=")+len("DETECTOR_DISTANCE="):].strip())
            elif l.startswith("!UNIT_CELL_A-AXIS="):
                self.a_axis = map(float, l[len("!UNIT_CELL_A-AXIS="):].strip().split())
            elif l.startswith("!UNIT_CELL_B-AXIS="):
                self.b_axis = map(float, l[len("!UNIT_CELL_B-AXIS="):].strip().split())
            elif l.startswith("!UNIT_CELL_C-AXIS="):
                self.c_axis = map(float, l[len("!UNIT_CELL_C-AXIS="):].strip().split())


            if re_column_info.search(l):
                l = l[1:].strip() # remove !
                column_names.extend(filter(lambda x:x!="", l.split(",")))

            if l.startswith("!END_OF_HEADER"):
                if not read_data:
                    break

                assert set(read_columns).issubset(set(column_names))

                read_data_flag = True
                read_indices = [i for i, c in enumerate(column_names) if c in read_columns]
                # Set up data structure
                for i in read_indices:
                    self.data[column_names[i]] = flex.double()
                continue

            if read_data_flag:
                if l.startswith("!END_OF_DATA"):
                    break

                sp = l.split()
                self.hkl.append(map(lambda x:int(x), sp[:3]))
                for i in xrange(3, len(sp)):
                    if i in read_indices:
                        self.data[column_names[i]].append(float(sp[i]))


        #print column_names, read_indices
        #print self.data

        self.column_names = column_names
        self._filein = filein
    # __init__ ()

    def get_column_names(self): return self.column_names

    def crystal_symmetry(self):
        return crystal.symmetry(unit_cell=self.unit_cell,
                                space_group_symbol=self.space_group_number)

    def arrays(self):
        ret = {}

        for key in self.data:
            arr = miller.array(miller_set=miller.set(crystal_symmetry=self.crystal_symmetry(),
                                                     indices=self.hkl,
                                                     anomalous_flag=False),
                               data=self.data[key])
            ret[key] = arr

        return ret
    # arrays()

    def i_obs(self, anomalous_flag=None):
        assert "IOBS" in self.data
        assert "SIGMA" in self.data

        array_info = miller.array_info(source_type="xds_integrate")#, wavelength=)
        return miller.array(miller_set=miller.set(crystal_symmetry=self.crystal_symmetry(),
                                                  indices=self.hkl,
                                                  anomalous_flag=anomalous_flag),
                            data=self.data["IOBS"],
                            sigmas=self.data["SIGMA"]).set_info(array_info).set_observation_type_xray_intensity()
    # i_obs()

    def write_peak_corrected(self, hklout):
        ofs = open(hklout, "w")

        i_peak = self.column_names.index("PEAK")
        i_iobs, i_sigma = self.column_names.index("IOBS"), self.column_names.index("SIGMA")

        data_flag = False

        for line in open(self._filein):
            if line.startswith('!END_OF_HEADER'):
                ofs.write(line)
                data_flag = True
            elif line.startswith("!END_OF_DATA"):
                ofs.write(line)
                break
            elif not data_flag:
                ofs.write(line)
            elif data_flag:
                sp = line.split()
                peak = float(sp[i_peak]) * .01
                sp[i_iobs] = "%.3e" % (float(sp[i_iobs]) * peak)
                sp[i_sigma] = "%.3e" % (float(sp[i_sigma]) * peak)
                ofs.write(" %s\n" % " ".join(sp))
    # write_peak_corrected()


if __name__ == "__main__":
    import sys

    reader(sys.argv[1], ["IOBS","SIGMA","PEAK"])
