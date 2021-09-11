#!/usr/bin/env cctbx.python
from __future__ import print_function
from __future__ import unicode_literals
import iotbx.phil
from cctbx import uctbx
import os

master_params_str = """
scain = None
 .type = path
scaout = None
 .type = path
d_min = None
 .type = float
d_max = None
 .type = float
"""

def run(params):
    if params.scain == params.scaout: # need to be more careful
        print("MUST BE DIFFERENT FILE!")
        return

    file_handle = open(params.scain)
    ofs = open(params.scaout, "w")

    line = file_handle.readline()
    ofs.write(line)
    if line.rstrip() != "    1":
        raise Exception("line 1: expecting '    1'")

    ofs.write(file_handle.readline()) # ignore line 2
    line_error = "line 3: expecting unit cell parameters and space group label"
    line = file_handle.readline()
    ofs.write(line)
    if (len(line) < 63 or line[60] != ' '):
        raise Exception("Cell reading error")

    uc_params = [float(line[i * 10 : (i + 1) * 10]) for i in range(6)]
    unit_cell = uctbx.unit_cell(uc_params)
    rejected = 0

    for line in file_handle:
        h, k, l = list(map(int, (line[:4], line[4:8], line[8:12])))
        d = unit_cell.d((h,k,l))
        accept = True
        if params.d_max is not None: accept &= (d <= params.d_max)
        if params.d_min is not None: accept &= (d >= params.d_min)
        if accept:
            ofs.write(line)
        else:
            rejected += 1

    print("%d reflections rejected." % rejected)
    print("See", params.scaout)
# run()

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    for arg in args:
        if os.path.isfile(arg) and params.scain is None:
            params.scain = arg

    if params.scain is None:
        print("Give scain=")
        quit()

    if params.scaout is None:
        params.scaout = os.path.splitext(os.path.basename(params.scain))[0] + "_cut.sca"

    run(params)
