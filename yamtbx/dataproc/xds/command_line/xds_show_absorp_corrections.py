#!/usr/bin/env yamtbx.python

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
master_params_str = """
bkgpix_cbf = "BKGPIX.cbf"
 .type = path
absorp_cbf = "ABSORP.cbf"
 .type = path
correct_lp = "CORRECT.LP"
 .type = path

output {
 prefix = "absorp"
  .type = str
}
"""

import re
import os
import numpy
import iotbx.phil
from yamtbx.dataproc import cbf

def read_correct_lp(lpin):
    re_x = re.compile("XMIN= *([0-9\.]+) XMAX= *([0-9\.]+) NXBIN= *([0-9]+)")
    ret = [None, None, None, []] # [xmin, xmax, nxbin, [positions,..]]
    read_flag = False

    for l in open(lpin):
        if "CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf" in l:
            read_flag = True
        elif read_flag and "XMIN=" in l:
            r = re_x.search(l)
            xmin, xmax, nxbin = r.groups()
            ret[:3] = float(xmin), float(xmax), int(nxbin)
        elif read_flag and "DETECTOR_SURFACE_POSITION=" in l:
            ret[3].append(list(map(int, l[l.index("=")+1:].split())))
        elif "CORRECTION PARAMETERS FOR THE STANDARD ERROR OF REFLECTION INTENSITIES" in l:
            break
    return ret
# read_correct_lp()

def read_absorp_cbf(cbfin):
    data, ndimfast, ndimmid = cbf.load_minicbf_as_numpy(cbfin)
    data = data.reshape(ndimmid, ndimfast)
    #fac = data[ipos, ix]
    return data
# read_absorp_cbf()

def read_bkgpix_cbf(cbfin, positions):
    data, nx, ny = cbf.load_minicbf_as_numpy(cbfin)
    data = data.reshape(ny, nx)
    npos = len(positions)

    dists = numpy.array((data,)*npos)

    gx, gy = numpy.meshgrid(range(nx), range(ny))

    for i, pos in enumerate(positions):
        print(i, pos)
        dists[i,] = (gx-pos[0])**2 + (gy-pos[1])**2 # squared distance to each position
        
    ret = dists.argmin(axis=0) + 1
    ret[data < 0] = -10
    return ret
# read_bkgpix_cbf()

def run(params):
    xmin, xmax, nxbin, positions = read_correct_lp(params.correct_lp)
    xstep = (xmax-xmin)/nxbin

    absorp = read_absorp_cbf(params.absorp_cbf)
    assert numpy.all(absorp > 13) # if 

    det = read_bkgpix_cbf(params.bkgpix_cbf, positions)

    datout = open("%s.dat"%params.output.prefix, "w")
    datout.write("ix xmin xmax ipos posx posy fac\n")

    for ix in range(nxbin):
        x1, x2 = xmin + ix*xstep, xmin + (ix+1)*xstep
        tmp = numpy.array(det)
        for i in range(len(positions)):
            tmp[tmp==i+1] = absorp[i,ix]
            datout.write("%2d %.2f %.2f %2d %d %d %5d\n" % (ix, x1, x2, i, positions[i][0], positions[i][1], absorp[i,ix]))

        tmp = numpy.array(tmp, dtype=numpy.int32)
        cbf.save_numpy_data_as_cbf(tmp.flatten(), tmp.shape[1], tmp.shape[0],
                                   "absorp_%.1f-%.1f"% (x1,x2),
                                   "%s_%.3d.cbf" % (params.output.prefix, ix))

    datout.close()
# run()

if __name__ == "__main__":
    import sys

    if "-h" in sys.argv or "--help" in sys.argv:
        print("""\
This script generates cbf files for visual inspection of absorption correction factors, which are recorded in ABSORP.cbf.
In addition, a dat file is generated for plotting purpose etc.

Parameters:
        """)
        iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
        quit()

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    run(params)
