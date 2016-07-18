#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

master_params_str = """\
cell = average *first
 .type = choice(multi=False)
sgnum = None
 .type = int
nbins = 9
 .type = int
output = "XSCALE.HKL"
 .type = str
d_min = None
 .type = float
workdir = "."
 .type = path
use_tmpdir_if_available = False
 .type = bool
cbf_to_dat = True
 .type = bool
"""

import os
import sys
import iotbx.phil
from yamtbx import util
from yamtbx.dataproc.xds import xscale

xscale_comm = "xscale_par"

def check_valid_xac(xac):
    if not os.path.isfile(xac): return False

    line = open(xac).readline()
    return "FORMAT=XDS_ASCII" in line
# check_valid_xac()

def get_xac_info(xac):
    ret = {}

    for l in open(xac):
        if l.startswith("!FORMAT=XDS_ASCII"): # !FORMAT=XDS_ASCII    MERGE=FALSE    FRIEDEL'S_LAW=FALSE
            ret["friedels_law"] = l[l.rindex("=")+1:].strip()
        if l.startswith("!INCLUDE_RESOLUTION_RANGE="):
            ret["resol_range"] = l[l.index("=")+1:].strip()
        elif l.startswith("!SPACE_GROUP_NUMBER="):
            ret["spgr_num"] = l[l.index("=")+1:].strip()
        elif l.startswith("!UNIT_CELL_CONSTANTS="):
            ret["cell"] = l[l.index("=")+1:].strip()
        elif l.startswith("!END_OF_HEADER"):
            break

    return ret
# get_xac_info()

def make_shells(d_max, d_min, nbins):
    step = ( 1./(d_min**2) - 1./(d_max**2) ) / float(nbins)
    start = 1./(d_max**2)
    rshells = " ".join(map(lambda i: "%.2f" % (start + i * step)**(-1./2), xrange(1, nbins+1)))
    
    return " RESOLUTION_SHELLS= %s\n" % rshells
# make_shells()

def run(params, xac_files):
    if len(xac_files) == 0:
        print "No XDS_ASCII.HKL files provided."
        return

    xscale_inp_head = "!MINIMUM_I/SIGMA= 3\n\n" 

    infos = {}
    d_max, d_min = 0, 100
    cells = []
    for xds_ascii in xac_files:
        info = get_xac_info(xds_ascii)
        infos[xds_ascii] = info

        resrng = map(float, info["resol_range"].split())
        d_max = max(d_max, resrng[0])
        d_min = min(d_min, resrng[1])
        cells.append(map(float, info["cell"].split()))

    if params.d_min is not None:
        d_min = max(params.d_min, d_min)

    if params.cell == "average":
        cell_sum = reduce(lambda x,y: map(lambda a: a[0]+a[1], zip(x,y)), cells)
        cell_mean = map(lambda x: x/float(len(cells)), cell_sum)

        if params.sgnum is not None: sgnum = str(params.sgnum)
        else: sgnum = infos[xac_files[0]]["spgr_num"]

        xscale_inp_head += " SPACE_GROUP_NUMBER= %s\n" % sgnum
        xscale_inp_head += " UNIT_CELL_CONSTANTS= %s\n" % " ".join(map(lambda x: "%.3f"%x, cell_mean))

    #if anomalous_flag is not None:
    #    xscale_inp_head += " FRIEDEL'S_LAW= %s\n" % ("FALSE" if anomalous_flag else "TRUE")
    
    xscale_inp_head += make_shells(d_max, d_min, params.nbins) + "\n"
    xscale_inp_head += " OUTPUT_FILE= %s\n\n" % params.output

    xscale_inp = os.path.join(params.workdir, "XSCALE.INP")
    inp_out = open(xscale_inp, "w")
    inp_out.write(xscale_inp_head)

    for xds_ascii in xac_files:
        inp_out.write("  INPUT_FILE= %s\n" % os.path.relpath(xds_ascii, params.workdir))
        inp_out.write("  ! INCLUDE_RESOLUTION_RANGE= %s\n\n" % infos[xds_ascii]["resol_range"])

    inp_out.close()

    xscale.run_xscale(os.path.join(params.workdir, "XSCALE.INP"),
                      cbf_to_dat=params.cbf_to_dat,
                      use_tmpdir_if_available=params.use_tmpdir_if_available)

# run()

if __name__ == "__main__":
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    xac_files = filter(check_valid_xac, cmdline.remaining_args)

    print "XDS_ASCII.HKL files given:"
    for f in xac_files:
        print " %s" % f
    print

    run(params, xac_files)
