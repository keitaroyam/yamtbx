#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

master_params_str = """\
lstin = None
 .type = path
dmin_lst = None
 .type = path
anomalous = False
 .type = bool
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
nproc = None
 .type = int

wfac1 = None
 .type = float
reference = bmax bmed bmin
 .type = choice(multi=False)
 .help = Reference selection for final scaling
frames_per_batch = None
 .type = int
 .help = affects NBATCH=. When 1, NBATCH= matches the number of frames in the dataset.
corrections = decay modulation absorption
 .type = choice(multi=True)
"""

import os
import sys
import math
import iotbx.phil
from yamtbx import util
from yamtbx.dataproc.xds import xscale
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII

xscale_comm = "xscale_par"

def check_valid_xac(xac):
    if not os.path.isfile(xac): return False

    line = open(xac).readline()
    return "FORMAT=XDS_ASCII" in line
# check_valid_xac()

def get_xac_info(xac, get_nframes=False):
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

    if get_nframes:
        frame_range = XDS_ASCII(f, read_data=False).get_frame_range()
        ret["nframes"] = frame_range[1] - frame_range[0]

    return ret
# get_xac_info()

def make_shells(d_max, d_min, nbins):
    step = ( 1./(d_min**2) - 1./(d_max**2) ) / float(nbins)
    start = 1./(d_max**2)
    rshells = " ".join(map(lambda i: "%.2f" % (start + i * step)**(-1./2), xrange(1, nbins+1)))
    
    return " RESOLUTION_SHELLS= %s\n" % rshells
# make_shells()

def prep_xscale_inp(wdir, xscale_inp_head, xac_files, infos, frames_per_batch=None, corrections=None, ref_idx=None):
    xscale_inp = os.path.join(wdir, "XSCALE.INP")
    inp_out = open(xscale_inp, "w")
    inp_out.write(xscale_inp_head)

    for i, xds_ascii in enumerate(xac_files):
        star = "*" if ref_idx == i else " "
        inp_out.write("  INPUT_FILE=%s%s\n" % (star, os.path.relpath(xds_ascii, params.workdir)))

        if "resol_range_user" in infos[xds_ascii]:
            inp_out.write("  INCLUDE_RESOLUTION_RANGE= %s\n" % infos[xds_ascii]["resol_range_user"])
        else:
            inp_out.write("  ! INCLUDE_RESOLUTION_RANGE= %s\n" % infos[xds_ascii]["resol_range"])

        if frames_per_batch:
            nbatch = int(math.ceil(infos[xds_ascii]["nframes"] / frames_per_batch))
            inp_out.write("    NBATCH= %d\n" % nbatch)
        if corrections:
            inp_out.write("    CORRECTIONS= %s\n" % " ".join(map(lambda s: s.upper(), corrections)))

        inp_out.write("\n")

    inp_out.close()
# prep_xscale_inp()

def run(params, xac_files):
    if len(xac_files) == 0:
        print "No XDS_ASCII.HKL files provided."
        return

    # Parse
    dmin_dict = {}
    if params.dmin_lst:
        for l in open(params.dmin_lst):
            sp = l.split()
            if len(sp) != 2: continue
            f, dmin = sp
            dmin_dict[f] = dmin

    xscale_inp_head = "!MINIMUM_I/SIGMA= 3\n\n" 
    if params.wfac1 is not None:
        xscale_inp_head += "WFAC1= %.3f\n" % params.wfac1

    if params.nproc:
        xscale_inp_head += "MAXIMUM_NUMBER_OF_PROCESSORS= %d\n" % params.nproc

    infos = {}
    d_max, d_min = 0, 100
    cells = []
    for xds_ascii in xac_files:
        info = get_xac_info(xds_ascii, get_nframes=params.frames_per_batch is not None)
        if xds_ascii in dmin_dict:
            dmax,dmin = info["resol_range"].split()
            info["resol_range_user"] = "%s %s" % (dmax, dmin_dict[xds_ascii])

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
  
    xscale_inp_head += make_shells(d_max, d_min, params.nbins) + "\n"
    xscale_inp_head += " OUTPUT_FILE= %s\n" % params.output
    xscale_inp_head += "  FRIEDEL'S_LAW= %s\n\n" % ("FALSE" if params.anomalous else "TRUE")

    prep_xscale_inp(params.workdir, xscale_inp_head, xac_files, infos, params.frames_per_batch, params.corrections)

    xscale.run_xscale(os.path.join(params.workdir, "XSCALE.INP"),
                      cbf_to_dat=params.cbf_to_dat,
                      use_tmpdir_if_available=params.use_tmpdir_if_available)

    if params.reference:
        print "Choosing reference data (reference=%s)" % params.reference
        ref_idx = xscale.decide_scaling_reference_based_on_bfactor(os.path.join(params.workdir, "XSCALE.LP"), params.reference, return_as="index")
        if ref_idx != 0:
            for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(params.workdir, f))
            prep_xscale_inp(params.workdir, xscale_inp_head, xac_files, infos, params.frames_per_batch, params.corrections, ref_idx=ref_idx)
            xscale.run_xscale(os.path.join(params.workdir, "XSCALE.INP"),
                              cbf_to_dat=params.cbf_to_dat,
                              use_tmpdir_if_available=params.use_tmpdir_if_available)

# run()

if __name__ == "__main__":
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    xac_files = filter(check_valid_xac, cmdline.remaining_args)
    if params.lstin:
        xac_files.extend(filter(check_valid_xac, util.read_path_list(params.lstin)))

    print "XDS_ASCII.HKL files given:"
    for f in xac_files:
        print " %s" % f
    print

    run(params, xac_files)
