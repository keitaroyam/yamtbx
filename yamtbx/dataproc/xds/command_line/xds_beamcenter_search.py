#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

master_params_str = """\
dx = 2
 .type = float
dy = 2
 .type = float
nx = 3
 .type = int
ny = 3
 .type = int
unit = *px mm
 .type = choice(multi=False)
workdir = orgxy_search
 .type = str
nproc = 1
 .type = int
"""
from yamtbx.dataproc.xds import modify_xdsinp
from yamtbx.dataproc.xds import files
from yamtbx.dataproc.xds import make_backup
from yamtbx.dataproc.xds import get_xdsinp_keyword
from yamtbx.dataproc.xds import revert_files
from yamtbx.util import call

import re
import iotbx.phil
from libtbx import easy_mp
import shutil
import os

def analyze_result(idxreflp):
    re_outof = re.compile("^ *([0-9]+) OUT OF *([0-9]+) SPOTS INDEXED.")
    outof = [0,1]

    for l in open(idxreflp):
        r = re_outof.search(l)
        if r:
            outof = int(r.group(1)), int(r.group(2))

    print("  Indexed: %d/%d (%.2f%%)" % (outof[0],outof[1], 100.*float(outof[0])/float(outof[1])))

def work(rootdir, xdsinp, orgxy):
    workdir = os.path.join(rootdir, "bs_x%+.1f_y%+.1f" % orgxy)
    inpdir = os.path.normpath(os.path.dirname(xdsinp))
    print(workdir)
    if os.path.exists(workdir):
        print("- exists. skipping.")
        return
    
    os.makedirs(workdir)
    shutil.copyfile(os.path.join(inpdir, "SPOT.XDS"), os.path.join(workdir, "SPOT.XDS"))
    shutil.copyfile(xdsinp, os.path.join(workdir, "XDS.INP"))
    modify_xdsinp(os.path.join(workdir, "XDS.INP"), inp_params=[("JOB", "IDXREF"),
                                                                ("ORGX", orgxy[0]),
                                                                ("ORGY", orgxy[1]),
                                                                ])
    call("xds", wdir=workdir)
# work()

def run(params):
    xdsinp = "XDS.INP"
    kwds = dict(get_xdsinp_keyword(xdsinp))
    orgx_org, orgy_org = list(map(float, (kwds["ORGX"], kwds["ORGY"])))

    dx, dy = params.dx, params.dy
    if params.unit == "mm":
        assert "QX" in kwds
        assert "QY" in kwds
        dx /= float(kwds["QX"])
        dy /= float(kwds["QY"])

    #backup_needed = files.generated_by_IDXREF + ("XDS.INP",)
    #bk_prefix = make_backup(backup_needed)

    orgxy_list = []
    for i in range(-params.nx, params.nx+1):
        for j in range(-params.ny, params.ny+1):
            orgxy_list.append((orgx_org + i * dx, orgy_org + j * dy))

    easy_mp.pool_map(fixed_func=lambda x: work(os.path.abspath(params.workdir), os.path.abspath(xdsinp), x),
                     args=orgxy_list,
                     processes=params.nproc)
    #for ret in results:
    #    print ret,
    #    analyze_result(ret[0]+"_IDXREF.LP")

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    run(params)
