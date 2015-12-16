#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

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
"""
from yamtbx.dataproc.xds import modify_xdsinp
from yamtbx.dataproc.xds import files
from yamtbx.dataproc.xds import make_backup
from yamtbx.dataproc.xds import get_xdsinp_keyword
from yamtbx.dataproc.xds import revert_files
from yamtbx.util import call

import re
import iotbx.phil

def analyze_result(idxreflp):
    re_outof = re.compile("^ *([0-9]+) OUT OF *([0-9]+) SPOTS INDEXED.")
    outof = [0,1]

    for l in open(idxreflp):
        r = re_outof.search(l)
        if r:
            outof = int(r.group(1)), int(r.group(2))

    print "  Indexed: %d/%d (%.2f%%)" % (outof[0],outof[1], 100.*float(outof[0])/float(outof[1]))

def run(params):
    xdsinp = "XDS.INP"
    kwds = dict(get_xdsinp_keyword(xdsinp))
    orgx_org, orgy_org = map(float, (kwds["ORGX"], kwds["ORGY"]))

    dx, dy = params.dx, params.dy
    if params.unit == "mm":
        assert "QX" in kwds
        assert "QY" in kwds
        dx /= float(kwds["QX"])
        dy /= float(kwds["QY"])

    backup_needed = files.generated_by_IDXREF + ("XDS.INP",)
    bk_prefix = make_backup(backup_needed)
    try:
        results = []
        for i in xrange(-params.nx, params.nx+1):
            for j in xrange(-params.ny, params.ny+1):
                work_name = "bs_x%+.2d_y%+.2d" % (i, j)
                orgx = orgx_org + i * dx
                orgy = orgy_org + j * dy
                print "Trying", orgx, orgy

                modify_xdsinp(xdsinp, inp_params=[("JOB", "IDXREF"),
                                                  ("ORGX", orgx),
                                                  ("ORGY", orgy),
                                                  ])
                call("xds")
                make_backup(backup_needed, work_name+"_")

                results.append([work_name, orgx, orgy])

        for ret in results:
                print ret,
                analyze_result(ret[0]+"_IDXREF.LP")
                

    finally:
        revert_files(backup_needed, bk_prefix)

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    run(params)
