#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
"""
Try several scaling parameters in CORRECT step of XDS, and make mtz files which are ready for refinement.

Usage:
  PHENIX_TRUST_OTHER_ENV=1 phenix.python xds_try_scale_params.py wfac1=1,1.2,0.8,...
"""


master_params_str = """\
wfac1 = [1.0]
    .type = floats
minpk = [75.0]
    .type = floats

mtz_free = None
    .type = path
"""
import os

import iotbx.phil

from yamtbx.dataproc.xds import *
from yamtbx.util import call
from yamtbx.dataproc.xds.command_line import xds2mtz
from yamtbx.dataproc.xds import files

def copy_testflag(mtzfree, mtzin):
    call(cmd="copy_free_R_flag.py -r %s %s" % (mtzfree, mtzin),
         wdir=os.path.dirname(mtzin))
# copy_testflag()

def get_digit(x):
    """
    x must be sorted.
    """
    assert len(x) > 0

    if len(x) == 1:
        m = str(x[0])
    else:
        m = str(min(map(lambda x,y:y-x, x[:-1], x[1:])))
    if "." in m:
        return len(m) - m.index(".") - 1
    else:
        return 0
# get_digit()

def run(params):
    xdsinp = "XDS.INP"
    workdir = os.getcwd()

    params.wfac1.sort()
    params.minpk.sort()
    digit_wfac1 = str(get_digit(params.wfac1))
    digit_minpk = str(get_digit(params.minpk))

    backup_needed = files.generated_by_CORRECT + ("XDS.INP",)

    bk_prefix = make_backup(backup_needed)

    try:
        for wfac1 in params.wfac1:
            for minpk in params.minpk:
                work_name = ("CORRECT_wfac1_%."+digit_wfac1+"f_minpk_%."+digit_minpk+"f")%(wfac1,minpk)
                modify_xdsinp(xdsinp, inp_params=[("JOB","CORRECT"),
                                                  ("WFAC1", wfac1),
                                                  ("MINPK", minpk),
                                                  ])
                call("xds")
                # make .mtz file
                xds2mtz.xds2mtz(os.path.abspath("XDS_ASCII.HKL"),
                                dir_name=os.path.join(workdir, "ccp4_"+work_name),
                                run_ctruncate=True, run_xtriage=True
                                )

                if params.mtz_free is not None:
                    copy_testflag(mtzfree=params.mtz_free,
                                  mtzin=os.path.join(workdir, "ccp4_"+work_name, "XDS_ASCII.mtz"))

                make_backup(backup_needed, work_name+"_")

    finally:
        # 6. Revert XDS.INP, etc.
        revert_files(backup_needed, bk_prefix)

# run()

if __name__ == "__main__":

    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()

    if params.mtz_free is not None:
        params.mtz_free = os.path.abspath(params.mtz_free)

    print "Parameters:"
    cmdline.work.format(python_object=params).show(out=sys.stdout, prefix=" ")

    if len(params.wfac1) == 1 and len(params.minpk) == 1:
        print "Specify params!"
        sys.exit(1)

    run(params)
