#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

import sys, os, optparse, shutil, re, subprocess, random, glob, string, numpy
from yamtbx.dataproc.xds.make_adx import make_adxfile
from yamtbx.dataproc.xds.xparm import get_xparm_from_integrate_lp, XPARM
from yamtbx.dataproc.xds import *
from yamtbx.dataproc.xds import files as xds_files
##
# Helper script for XDS to see predictions of given frame.
#
#
# Scheme:
# 1. Backup XDS.INP, INTEGRATE.LP, INTEGRATE.HKL, FRAME.cbf
# 2. Modify XDS.INP
#    - JOB= INTEGRATE
#    - DATA_RANGE= 
#    - MINPK= 0
# 3. Run xds
# 4. Rename FRAME_%d.cbf, INTEGRATE_%d.{HKL,LP}
# 5. Make .adx file
# 6. Revert XDS.INP, INTEGRATE.LP, INTEGRATE.HKL, FRAME.cbf
# 
# * What you get are: FRAME_%d.cbf, INTEGRATE_%d.{HKL, LP, adx}
# * You can see the prediction using adxv FRAME_%d.cbf INTEGRATE_%d.adx
#

backup_needed = xds_files.generated_by_INTEGRATE + ("XPARM.XDS", "XDS.INP")

def check_needed_files(needed_files, wdir):
    not_exists = [ f for f in needed_files if not os.path.isfile(os.path.join(wdir, f)) ]

    if not_exists != []:
        print("We need these files!")
        print("  " + ",".join(not_exists))
        print()

    return len(not_exists) == 0
# check_needed_files()

def run(param_source, frame_num, wdir, sigmar=None, sigmab=None, need_adx=True, out_prefix=""):
    xparm_str = None

    # Get XPARM.XDS strings from XPARM.XDS, GXPARM.XDS, or INTEGRATE.LP (user-specified)
    try:
        XPARM(param_source)
        xparm_str = open(param_source).read()
    except:
        xparm_str = get_xparm_from_integrate_lp(param_source, frame_num)

    # Check all needed files exist
    if not check_needed_files(xds_files.needed_by_INTEGRATE+("XDS.INP",), wdir):
        return

    xdsinp = os.path.join(wdir, "XDS.INP")

    # 1. Backup XDS.INP, etc. (Make copies; not renaming)
    bk_prefix = make_backup(backup_needed, wdir=wdir)
    for f in xds_files.generated_by_INTEGRATE:
        if os.path.isfile(os.path.join(wdir, f)):
            os.remove(os.path.join(wdir, f))

    open(os.path.join(wdir, "XPARM.XDS"), "w").write(xparm_str)

    made_files = []

    try:
        # 2. Modify XDS.INP
        modify_params = ("JOB", "INTEGRATE"), ("MINPK", "0"), ("DATA_RANGE", "%d %d"%(frame_num,frame_num))
        if sigmar is not None: modify_params += ("REFLECTING_RANGE_E.S.D.", "%.4f"%sigmar), 
        if sigmab is not None: modify_params += ("BEAM_DIVERGENCE_E.S.D.", "%.4f"%sigmab),
            
        modify_xdsinp(xdsinp, modify_params)

        # 3. Run xds
        p = subprocess.Popen("xds", cwd=wdir)
        p.wait()
    
        # 4&5. Rename FRAME_%d.cbf, INTEGRATE_%d.{HKL,LP} and make .adx 
        for f in ("FRAME.cbf", "INTEGRATE.LP", "INTEGRATE.HKL"):
            if not os.path.isfile(os.path.join(wdir, f)):
                continue
            sp = os.path.splitext(f)
            dst = os.path.join(wdir, "%s%s_%.4d%s" % (out_prefix, sp[0], frame_num, sp[1]))
            os.rename(os.path.join(wdir, f), dst)
            made_files.append(dst)

        if need_adx:
            made_files.append(make_adxfile(os.path.join(wdir, "%sINTEGRATE_%.4d.HKL" % (out_prefix, frame_num))))
            
    finally:
        # 6. Revert XDS.INP, etc.
        revert_files(backup_needed, bk_prefix, wdir=wdir)

    return made_files
# run()

if __name__ == "__main__":

    parser = optparse.OptionParser(usage="usage: %prog [options] setting-parameter-file")
    parser.add_option("--wdir", action="store", dest="wdir", default=".", type=str,
                      help="Working directory")
    parser.add_option("--frame","-f", action="store", dest="frame", type=int,
                      help="Frame number")

    (opts, args) = parser.parse_args(sys.argv)

    if not opts.frame:
        parser.print_help()
        sys.exit()

    if len(args) > 1 and os.path.isfile(args[1]):
        param_source = args[1]
    else:
        param_source = os.path.join(opts.wdir, "XPARM.XDS")

    run(param_source, opts.frame, opts.wdir)

    print() 
    print() 
    print("Launch:")
    print("adxv FRAME_%.4d.cbf INTEGRATE_%.4d.adx" % (opts.frame, opts.frame))
