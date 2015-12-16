#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
"""
How well can XDS estimate the full-intensity from partially recorded reflection on single frame?

This script does:
  1. Run XDS integration for a range of frames including a target frame.
  2. Run XDS integration for the target frame.
  3. Compare the intensities.

TODO:
  - Need to use the same REFLECTING_RANGE_E.S.D.=, REFLECTING_RANGE=, BEAM_DIVERGENCE_E.S.D.=, BEAM_DIVERGENCE= in two XDS runs.

To launch this script,
PHENIX_TRUST_OTHER_ENV=1 phenix.python /.../xds_single_frame_integration_compare.py  -f 5 -o 2
"""

import subprocess
from yamtbx.xds import *
from yamtbx.xds.xparm import get_xparm_from_integrate_lp, XPARM
from yamtbx.xds import integrate_hkl_as_flex

def compare_integrate_hkls(result_single, result_range):
    from cctbx import miller
    from cctbx.array_family import flex

    def calc_r_fac(data1, data2):
        return flex.sum(flex.abs(data1 - data2))/flex.sum(data1)
    # calc_r_fac()

    def calc_cc(data1, data2):
        corr = flex.linear_correlation(data1, data2)
        assert corr.is_well_defined()
        return corr.coefficient()
    # calc_cc()

    read_columns = ["IOBS","SIGMA","PEAK"]
    reader_single = integrate_hkl_as_flex.reader(result_single, read_columns)
    reader_range = integrate_hkl_as_flex.reader(result_range, read_columns)

    assert len(reader_single.hkl) == len(set(reader_single.hkl)) # No duplicate hkls!
    assert len(reader_range.hkl) == len(set(reader_range.hkl)) # No duplicate hkls!

    data_single = reader_single.arrays()
    data_range = reader_range.arrays()

    # Take common sets
    pairs = miller.match_indices(reader_single.hkl, reader_range.hkl).pairs()

    # XXX. Make sure two arrays are sorted (having the same order)!
    for k in read_columns:
        data_single[k] = data_single[k].select(pairs.column(0))
        data_range[k] = data_range[k].select(pairs.column(1))

    ## Sort by PEAK
    #perm = flex.sort_permutation(data=data_single["PEAK"].data(), reverse=False)
    
    # Output to terminal
    N = 20
    print "      PEAK                          r_fac  cc"
    for i in xrange(N):
        lrange, rrange = i * 100./N, (i+1)*100./N
        sel = lrange < data_single["PEAK"].data()
        sel &= data_single["PEAK"].data() <= rrange
        sel &= 0.95 < data_range["PEAK"].data() 

        sel_I_range = data_range["IOBS"].select(sel)
        sel_I_single = data_single["IOBS"].select(sel)

        r_fac = calc_r_fac(sel_I_range.data(), 
                           sel_I_single.data())

        cc = calc_cc(sel_I_range.data(), 
                     sel_I_single.data())

        for r, s, p in zip(sel_I_range[:10], sel_I_single[:10], data_single["PEAK"].select(sel).data()):
            print "  ", r, s, p

        print "%6.2f .. %6.2f [nref= %10d] %.4f %.4f" % (lrange, rrange, sum(sel), r_fac, cc)

    # Output to file
    ofs = open("result.dat", "w")
    ofs.write("H K L single.I single.SIGMA single.PEAK range.I range.SIGMA range.PEAK\n")
    for values in zip(data_single["IOBS"].indices(),
                                           data_single["IOBS"].data(), data_single["SIGMA"].data(), data_single["PEAK"].data(),
                                           data_range["IOBS"].data(), data_range["SIGMA"].data(), data_range["PEAK"].data()):
        ofs.write("%d %d %d %f %f %f %f %f %f\n" % ((values[0][0], values[0][1], values[0][2])+values[1:]))

    
# compare_integrate_hkls()

def run_xds(frame, offset, fix_setting_params=True):
    xdsinp = "XDS.INP"
    backup_needed = ("XPARM.XDS", "XDS.INP", "INTEGRATE.LP", "INTEGRATE.HKL", "FRAME.cbf")

    # Back up files
    bk_prefix = make_backup(backup_needed)
    
    ####
    # Run first integration
    modify_xdsinp(xdsinp, [("JOB", "INTEGRATE"), ("MINPK", "1"), 
                           ("DATA_RANGE", "%d %d"%(frame-offset,frame+offset)),
                           ("DELPHI", "%d"%(2*offset+1))]) # XXX. Assuming no osc_range > 1..

    p = subprocess.Popen("xds")
    p.wait()
    os.rename("INTEGRATE.LP", "INTEGRATE_around%d.LP" % frame)
    os.rename("INTEGRATE.HKL", "INTEGRATE_around%d.HKL" % frame)

    revert_files(backup_needed, bk_prefix)
    
    # Back up files for second run
    bk_prefix = make_backup(backup_needed)

    ####
    # Run second integration (single frame).
    if fix_setting_params: # Take parameters from the first run!!
        modify_xdsinp(xdsinp, [("JOB", "INTEGRATE"), ("MINPK", "1"), 
                               ("DATA_RANGE", "%d %d"%(frame,frame)),
                               ("DELPHI", "%d"%(2*offset+1)), # XXX. Assuming no osc_range > 1..
                               ("REFINE(INTEGRATE)", "")])
    
        xparm_str = get_xparm_from_integrate_lp("INTEGRATE_around%d.LP" % frame, frame)
        open("XPARM.XDS", "w").write(xparm_str)
    else:
        modify_xdsinp(xdsinp, [("JOB", "INTEGRATE"), ("MINPK", "1"), 
                               ("DATA_RANGE", "%d %d"%(frame,frame)),
                               ("DELPHI", "%d"%(2*offset+1))]) # XXX. Assuming no osc_range > 1..
    
    p = subprocess.Popen("xds")
    p.wait()

    os.rename("INTEGRATE.LP", "INTEGRATE_only%d.LP" % frame)
    os.rename("INTEGRATE.HKL", "INTEGRATE_only%d.HKL" % frame)

    revert_files(backup_needed, bk_prefix)

def run(frame, offset):

    # Run XDS twice
    run_xds(frame, offset)

    # ..then compare INTEGRATE.HKL
    # 1. take common sets of them, (make sure no duplicate hkls!!)
    # 2. sort by PEAK and calculate R-factor or CC (or chi**2?)
    compare_integrate_hkls("INTEGRATE_only%d.HKL" % frame,
                           "INTEGRATE_around%d.HKL" % frame)
    
if __name__ == "__main__":
    import optparse,sys

    parser = optparse.OptionParser(usage="usage: %prog [options] setting-parameter-file")
    parser.add_option("--frame","-f", action="store", dest="frame", type=int,
                      help="Target frame number")
    parser.add_option("--offset","-o", action="store", dest="offset", type=int, default=2,
                      help="Offset (number of frames before/after the target frame")

    (opts, args) = parser.parse_args(sys.argv)

    if not opts.frame:
        parser.print_help()
        sys.exit()

    run(frame=opts.frame, offset=opts.offset)
