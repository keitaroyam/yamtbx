#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
"""
Find data sets, make directories for XDS runs, and prepare XDS.INP.
"""

from yamtbx.dataproc import dataset
from yamtbx.util import call

import os

def split_runs_based_on_pf_logfiles(logfiles):
    def get_min_max_frames(logfile):
        from yamtbx.dataproc.dataset import re_pref_num_ext
        nums = []
        for l in open(logfile):
            if len(l.split()) == 4:
                r = re_pref_num_ext.search(os.path.basename(l.split()[-1]))
                print l[-1]
                if r:
                    nums.append(int(r.group(2)))
        return min(nums), max(nums)

    res = [ get_min_max_frames(logfile) for logfile in logfiles]
    res.sort()
    return res

def run(search_dir, pf_logfiles=None):
    for img_template, min_frame, max_frame in dataset.find_data_sets(search_dir,skip_symlinks=False):

        print "Dataset found:"
        print " NAME_TEMPLATE_OF_DATA_FRAMES= %s" % img_template
        print " DATA_RNAGE= %d %d" % (min_frame, max_frame)

        pdir, ppath = os.path.split(img_template)
        if not os.path.exists(os.path.join(pdir, "xds")):
            os.mkdir(os.path.join(pdir, "xds"))
        else:
            assert os.path.isdir(os.path.join(pdir, "xds"))

        prefix = ppath[:ppath.index("?")]
        if not prefix.endswith("_"):
            prefix += "_"

        if pf_logfiles is None:
            wdir = os.path.join(pdir, "xds", "xds_%s%d-%d" % (prefix, min_frame, max_frame))
            if os.path.isfile(os.path.join(wdir, "XDS.INP")):
                print " Already exist:", os.path.join(wdir, "XDS.INP")
                continue

            os.mkdir(wdir)

            cmd = 'generate_XDS.INP "%s"' % os.path.relpath(img_template, wdir)
            print " launching", cmd
            call(cmd=cmd, wdir=wdir)
        else:
            for minf, maxf in split_runs_based_on_pf_logfiles(pf_logfiles):
                if maxf < min_frame or max_frame < minf:
                    continue
                if minf < min_frame:
                    minf = min_frame
                if maxf > max_frame:
                    maxf = max_frame

                wdir = os.path.join(pdir, "xds", "xds_%s%d-%d" % (prefix, minf, maxf))

                if os.path.isfile(os.path.join(wdir, "XDS.INP")):
                    print "Already exist:", os.path.isfile(os.path.join(wdir, "XDS.INP"))
                    continue

                os.mkdir(wdir)

                cmd = 'generate_XDS.INP "%s"' % os.path.relpath(img_template, wdir)
                print " launching", cmd
                call(cmd=cmd, wdir=wdir)
                # XXX dirty hacks..
                inp = open(os.path.join(wdir, "XDS.INP")).readlines()
                ofs = open(os.path.join(wdir, "XDS.INP"), "w")
                for l in inp:
                    if l.startswith("DATA_RANGE="):
                        l = "DATA_RANGE= %d %d\n" % (minf, maxf)
                    if l.startswith("SPOT_RANGE="):
                        l = "SPOT_RANGE= %d %d\n" % (minf, (maxf-minf)/2+minf)
                    ofs.write(l)
                ofs.close()
        print

# run()

if __name__ == "__main__":
    import sys

    pf_logfiles = None

    if len(sys.argv) > 1:
        wdir = sys.argv[1]
        if len(sys.argv) > 2:
            pf_logfiles = sys.argv[2:]
    else:
        wdir = os.getcwd()

    run(wdir, pf_logfiles)
