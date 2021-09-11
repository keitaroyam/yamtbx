"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
"""
Use this script after xscale_cc_against_merged.py
"""
from __future__ import print_function
from __future__ import unicode_literals

from collections import OrderedDict
import os
import shutil
from yamtbx.dataproc.xds import xds_ascii

def read_input(tabin, cc_min):
    ifs = open(tabin)
    header = ifs.readline().split()
    assert header[-1] == "cc"
    i_frame = header.index("frame")
    i_file = header.index("file")

    ret = OrderedDict()

    for l in ifs:
        sp = l.split()
        cc = float(sp[-1])
        if cc != cc or cc < cc_min: # NaN or lower than limit
            continue

        f = sp[i_file]

        if i_frame >= 0:
            n = int(sp[i_frame])
            ret.setdefault(f,[]).append(n)
        else:
            ret[f] = None # use all

    return ret
# read_input()

def run(tabin, cc_min):
    lst = read_input(tabin, cc_min)
    inp_ofs = open("XSCALE.INP", "w")

    topdir = os.path.dirname(os.path.commonprefix([os.path.abspath(x) for x in list(lst.keys())]))

    print("Cut frames with CC < %.4f" % cc_min)

    for f in lst:
        useframes = lst[f]
        xac = xds_ascii.XDS_ASCII(f)

        hklout = os.path.join("files", os.path.relpath(f, topdir))
        if not os.path.exists(os.path.dirname(hklout)): os.makedirs(os.path.dirname(hklout))

        if useframes is None or set(useframes).issuperset(set(range(min(xac.iframe), max(xac.iframe)))):
            print("Just copying to %s" % hklout)
            shutil.copyfile(f, hklout)
        else:
            sel = xac.iframe == useframes[0]
            for x in useframes[1:]: sel |= xac.iframe == x

            # XXX we may need to check I/sigma and throw < -3..
            if sum(sel) < 10:
                print("Skipping %s with only %6d/%6d" % (hklout, sum(sel), len(sel)))
                continue

            print("Saving %s %6d/%6d" % (hklout, sum(sel), len(sel)))
            xac.write_selected(sel, hklout)

        inp_ofs.write("INPUT_FILE= %s\n" % hklout)
# run()

if __name__ == "__main__":
    import sys

    tabin = sys.argv[1]
    cc_min = float(sys.argv[2])
    run(tabin, cc_min)
    
