#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc import cbf

def run(files, cbfout):
    merged = None
    for i, f in enumerate(files):
        repl = -10 * (i+1)
        print "%s %d" % (f, repl)

        data, ndimfast, ndimmid = cbf.load_minicbf_as_numpy(f)
        if i == 0:
            merged = data.copy()
            continue

        merged[data==-10] = 65540 # repl # for adxv visualization. only assuming two files.

    cbf.save_numpy_data_as_cbf(merged, ndimfast, ndimmid, "merged_predictions", cbfout)



if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print "Usage: %s FRAME.1.cbf FRAME.2.cbf .." % sys.argv[0]
        quit()

    run(sys.argv[1:], "FRAME_merged.cbf")
