#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc import cbf
import os

def run(cbfin, cbfout, repl=None):
    data, ndimfast, ndimmid = cbf.load_minicbf_as_numpy(cbfin)
    max_I = max(data)

    if repl is None:
        repl = max_I+100

    data[data<0] = repl
    cbf.save_numpy_data_as_cbf(data, ndimfast, ndimmid, "negative_replaced", cbfout)

    print "Launch:"
    print "adxv -overload %d %s" % (max_I, cbfout)
# run()

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print "Usage: %s ?????.cbf [replace-number]" % sys.argv[0]
        quit()

    repl = None
    if len(sys.argv) > 2:
        repl = int(sys.argv[2])

    run(sys.argv[1], 
        os.path.splitext(os.path.basename(sys.argv[1]))[0]+"_nr.cbf",
        repl)

