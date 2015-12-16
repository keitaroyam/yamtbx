#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import os
from yamtbx.dataproc.xds.xparm import get_xparm_from_integrate_lp

def run(integrate_lp, frame_num):
    xparm_str = get_xparm_from_integrate_lp(integrate_lp, frame_num)
    open("XPARM.XDS.%.4d"%frame_num, "w").write(xparm_str)
# run()

if __name__ == "__main__":
    import sys

    integrate_lp = sys.argv[1]
    frame = int(sys.argv[2])
    run(integrate_lp, frame)
