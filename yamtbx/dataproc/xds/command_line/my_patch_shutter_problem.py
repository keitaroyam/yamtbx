#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
from yamtbx.dataproc.xds import xds_ascii
import math

def run(xds_ascii_in):
    reader = xds_ascii.XDS_ASCII(xds_ascii_in, sys.stdout)
    #reader.remove_rejected()

    data = [[] for i in range(10)]

    for iobs, zd in zip(reader.iobs, reader.zd):
        dec = math.modf(zd)[0]
        if dec < 0:
            dec += 1

        idx = int(math.ceil((dec-0.05)/0.1))
        if idx > 9:
            idx = 0

        data[idx].append(iobs)

    iobs_overall = sum(map(lambda x:sum(x), data)) / sum(map(lambda x:len(x), data))

    print("   INTERVAL     NUMBER      INTENSITY      FACTOR")
    for i, iobs in enumerate(data):
        iobs_interval = old_div(sum(iobs),len(iobs))
        print(" % .2f % .2f %9d %14.3f %11.3f"%(-0.05+i*.1, -0.05+i*.1+.1, len(iobs), iobs_interval, iobs_overall/iobs_interval))

    print()
    print(" % .2f % .2f %9d %14.3f %11.3f"%(-0.05, 0.95, reader.iobs.size(), iobs_overall, 1.))

if __name__ == "__main__":
    ascin = sys.argv[1]
    run(ascin)
