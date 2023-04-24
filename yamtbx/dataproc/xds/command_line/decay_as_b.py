"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from yamtbx.dataproc import cbf
import numpy
import re

"""
header example (Hlg):

 INPUT_FILE=INTEGRATE.HKL
 REC. CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 XMIN=     3.5 XMAX=   379.5 NXBIN=   38
 YMIN= 0.00042 YMAX= 0.15998 NYBIN=   20
"""

re_x = re.compile("XMIN= *([0-9\.]+) *XMAX= *([0-9\.]+) *NXBIN= *([0-9]+)")
re_y = re.compile("YMIN= *([0-9\.]+) *YMAX= *([0-9\.]+) *NYBIN= *([0-9]+)")

def run(cbfin):
    header, data, nslow, nfast = cbf.load_xds_special(cbfin)
    assert "REC. CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION" in header
    #print header
    data = numpy.array(data, dtype=float) / 1000
    data = data.reshape((nslow, nfast))

    xmin, xmax, nxbin = re_x.search(header).groups() # batch number
    ymin, ymax, nybin = re_y.search(header).groups() # 1/d^2
    xmin, xmax, ymin, ymax = list(map(float, (xmin, xmax, ymin, ymax)))
    nxbin, nybin = list(map(int, (nxbin, nybin)))

    xstep = (xmax-xmin)/nxbin
    ystep = (ymax-ymin)/nybin


    print("1/d^2 scale batch")
    for x in range(nfast):
        for j, d in enumerate(data[:,x]):
            print(ymin+(j+.5)*ystep, d, x)

    # TODO fit B here.

# run()

if __name__ == "__main__":
    import sys

    run(sys.argv[1])
