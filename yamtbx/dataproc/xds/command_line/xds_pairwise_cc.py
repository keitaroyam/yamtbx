"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
import sys

from cctbx.array_family import flex
from cctbx import miller

from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII

def run(lstin):
    data = []
    for l in open(lstin):
        xdsasc = l.strip()
        xa = XDS_ASCII(xdsasc, sys.stdout, i_only=True)
        ma = miller.array(miller_set=xa.as_miller_set(anomalous_flag=False),
                          data=xa.iobs)
        data.append((xdsasc, ma))

    print("index filename")
    for i, d in enumerate(data):
        print(i, d[0])

    print("i j n.i n.j n.common cc")
    for i in range(len(data)-1):
        for j in range(i+1, len(data)):
            di, dj = data[i][1].common_sets(data[j][1], assert_is_similar_symmetry=False)
            print(i, j, data[i][1].data().size(), data[j][1].data().size(), end=' ') 
            if len(di.data()) == 0:
                print(0, "nan")
            else:
                corr = flex.linear_correlation(di.data(), dj.data())
                assert corr.is_well_defined()
                cc =  corr.coefficient()
                print(len(di.data()), cc)
# run()

if __name__ == "__main__":
    lst = sys.argv[1]
    run(lst)
