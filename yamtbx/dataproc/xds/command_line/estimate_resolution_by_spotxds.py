"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

# The same method in CrystFEL
# libcrystfel/src/integration.c
#  estimate_resolution()

from yamtbx.dataproc.xds.idxreflp import SpotXds

def run(spot_xds, xparm_in):
    sx = SpotXds(spot_xds)
    sx.set_xparm(xparm_in)

    idxed = sx.indexed_and_unindexed_on_detector(with_resolution=True)["indexed"]
    acc = map(lambda x:1./x[-1], idxed) # list of 1/d

    if len(acc) < 3:
        print "WARNING: Too few peaks to estimate resolution."
        return 0

    # outlier removal
    acc.sort()
    n = len(acc)//50
    if n < 2: n = 2
    max_res = acc[(len(acc)-1)-n]

    return 1./max_res
# run()

if __name__ == "__main__":
    import sys
    spot_xds = sys.argv[1]
    xparm_in = sys.argv[2]
    print run(spot_xds, xparm_in)
