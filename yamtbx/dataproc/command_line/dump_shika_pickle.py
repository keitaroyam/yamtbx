"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import pickle


def run(pklin):
    stats = pickle.load(open(pklin, "rb"))

    print "Files in this pickle:"
    for f in stats:
        print "  %s" % f

    print
    print

    for f in stats:
        print f
        print "=" * len(f)
        stat = stats[f]
        print " Detector=", stat.detector
        print " Gonio=", stat.gonio
        print " Grid coord=", stat.grid_coord
        print " Params=", "rawname:",stat.params.distl.image  #dir(stat.params)
        print " Scan info=", "wavelen:", stat.scan_info.wavelength
        print " Spots=", stat.spots.get_n_spots("all")
        print " Thumb pos mag=", stat.thumb_posmag
        print

if __name__ == "__main__":
    import sys
    pklin = sys.argv[1]
    run(pklin)
