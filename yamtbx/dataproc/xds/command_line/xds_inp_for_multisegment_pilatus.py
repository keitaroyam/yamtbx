"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.command_line.pilatus_for_crystfel import make_geom_decomposed_for_raw
from yamtbx.dataproc.XIO import XIO
from yamtbx.dataproc import crystfel
from yamtbx.dataproc.crystfel.command_line.geom_for_xds import geom_to_xdsinp_str

def run(img_in):
    #im = XIO.Image(img_in)
    #make_geom_decomposed_for_raw(im.header, (2527, 2463), "junk.geom")

    #geom = crystfel.geom.Geomfile("junk.geom")
    #print geom_to_xdsinp_str(geom)
    shape = (2527, 2463)
    panel = 0

    if shape == (2527, 2463):
        for j in xrange(12):
            for i in xrange(5):
                panel += 1
                print """\
 SEGMENT= %(min_fs)d %(max_fs)d %(min_ss)d %(max_ss)d ! %(panel).2d
  DIRECTION_OF_SEGMENT_X-AXIS= 1.000000 0.000000 0
  DIRECTION_OF_SEGMENT_Y-AXIS= 0.000000 1.000000 0
  SEGMENT_ORGX= 0
  SEGMENT_ORGY= 0
  SEGMENT_DISTANCE= 0
  REFINE_SEGMENT=\
""" % dict(panel=panel, 
           min_ss=j*212+1, max_ss=j*212+195,
           min_fs=i*494+1, max_fs=i*494+487
           )



if __name__ == "__main__":
    import sys
    run(sys.argv[1])
