from __future__ import print_function
from __future__ import unicode_literals
from yamtbx.dataproc import crystfel
import numpy

def geom_to_xdsinp_str(geom):
    ret = ""

    for sid in sorted(geom.sensors):
        s = geom.sensors[sid]
        ret += " SEGMENT= %d %d %d %d ! %s\n" % (s.min_fs+1, s.max_fs+1,
                                                 s.min_ss+1, s.max_ss+1,
                                                 sid)
        ret += "  DIRECTION_OF_SEGMENT_X-AXIS= %f %f 0\n" % tuple(s.fs)
        ret += "  DIRECTION_OF_SEGMENT_Y-AXIS= %f %f 0\n" % tuple(s.ss)
        #print "XXX=", numpy.cross(s.fs, s.ss)

        orgx, orgy = -numpy.dot(numpy.linalg.inv(numpy.vstack((s.fs, s.ss)).transpose()), (s.corner_x, s.corner_y))

        #print " SEGMENT_ORGX= %f" % s.corner_x
        #print " SEGMENT_ORGY= %f" % s.corner_y
        ret += "  SEGMENT_ORGX= %f\n" % (orgx+s.min_fs+1)
        ret += "  SEGMENT_ORGY= %f\n" % (orgy+s.min_ss+1)
        ret += "  SEGMENT_DISTANCE= %f\n" % (float(geom.clen)*1.e3)
        ret += "  REFINE_SEGMENT= \n"

    return ret
# geom_to_xdsinp_str()

def run(geom_in):
    geom = crystfel.geom.Geomfile(geom_in)
    print(geom_to_xdsinp_str(geom))
# run()

if __name__ == "__main__":
    import sys
    geom_in = sys.argv[1]

    run(geom_in)
