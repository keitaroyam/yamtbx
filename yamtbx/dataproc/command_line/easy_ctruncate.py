#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.util import call
import iotbx.file_reader
import os
def run(hklin):
    arrays = filter(lambda x:x.is_xray_intensity_array(),
                    iotbx.file_reader.any_file(hklin).file_server.miller_arrays)

    # TODO Copy non-intensity arrays to new mtz!

    if len(arrays) == 0:
        print "No intensity array"
        return

    print "Intensity arrays:"
    for ar in arrays:
        print "", ar.info().label_string()
    print

    colin, colano = "", ""

    if len(arrays) == 1:
        coltmp = arrays[0].info().label_string()
        if arrays[0].anomalous_flag():
            colano = "/*/*/[%s]" % coltmp.replace(",merged","")
        else:
            colin = "/*/*/[%s]" % coltmp
    elif len(arrays) == 2 and sum(map(lambda x: x.anomalous_flag(), arrays)) == 1:
        tmpano = filter(lambda x: x.anomalous_flag(), arrays)[0]
        tmpnoano = filter(lambda x: not x.anomalous_flag(), arrays)[0]
        colano = "/*/*/[%s]" % tmpano.info().label_string().replace(",merged","")
        colin = "/*/*/[%s]" % tmpnoano.info().label_string()
    else:
        print "What should I do.. exiting."
        return

    hklout = os.path.splitext(os.path.basename(hklin))[0] + "_ctruncate.mtz"
    cmd = "ctruncate -hklin %s -hklout %s " % (hklin, hklout)
    if colin != "":
        cmd += "-colin '%s' " % colin
    if colano != "":
        cmd += "-colano '%s' " % colano

    call(cmd=cmd,
         stdout=open("ctruncate_%s.log" % os.path.splitext(os.path.basename(hklin))[0], "w")
         )

if __name__ == "__main__":
    import sys

    hklin = sys.argv[1]
    run(hklin)
