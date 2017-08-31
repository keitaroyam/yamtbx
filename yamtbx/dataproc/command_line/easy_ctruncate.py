#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import iotbx.phil
import iotbx.file_reader
from yamtbx.util import call
from yamtbx.util import get_temp_filename
import os

master_params_str = """ 
hklin = None
 .type = path
 .help = "Input file"
hklout = None
 .type = path
 .help = "Output file"
logout = None
 .type = path
 .help = "Logfile name"
d_min = None
 .type = float
d_max = None
 .type = float
"""

def run(params):
    i_arrays = filter(lambda x:x.is_xray_intensity_array(),
                      iotbx.file_reader.any_file(params.hklin).file_server.miller_arrays)
    f_arrays = filter(lambda x:x.is_xray_amplitude_array(),
                      iotbx.file_reader.any_file(params.hklin).file_server.miller_arrays)

    # TODO Copy non-intensity arrays to new mtz!

    if not i_arrays and not f_arrays:
        print "No observation data"
        return

    opts = ""

    if i_arrays:
        arrays = i_arrays
        print "Intensity arrays:"
    else:
        arrays = f_arrays
        opts = " -amplitudes"
        print "No intensity arrays. Using amplitude arrays instead:"

    for ar in arrays:
        print "", ar.info().label_string()
    print

    colin, colano = "", ""

    ano = filter(lambda x: x.anomalous_flag(), arrays)
    noano = filter(lambda x: not x.anomalous_flag(), arrays)

    if ano:
        colano = "/*/*/[%s]" % ano[0].info().label_string().replace(",merged","")
    if noano:
        colin = "/*/*/[%s]" % noano[0].info().label_string()

    hklin = params.hklin
        
    if (params.d_min, params.d_max).count(None) < 2:
        hklin = get_temp_filename(suffix=".mtz")
        d_min, d_max = 0, 10000
        if params.d_min is not None: d_min = params.d_min
        if params.d_max is not None: d_max = params.d_max
        call(cmd="mtzutils",
             arg="hklin %s hklout %s" % (params.hklin, hklin),
             stdin="resolution %f %f" % (d_min, d_max))

    cmd = "ctruncate -hklin %s -hklout %s " % (hklin, params.hklout)
    if colin != "":
        cmd += "-colin '%s' " % colin
    if colano != "":
        cmd += "-colano '%s' " % colano

    cmd += opts

    print cmd

    call(cmd=cmd,
         stdout=open(params.logout, "w")
         )

    if hklin != params.hklin:
        os.remove(hklin)
        log_content = open(params.logout).read()
        open(params.logout, "w").write(log_content.replace(hklin, params.hklin))
        print "NOTE: %s was modified to show the original input file name" % params.logout 
        
# run()

def run_from_args(argv):
    if not argv or "-h" in argv or "--help" in argv:
        iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
        return

    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    for arg in args:
        if not os.path.isfile(arg):
            print "File not found: %s" % arg
            return

        if params.hklin is None and arg.endswith(".mtz"):
            params.hklin = arg

    if params.hklin is None:
        print "Please give mtz file."
        return

    if not params.hklout: params.hklout = os.path.splitext(os.path.basename(params.hklin))[0] + "_ctruncate.mtz"
    if not params.logout: params.logout = "ctruncate_%s.log" % os.path.splitext(os.path.basename(params.hklin))[0]
        
    run(params)
# run_from_args()

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
