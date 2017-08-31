#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.util import call
from yamtbx.util.xtal import format_unit_cell
import iotbx.file_reader
import iotbx.scalepack.merge
import iotbx.shelx.hklf
from iotbx import crystal_symmetry_from_any
import os

def check_symm(xs_hkl, xs_pdb):
    pdb_laue = xs_pdb.space_group().build_derived_reflection_intensity_group(False)
    hkl_laue = xs_hkl.space_group().build_derived_reflection_intensity_group(False)
    if pdb_laue != hkl_laue or not xs_pdb.unit_cell().is_similar_to(xs_hkl.unit_cell()):
        print "WARNING! Incompatible symmetry (symmetry mismatch or too different unit cell)"
        print " hklin:"
        xs_hkl.show_summary(prefix="  ")
        print " pdbin:"
        xs_pdb.show_summary(prefix="  ")


def run(hklin, pdbin):
    arrays = iotbx.file_reader.any_file(hklin).file_server.miller_arrays
    i_arrays = filter(lambda x:x.is_xray_intensity_array() and x.anomalous_flag(),
                      arrays)
    f_arrays = filter(lambda x:x.is_xray_amplitude_array() and x.anomalous_flag(),
                      arrays)

    if not i_arrays and not f_arrays:
        print "No anomalous observation data"
        return

    wdir = "anode"
    if os.path.exists(wdir):
        print "%s already exists. quiting." % wdir
        return

    os.mkdir(wdir)

    xs = crystal_symmetry_from_any.extract_from(pdbin)

    sh_out = open(os.path.join(wdir, "run_anode.sh"), "w")
    sh_out.write("#!/bin/sh\n\n")
    sh_out.write("shelxc anode <<+ > shelxc.log 2>&1\n")
    sh_out.write("cell %s\n" % format_unit_cell(xs.unit_cell()))
    sh_out.write("spag %s\n" % str(xs.space_group_info()).replace(" ",""))

    if i_arrays:
        obs_array = i_arrays[0]
        infile = "%s.hkl" % os.path.splitext(os.path.basename(hklin))[0]
        in_opt = "%s" % infile
        print "Using intensity array:", obs_array.info().label_string()
    else:
        obs_array = f_arrays[0]
        infile = "%s_f.hkl" % os.path.splitext(os.path.basename(hklin))[0]
        in_opt ="-f %s" % infile
        print "No intensity arrays. Using amplitude arrays instead:", obs_array.info().label_string()
        
    sh_out.write("! data from %s : %s\n" % (os.path.abspath(hklin), obs_array.info().label_string()))
    obs_array.crystal_symmetry().show_summary(sh_out, prefix="! ")
    check_symm(obs_array.crystal_symmetry(), xs)
    sh_out.write("sad %s\n" % in_opt)
    iotbx.shelx.hklf.miller_array_export_as_shelx_hklf(obs_array, open(os.path.join(wdir, infile), "w"),
                                                       normalise_if_format_overflow=True)
    sh_out.write("+\n\n")
    sh_out.write('ln -s "%s" anode.pdb\n\n' % os.path.relpath(pdbin, wdir))
    sh_out.write("anode anode\n")
    sh_out.close()

    call(cmd="sh", arg="./run_anode.sh", wdir=wdir)

    print "\nDone. See %s/" % wdir
# run()

if __name__ == "__main__":
    import sys

    hklin = sys.argv[1]
    pdbin = sys.argv[2]
    run(hklin, pdbin)
