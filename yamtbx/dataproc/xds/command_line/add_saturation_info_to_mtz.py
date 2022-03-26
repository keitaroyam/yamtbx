#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
"""
Add staturation flags (SATURATED_PART, SATURATED_ALL) to mtz file.
Read MAXC column from INTEGRATE.HKL and compare it with user-specified overload= value.
For a unique reflection, if more than one (or, all) observed reflection is saturated, the part flag (or, all) is set True.
First of all, read XDS_ASCII.HKL and reject reflections with negative sigma!!

Usage: cctbx.python add_saturation_info_to_mtz.py XDS_ASCII_copy_free.mtz ../INTEGRATE.HKL ../XDS_ASCII.HKL overload=65535 hklout=test2.mtz
"""
from __future__ import print_function
from __future__ import unicode_literals

master_params_str = """\
hklin = None
    .type = path
    .help = MTZ file (input)
hklout = None
    .type = path
    .help = MTZ file (output)
integrate_hkl = None
    .type = path
    .help = INTEGRATE.HKL file
xds_ascii= None
    .type = path
    .help = XDS_ASCII.HKL file
overload = None
    .type = int
    .help = Overload value (65535 for 16 bit)
"""

import iotbx.phil
import iotbx.mtz
from cctbx import miller
from cctbx.array_family import flex
from yamtbx.dataproc.xds import integrate_hkl_as_flex
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII

def run(params, log_out):
    xa = XDS_ASCII(params.xds_ascii, log_out)
    rejected_array = miller.array(miller_set=miller.set(crystal_symmetry=xa.symm,
                                                        indices=xa.indices,
                                                        anomalous_flag=False),
                                  data=xa.sigma_iobs < 0)
    xa_zd = miller.array(miller_set=miller.set(crystal_symmetry=xa.symm,
                                               indices=xa.indices,
                                               anomalous_flag=False),
                         data=xa.zd)

    # Read ZCAL, not ZOBS, because ZOBS (and XOBS, YOBS) can be zero (in case unobserved).
    integ_data = integrate_hkl_as_flex.reader(params.integrate_hkl, ["MAXC","ZCAL"]).arrays()
    maxc_array, integ_zcal = integ_data["MAXC"], integ_data["ZCAL"]

    assert integ_zcal.unit_cell().is_similar_to(xa_zd.unit_cell()) # two set of indices should be comparable.

    overload_flags = maxc_array.customized_copy(data=maxc_array.data() == params.overload)
    print("Overloaded observations in INTEGRATE.HKL:", overload_flags.data().count(True))
    print("Rejected (sigma<0) observations in XDS_ASCII.HKL:", rejected_array.data().count(True))
    # common_sets() does not work correctly for unmerged data!

    rejected_zd = xa_zd.select(rejected_array.data())

    #reject_indices = flex.bool([False for i in xrange(overload_flags.size())])

    print("making indices...........")
    import yamtbx_utils_ext
    integ_zcal = integ_zcal.sort(by_value="packed_indices") # Must be sorted before C++ function below!!
    reject_indices = yamtbx_utils_ext.make_selection_for_xds_unmerged(rejected_zd.indices(),
                                                                      rejected_zd.data(),
                                                                      integ_zcal.indices(),
                                                                      integ_zcal.data(),
                                                                      3.)
    """
    # This loop is too slow!
    for i in xrange(rejected_zd.size()):
        sel = integ_zcal.indices() == rejected_zd.indices()[i]
        sel &= (integ_zcal.data() - rejected_zd.data()[i]) < 3
        reject_indices.set_selected(sel, True)
        print i, rejected_zd.size(), sel.count(True)
        """
    """
    # This loop is also too slow!
    for j in xrange(integ_zcal.size()): # j: INTEGRATE.HKL
        if rejected_zd.indices()[i] != integ_zcal.indices()[j]:
            continue
        if abs(rejected_zd.data()[i] - integ_zcal.data()[j]) < 3: # within 3 frames.. OK?
            reject_indices[j] = True
    """

    print("Found rejected observations in INTEGRATE.HKL:", reject_indices.count(True))
    overload_flags.data().set_selected(reject_indices, False) # Set 'Un-overloaded'
    print("Remained overloaded observations:", overload_flags.data().count(True))

    overload_flags_partial = overload_flags.map_to_asu().merge_equivalents(incompatible_flags_replacement=True).array()
    overload_flags_all = overload_flags.map_to_asu().merge_equivalents(incompatible_flags_replacement=False).array()

    mtz_object = iotbx.mtz.object(params.hklin). \
        add_crystal("crystal", "project", overload_flags_all.unit_cell()). \
        add_dataset(name="dataset", wavelength=0). \
        add_miller_array(miller_array=overload_flags_all, column_root_label="SATURATED_ALL"). \
        add_miller_array(miller_array=overload_flags_partial, column_root_label="SATURATED_PART"). \
        mtz_object()
    mtz_object.write(file_name=params.hklout)


if __name__ == "__main__":
    import sys
    import os

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()

    for arg in cmdline.remaining_args:
        if params.hklin is None and arg.endswith(".mtz"):
            params.hklin = arg
        elif params.integrate_hkl is None and "INTEGRATE" in arg:
            params.integrate_hkl = arg
        elif params.xds_ascii is None and "XDS_ASCII" in arg:
            params.xds_ascii = arg

    print("Paramters:")
    cmdline.work.format(python_object=params).show(out=sys.stdout, prefix=" ")
    print()

    if None in (params.hklin, params.xds_ascii, params.integrate_hkl, params.overload):
        print("Missing information!")
        sys.exit(1)

    if params.hklout is None:
        params.hklout = os.path.splitext(os.path.basename(params.hklin))[0] + "_olflag.mtz"

    run(params, sys.stdout)
