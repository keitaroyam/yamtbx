#!/usr/bin/env phenix.python

import sys
import os
import iotbx.mtz
from cctbx import crystal
from cctbx import r_free_utils
from iotbx.reflection_file_editor import get_best_resolution
from iotbx.reflection_file_utils import make_joined_set

def run(mtz, mtz_out, fraction, flag_name=None, ccp4=True, use_lattice_symmetry=True, n_shells=20, log_out=sys.stdout):

    # Open mtz
    miller_arrays = iotbx.mtz.object(mtz).as_miller_arrays()
    print >>log_out, "Opening", mtz
    print >>log_out, " Using information from"
    miller_arrays[0].show_summary(log_out, " ")

    input_symm = crystal.symmetry(
        unit_cell=miller_arrays[0].unit_cell(),
        space_group_info=miller_arrays[0].space_group().info(),
        assert_is_compatible_unit_cell=False,
        force_compatible_unit_cell=False)

    d_max, d_min = get_best_resolution(miller_arrays, input_symm)
    
    print >>log_out, " d_max, d_min=", d_max, d_min
    print >>log_out, " Symm:", input_symm.space_group_info(), input_symm.unit_cell()
    print >>log_out


    # Extend flag
    complete_set = make_joined_set(miller_arrays).complete_set()

    if use_lattice_symmetry:
        from cctbx.sgtbx import lattice_symmetry

        print >>log_out, "Lattice symmetry:"
        cb_op_to_niggli = complete_set.change_of_basis_op_to_niggli_cell()
        tmp_ma = complete_set.change_basis( cb_op_to_niggli )
        lattice_group = lattice_symmetry.group(tmp_ma.unit_cell(), max_delta=5.0)
        tmp_ma.show_summary(log_out, " ")
        print >>log_out

    new_r_free_array = complete_set.generate_r_free_flags(fraction=fraction,
                                                          max_free=None,
                                                          lattice_symmetry_max_delta=5.0,
                                                          use_lattice_symmetry=use_lattice_symmetry,
                                                          n_shells=n_shells)

    new_r_free_array.show_r_free_flags_info(out=log_out)

    if ccp4:
        new_r_free_array = new_r_free_array.customized_copy(data=r_free_utils.export_r_free_flags_for_ccp4(flags=new_r_free_array.data(), test_flag_value=True))


    print >>log_out

    # Write mtz file
    mtz_object = iotbx.mtz.object(mtz).add_crystal("crystal", "project", new_r_free_array.unit_cell()). \
        add_dataset(name="dataset", wavelength=0). \
        add_miller_array(miller_array=new_r_free_array, column_root_label=flag_name).mtz_object()
    mtz_object.write(file_name=mtz_out)

    print >>log_out
    print >>log_out, "Written:", mtz_out
    print >>log_out
# run()

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser(usage="usage: %prog [options] mtzfile",
                                   description="Create Free R flag.")

    parser.add_option("--seed", action="store", dest="seed", type=int, help="")
    parser.add_option("--ccp4", action="store_true", dest="ccp4_style", help="CCP4 Style")
    parser.add_option("-p", "--fraction", action="store", dest="fraction", default=0.05, type=float, help="Fraction")
    parser.add_option("-f", "--flag", action="store", dest="flag", default="FreeR_flag", help="Flag name")
    parser.add_option("--mtz-out","-o", action="store", dest="mtz_out", type=str, help="output MTZ file")

    (opts, args) = parser.parse_args(sys.argv)

    if len(args) < 2:
        print parser.print_help()
        quit()

    
    if opts.seed is not None:
        random.seed(opts.seed)

    if opts.mtz_out is None:
        opts.mtz_out = os.path.splitext(os.path.basename(args[1]))[0] + "_free.mtz"

    run(mtz=args[1], mtz_out=opts.mtz_out, flag_name=opts.flag, ccp4=opts.ccp4_style, fraction=opts.fraction)
