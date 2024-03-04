#!/usr/bin/env phenix.python

from __future__ import print_function
from __future__ import unicode_literals
import sys, os, optparse
import iotbx.mtz
from cctbx import crystal
from cctbx import r_free_utils
from iotbx.reflection_file_editor import is_rfree_array, get_best_resolution, get_original_array_types
from iotbx.reflection_file_utils import get_r_free_flags_scores, make_joined_set

def create_mtz_dataset(miller_arrays):
    def get_label(m):
        l = m.info().labels[0]
        if l.endswith(("(+)", "(-)")):
            l = l[0:-3]
        return l
    # get_label()
    mtz_dataset = miller_arrays[0].as_mtz_dataset(column_root_label=get_label(miller_arrays[0]))
    for arr in miller_arrays[1:]:
        mtz_dataset.add_miller_array(miller_array=arr, column_root_label=get_label(arr))
    return mtz_dataset
# create_mtz_dataset()

def run(mtz, mtz_out, mtz_ref, flag_name=None, flag_value=None):

    ref_arrays = iotbx.mtz.object(mtz_ref).as_miller_arrays()
    print("Opening reference:", mtz_ref)

    # Get flag array
    flag_array = None
    if flag_name is None:
        flags = [x for x in ref_arrays if is_rfree_array(x, x.info())]
        if len(flags) == 0:
            print(" No R free flags like column found.")
            quit()
        elif len(flags) > 1:
            print(" More than one column which looks like R free flag:")
            for f in flags:
                print(" ", f.info().label_string())
            quit()
        else:
            flag_name = flags[0].info().label_string()
            flag_array = flags[0]
            print(" Guessing R free flag:", flag_name)
    else:
        flags = [x for x in ref_arrays if flag_name==x.info().label_string()]
        if len(flags) == 0:
            print(" Specified flag name not found:", flag_name)
            quit()
        else:
            print(" Use specified flag:", flag_name)
            flag_array = flags[0]
            
    # Get flag number
    if flag_value is None:
        flag_scores = get_r_free_flags_scores(miller_arrays=[flag_array], test_flag_value=flag_value)
        flag_value = flag_scores.test_flag_values[0]
        print(" Guessing flag number:", flag_value)
    else:
        print(" Specified flag number:", flag_value)

    print(" d_max, d_min=", get_best_resolution([flag_array], flag_array.crystal_symmetry()))
    print(" Symm:", flag_array.space_group().info(), flag_array.unit_cell())
    print()

    # Open mtz
    miller_arrays = iotbx.mtz.object(mtz).as_miller_arrays()
    print("Opening", mtz)

    if flag_name in [arr.info().label_string() for arr in miller_arrays]:
        print("Error: The column %s already exists in the mtz file: %s" % (flag_name, mtz))
        quit()

    print(" Using information from", miller_arrays[0].info().label_string())
    input_symm = crystal.symmetry(
        unit_cell=miller_arrays[0].unit_cell(),
        space_group_info=miller_arrays[0].space_group().info(),
        assert_is_compatible_unit_cell=False,
        force_compatible_unit_cell=False)

    d_max, d_min = get_best_resolution(miller_arrays, input_symm)
    
    print(" d_max, d_min=", d_max, d_min)
    print(" Symm:", input_symm.space_group_info(), input_symm.unit_cell())
    print()

    # Extend flag
    complete_set = make_joined_set(miller_arrays).complete_set()
    r_free_flags = r_free_utils.extend_flags(
        r_free_flags=flag_array.customized_copy(crystal_symmetry=input_symm),
        test_flag_value=flag_value,
        array_label=flag_name,
        complete_set=complete_set,
        accumulation_callback=None,
        preserve_input_values=True,
        d_max=d_max,
        d_min=d_min,
        log=sys.stdout).common_set(complete_set) #resolution_filter(d_min=d_min-0.01)

    print()

    print(r_free_flags.customized_copy(data=r_free_flags.data()==flag_value).show_r_free_flags_info())


    # Write mtz file
    #mtz_dataset = create_mtz_dataset(miller_arrays)
    #mtz_dataset.add_miller_array(miller_array=r_free_flags, column_root_label=flag_name)
    #mtz_object = mtz_dataset.mtz_object()
    mtz_object = iotbx.mtz.object(mtz).add_crystal("crystal", "project", r_free_flags.unit_cell()). \
        add_dataset(name="dataset", wavelength=0). \
        add_miller_array(miller_array=r_free_flags, column_root_label=flag_name).mtz_object()
    #mtz_object.add_history("copy and extend test flag from "+mtz_ref)
    #mtz_object.show_summary(out=sys.stdout, prefix="  ")
    mtz_object.write(file_name=str(mtz_out))

    print()
    print("Writing:", mtz_out)
    print() 
# run()

if __name__ == "__main__":
    parser = optparse.OptionParser(usage="usage: %prog [options] mtzfile",
                                   description="Import and extend Free R flag.")

    parser.add_option("--mtz-ref","-r", action="store", dest="mtz_ref", type=str, help="reference MTZ file")
    parser.add_option("--flag","-f", action="store", dest="flag", type=str, help="flag name: e.g. FreeR_flag. If unspecified, it will be automatically selected.")
    parser.add_option("--flag-value","-v", action="store", dest="flag_value", type=int, help="flag value. If unspecified, it will be automatically selected.")
    parser.add_option("--seed", action="store", dest="seed", type=int, help="")
    parser.add_option("--mtz-out","-o", action="store", dest="mtz_out", type=str, help="output MTZ file")
    (opts, args) = parser.parse_args(sys.argv)

    if len(args) < 2 or opts.mtz_ref is None:
        print(parser.print_help())
        quit()

    if not os.path.isfile(args[1]):
        print("File not found:", args[1])
        quit()

    if not os.path.isfile(opts.mtz_ref):
        print("File not found:", opts.mtz_ref)
        quit()

    
    if opts.seed is not None:
        random.seed(opts.seed)

    if opts.mtz_out is None:
        opts.mtz_out = os.path.splitext(os.path.basename(args[1]))[0] + "_copy_free.mtz"

    if os.path.isfile(opts.mtz_out):
        print("File already exists:", opts.mtz_out)
        print("Please remove the file or change the output name with -o option")
        quit()
    
    run(mtz=args[1], mtz_ref=opts.mtz_ref, mtz_out=opts.mtz_out, flag_name=opts.flag, flag_value=opts.flag_value)
