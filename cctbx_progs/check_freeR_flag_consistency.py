#!/usr/bin/env phenix.python
from __future__ import print_function
from __future__ import unicode_literals
import sys, optparse
import iotbx.mtz
from iotbx.reflection_file_editor import is_rfree_array
from iotbx.reflection_file_utils import get_r_free_flags_scores

def guess_flag_name_and_value(mtzin, name=None, val=None):
    arrays = iotbx.mtz.object(mtzin).as_miller_arrays()
    if name is None:
        flags = [x for x in arrays if is_rfree_array(x, x.info())]
        name = flags[0].info().label_string()
        flag_array = flags[0]
    else:
        flags = [x for x in arrays if name==x.info().label_string()]
        flag_array = flags[0]
    
    if val is None:
        flag_scores = get_r_free_flags_scores(miller_arrays=[flag_array], test_flag_value=val)
        val = flag_scores.test_flag_values[0]

    return name, val
# guess_flag_name_and_value()

def compare(mtzin_list, label_list, flag_list):
    ##
    # Check test set index consistency.
    # If there're test set indices which aren't consistent among given mtz files, 
    # this function displays the indices.
    #
    # @param mtzin_list list of mtz file names
    # @param label_list list of labels for test flag
    # @param flag_list  list of flag numbers which indicate test set
    #

    assert( len(mtzin_list) == len(label_list) == len(flag_list) )

    for m, l, f in zip(mtzin_list, label_list, flag_list):
        print(m, l, f)

    # number of mtz files
    num = len(mtzin_list)
    
    # list of mtz objects
    mtz_objs = [ iotbx.mtz.object(file_name=mtzin) for mtzin in mtzin_list ]
    
    # list of mtz columns (free-R)
    columns  = [ mtz_obj.get_column(label=label) for mtz_obj, label in zip(mtz_objs, label_list) ]

    # list of hkl indices
    hkls     = [ mtz_obj.extract_miller_indices() for mtz_obj in mtz_objs ]

    
    for_comp = {} # {(h,k,l): [flag_values for each file]}

    for i, (hkl, column) in enumerate(zip(hkls, columns)):
        for (h,k,l), flag in zip(hkl, column.extract_values()):
            for_comp.setdefault((h,k,l), [None]*num)[i] = flag

    n_match, n_notmatch = 0, 0

    # Compare flag numbers

    print("Inconsistent indices:")
    print()

    for (h,k,l), flags in sorted(for_comp.items()):
        bool_list = [ flag == ref for flag,ref in zip(flags, flag_list) ]

        # If more than 0 indices of test set included and they aren't consistent among files
        if 0 < sum(bool_list):
            if sum(bool_list) < len(bool_list):
                n_notmatch += 1
                print("%3d %3d %3d" % (h, k, l), " ".join([str(x) for x in flags]))
            else:
                n_match += 1

    print()
    print("SUMMARY")
    digit = str(len(str(max(n_match, n_notmatch))))
    print(("%"+digit+"d match") % n_match)
    print(("%"+digit+"d not match") % n_notmatch)


if __name__ == "__main__":
    
    parser = optparse.OptionParser(usage="usage: %prog -m hoge.mtz -l FreeR_flag -f 0 -m fuga.mtz -l R-free-flags -f 1",
                                   description="Check test set index consistency")

    parser.add_option("--mtz","-m", action="append", dest="mtz_list", type=str, 
                      help="MTZ files.")
    parser.add_option("--label","-l", action="append", dest="label_list", type=str, 
                      help="MTZ column labels. e.g. R-free-flags, FreeR_flag.")
    parser.add_option("--flag","-f", action="append", dest="flag_list", type=int,
                      help="Flag numbers which indicate.")

    (opts, args) = parser.parse_args(sys.argv)

    if not opts.mtz_list:
        parser.print_help()
        quit()

    if None in (opts.label_list, opts.flag_list):
        opts.label_list, opts.flag_list = [], []
        for i, mtz in enumerate(opts.mtz_list):
            name, val = guess_flag_name_and_value(mtz)
            opts.label_list.append(name)
            opts.flag_list.append(val)
    

    if not( len(opts.mtz_list) == len(opts.label_list) == len(opts.flag_list) ):
        print("Number doesn't match")
        exit(1)



    compare(mtzin_list= opts.mtz_list,
            label_list= opts.label_list,
            flag_list= opts.flag_list
            )


