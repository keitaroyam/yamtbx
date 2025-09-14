#!/usr/bin/env phenix.python
"""
This is simplified version of iotbx.reflection_file_editor.
"""
from __future__ import print_function
from __future__ import unicode_literals
import sys, os, optparse
import string
import re
import iotbx.mtz
from cctbx.array_family import flex

from iotbx.reflection_file_editor import guess_array_output_labels


# Original is in iotbx.reflection_file_editor
def get_original_array_types (mtz_file, original_labels) :
    array_types = ""
    mtz_columns = mtz_file.column_labels()
    mtz_types = mtz_file.column_types()
    mtz_crossref = dict(list(zip(mtz_columns, mtz_types)))
    for label in original_labels :
        array_types += mtz_crossref[label]
    return array_types

# Original is in iotbx.reflection_file_editor
def add_array_to_mtz_dataset (mtz_dataset, output_array, fake_label, column_types) :
    if mtz_dataset is None :
        mtz_dataset = output_array.as_mtz_dataset(column_root_label=fake_label,
                                                  column_types=column_types)
    else :
        mtz_dataset.add_miller_array(miller_array=output_array,
                                     column_root_label=fake_label,
                                     column_types=column_types)

    return mtz_dataset

def run(mtz, bs):
    # Open mtz
    mtz_file = iotbx.mtz.object(mtz)
    miller_arrays = mtz_file.as_miller_arrays()
    print("Opening", mtz)

    for b in bs:
        mtz_out = os.path.splitext(os.path.basename(args[1]))[0] + "_b%.2f.mtz" % b
        mtz_dataset = None
        labels = ["H", "K", "L"]
        
        for i, ar in enumerate(miller_arrays):
            fake_label = 2 * string.ascii_uppercase[i]
            for lab in guess_array_output_labels(ar):
                labels.append(lab)
            array_types = get_original_array_types(mtz_file, ar.info().labels)
            default_types = iotbx.mtz.default_column_types(ar)
            if len(default_types) == len(array_types):
                column_types = array_types
            else:
                column_types = None

            if ar.is_xray_data_array() or ar.is_complex_array():
                fac = 2. if ar.is_xray_intensity_array() else 1.
                print("Applying B=%.2f to %s" % (b*fac, ar.info()))
                k = flex.exp(-ar.d_star_sq().data() * b*fac)
                ar = ar.customized_copy(data=ar.data()*k, sigmas=ar.sigmas()*k if ar.sigmas() else None)

            mtz_dataset = add_array_to_mtz_dataset(mtz_dataset, ar, fake_label,
                                                   column_types)

        # Decide labels and write mtz file
        mtz_object = mtz_dataset.mtz_object()
        invalid_chars = re.compile(r"[^A-Za-z0-9_\-+\(\)]")

        used = dict([ (label, 0) for label in labels ])

        for i, column in enumerate(mtz_object.columns()):
            if column.label() != labels[i] :
                label = labels[i]
                original_label = label

                assert used[label] == 0

                try:
                    column.set_label(label)
                except RuntimeError as e:
                    if ("new_label is used already" in str(e)) :
                        col_names = [ col.label() for col in mtz_object.columns() ]
                        raise RuntimeError(("Duplicate column label '%s': current labels "+
                                            "are %s; user-specified output labels are %s.") %
                                           (label, " ".join(col_names), " ".join(labels)))
                else:
                    used[original_label] += 1

        mtz_object.write(file_name=mtz_out)

        print()
        print("Writing:", mtz_out)
        print()
# run()

if __name__ == "__main__":
    parser = optparse.OptionParser(usage="usage: %prog [options] mtzfile")

    parser.add_option("-b", action="append", dest="b", type=float, help="like 2.5,3")
    (opts, args) = parser.parse_args(sys.argv)

    if len(args) < 2 or opts.b is None:
        print(parser.print_help())
        quit()

    if not os.path.isfile(args[1]):
        print("File not found:", args[1])
        quit()

    run(mtz=args[1], bs=opts.b)
