#!/usr/bin/env phenix.python

import os
import string
import re
import math
import iotbx.mtz
from cctbx.array_family import flex

from iotbx.reflection_file_editor import guess_array_output_labels


# Original is in iotbx.reflection_file_editor
def get_original_array_types (mtz_file, original_labels) :
    array_types = ""
    mtz_columns = mtz_file.column_labels()
    mtz_types = mtz_file.column_types()
    mtz_crossref = dict(zip(mtz_columns, mtz_types))
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

def aniso_res_cut_selection(arr, astarlim, bstarlim, cstarlim):
    """
    cctbx.miller.array.ellipsoidal_truncation_by_sigma()
    """

    ehm = 1./astarlim
    ekm = 1./bstarlim
    elm = 1./cstarlim
    selection = flex.bool(arr.indices().size(), False)
    #selection = flex.bool(arr.indices().size(), True)
    data = arr.data()
    sigmas = arr.sigmas()
    for i_mi, mi in enumerate(arr.indices()):
      rsv = arr.unit_cell().reciprocal_space_vector(mi)
      r = math.sqrt((rsv[0]/ehm)**2 + (rsv[1]/ekm)**2 + (rsv[2]/elm)**2)
      if(r<=1): selection[i_mi] = True
      #if(r>1 and data[i_mi]/sigmas[i_mi]<sigma_cutoff): selection[i_mi] = False
    return selection

def run(mtzin, rescut, mtzout):
    # Open mtz
    mtz_file = iotbx.mtz.object(mtzin)
    miller_arrays = mtz_file.as_miller_arrays()
    print "Opening", mtzin

    mtz_dataset = None
    labels = ["H", "K", "L"]

    for i, ar in enumerate(miller_arrays):
        d_spacings = ar.d_spacings().data()
        sel = aniso_res_cut_selection(ar, *rescut)

        print "%d reflections removed from %s" % (sum(~sel), ar.info().label_string())
        fake_label = 2 * string.uppercase[i]
        for lab in guess_array_output_labels(ar):
            labels.append(lab)
        array_types = get_original_array_types(mtz_file, ar.info().labels)
        default_types = iotbx.mtz.default_column_types(ar)
        if len(default_types) == len(array_types):
            column_types = array_types
        else:
            column_types = None

        mtz_dataset = add_array_to_mtz_dataset(mtz_dataset, ar.select(sel), fake_label,
                                               column_types)


    # Decide labels and write mtz file
    mtz_object = mtz_dataset.mtz_object()
    invalid_chars = re.compile("[^A-Za-z0-9_\-+\(\)]")

    used = dict([ (label, 0) for label in labels ])

    for i, column in enumerate(mtz_object.columns()):
        if column.label() != labels[i] :
            label = labels[i]
            original_label = label

            assert used[label] == 0

            try:
                column.set_label(label)
            except RuntimeError, e:
                if ("new_label is used already" in str(e)) :
                    col_names = [ col.label() for col in mtz_object.columns() ]
                    raise RuntimeError(("Duplicate column label '%s': current labels "+
                                        "are %s; user-specified output labels are %s.") %
                                       (label, " ".join(col_names), " ".join(labels)))
            else:
                used[original_label] += 1

    mtz_object.write(file_name=mtzout)

    print
    print "Writing:", mtzout
    print
# run()

if __name__ == "__main__":
    import sys

    mtzin = sys.argv[1]
    res_cut = map(float, sys.argv[2].split(",")) # resolution cutoff along a*,b*,c*
    mtzout = os.path.splitext(os.path.basename(mtzin))[0] + "_anisocutoff.mtz"

    run(mtzin= mtzin, rescut= res_cut, mtzout=mtzout)

