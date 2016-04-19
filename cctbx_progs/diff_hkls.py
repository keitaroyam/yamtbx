"""
Give two mtz files (mtz_1 and mtz_2).
Reflections only in mtz_1 will be written as only_in_1.mtz. only_in_2.mtz will be written in the same way.
"""

import mmtbx.utils
import iotbx.phil
import iotbx.reflection_file_utils
import iotbx.mtz
from cctbx import miller

master_params_str = """\
hklin_1 = None
    .type=str
    .input_size = 160
    .short_caption = file name of mtz 1
hklin_2 = None
    .type=str
    .input_size = 160
    .short_caption = file name of mtz 2
labels_1 = None
    .type=str
    .input_size = 160
    .short_caption = Data labels for mtz 1
labels_2 = None
    .type=str
    .input_size = 160
    .short_caption = Data labels for mtz 2
"""

def get_data(mtzin, labels):
    mtzobj = iotbx.mtz.object(file_name=mtzin)
    selected = filter(lambda x:x.info().label_string() == labels, mtzobj.as_miller_arrays())

    if len(selected) < 1:
        print mtzin, "does not have", labels
        print " Possible labels:", [x.info().label_string() for x in mtzobj.as_miller_arrays()]

    return selected[0]
# get_data()

if __name__ == "__main__":
    import sys

    parsed = iotbx.phil.parse(master_params_str, process_includes=True)


    processed_args = mmtbx.utils.process_command_line_args(args = sys.argv[1:],
                                                           log = sys.stdout,
                                                           master_params = parsed)

    working_phil = processed_args.params
    params = working_phil.extract()

    if params.hklin_1 is None and params.hklin_2 is None:
        if len(processed_args.reflection_file_names) != 2:
            print "Exactly two mtz files must be given."
            sys.exit(1)
        params.hklin_1, params.hklin_2 = processed_args.reflection_file_names

    working_phil = parsed.format(python_object=params)
    print "Parameters to compute maps:"
    working_phil.show(out = sys.stdout, prefix=" ")


    data_1 = get_data(params.hklin_1, params.labels_1)
    data_2 = get_data(params.hklin_2, params.labels_2)

    matches = miller.match_indices(data_1.indices(), data_2.indices())

    only_1 = data_1.select(matches.singles(0))
    only_2 = data_2.select(matches.singles(1))

    only_1.as_mtz_dataset(column_root_label=params.labels_1.split(",")[0]).mtz_object().write("only_in_1.mtz")
    only_2.as_mtz_dataset(column_root_label=params.labels_2.split(",")[0]).mtz_object().write("only_in_2.mtz")
