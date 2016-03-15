master_params_str = """\
output = *mtz *sca
 .type = choice(multi=True)
prefix = None
 .type = str
 .help = output prefix
dmin = None
 .type = float
dmax = None
 .type = float
"""

import os
from yamtbx.dataproc import crystfel
import iotbx.scalepack.merge
import iotbx.mtz
import iotbx.phil

def run(args):
    
    cmdline = iotbx.phil.process_command_line(args=args,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    symm_source, hklin = cmdline.remaining_args

    hklfile = crystfel.hkl.HKLfile(symm_source=symm_source,
                                   hklin=hklin)

    hklfile.set_resolution(d_min=params.dmin, d_max=params.dmax)

    if params.prefix is None:
        params.prefix = os.path.splitext(os.path.basename(hklin))[0]

    if "sca" in params.output:
        iotbx.scalepack.merge.write(file_name="%s.sca" % (params.prefix),
                                    miller_array=hklfile.array)

    if "mtz" in params.output:
        mtz_dataset = hklfile.array.as_mtz_dataset(column_root_label="I")
        mtz_dataset.add_miller_array(miller_array=hklfile.redundancies, column_root_label="MULT")
        mtz_object = mtz_dataset.mtz_object()
        mtz_object.write(file_name="%s.mtz" % (params.prefix))
# run()

if __name__ == "__main__":
    import sys
    run(sys.argv[1:])
