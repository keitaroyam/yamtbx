import iotbx.file_reader
from cctbx.array_family import flex
import iotbx.phil

master_params_str = """
d_min = None
 .type = float
d_max = None
 .type = float
label = None
 .type = str
"""

def run(mtzin, params):
    arrays = iotbx.file_reader.any_file(mtzin).file_server.miller_arrays

    if len(arrays) == 1:
        array = arrays[0]
    else:
        sel = filter(lambda x: x.info().label_string()==params.label, arrays)
        if len(sel) != 1:
            print "Possible choices for label=:"
            for a in arrays:
                print " %s" % a.info().label_string()
            return
        else:
            array = sel[0]

    array = array.resolution_filter(d_min=params.d_min, d_max=params.d_max)
    
    print flex.mean(array.data().as_double())

if __name__ == "__main__":
    import sys
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    run(args[0], params)
