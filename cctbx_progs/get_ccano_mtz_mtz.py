import iotbx.file_reader
from cctbx.array_family import flex
import iotbx.phil

master_params_str = """\
d_min = None
 .type = float
d_max = None
 .type = float
"""

def get_ccano(a1, a2):
    a1, a2 = a1.anomalous_differences(), a2.anomalous_differences()
    a1, a2 = a1.common_sets(a2)
    corr = flex.linear_correlation(a1.data(), a2.data())
    assert corr.is_well_defined()
    return corr.coefficient()
# get_ccano()

def remove_phase(arrays):
    for i in xrange(len(arrays)):
        if arrays[i].is_complex_array():
            info = arrays[i].info()
            arrays[i] = arrays[i].amplitudes().set_observation_type_xray_amplitude()
            arrays[i] = arrays[i].set_info(info)

def run(hklin1, hklin2, params):
    arrays1 = filter(lambda x:x.anomalous_flag(), iotbx.file_reader.any_file(hklin1).file_server.miller_arrays)
    arrays2 = filter(lambda x:x.anomalous_flag(), iotbx.file_reader.any_file(hklin2).file_server.miller_arrays)
    remove_phase(arrays1)
    remove_phase(arrays2)

    arrays1 = map(lambda x:x.resolution_filter(d_min=params.d_min, d_max=params.d_max).set_info(x.info()), arrays1)
    arrays2 = map(lambda x:x.resolution_filter(d_min=params.d_min, d_max=params.d_max).set_info(x.info()), arrays2)

    print "array_1: %s" % hklin1
    for a in arrays1:
        print "", a.info().label_string(), a.observation_type()
    print
    print "array_2: %s" % hklin2
    for a in arrays2:
        print "", a.info().label_string(), a.observation_type()
    print

    for a1 in arrays1:
        for a2 in arrays2:
            is_sq1, is_sq2 = False, False
            if a1.is_xray_intensity_array() and a2.is_xray_amplitude_array():
                a2 = a2.as_intensity_array().set_info(a2.info())
                is_sq2 = True
            elif a2.is_xray_intensity_array() and a1.is_xray_amplitude_array():
                a1 = a1.as_intensity_array().set_info(a1.info())
                is_sq1 = True

            try:
                ccano = get_ccano(a1, a2)
            except:
                ccano = None
            if ccano is not None:
                print a1.info(), "^2" if is_sq1 else "", "(%.2f, %.2f)"%a1.d_max_min()
                print a2.info(), "^2" if is_sq2 else "",  "(%.2f, %.2f)"%a2.d_max_min()
                print "CCano=", ccano
                print

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    hkl1, hkl2 = args
    run(hkl1, hkl2, params)
