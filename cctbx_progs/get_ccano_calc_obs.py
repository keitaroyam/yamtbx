"""
Computes anomalous CC(DFmodel, DFobs). Assuming user provides phenix.refine'd mtz file.
"""
from __future__ import print_function
from __future__ import unicode_literals
import iotbx.file_reader
from cctbx.array_family import flex

def ccano(hkl1, hkl2):
    hkl1, hkl2 = hkl1.as_intensity_array(), hkl2.as_intensity_array()

    corr = flex.linear_correlation(hkl1.anomalous_differences().data(),
                                   hkl2.anomalous_differences().data())
    assert corr.is_well_defined()
    return corr.coefficient()
# ccano()

def run(mtzin):
    arrays = iotbx.file_reader.any_file(mtzin).file_server.miller_arrays

    ano_arrays = [x for x in arrays if x.anomalous_flag()]
    print([x.info().labels for x in ano_arrays])

    f_model = [x for x in ano_arrays if "model" in x.info().label_string()][0]
    f_model = f_model.as_amplitude_array()
    f_obs = [x for x in ano_arrays if "obs" in x.info().label_string()][0]
    f_model, f_obs = f_model.common_sets(f_obs)

    binner = f_obs.setup_binner(n_bins=20)
    for i_bin in binner.range_used():
        sel = (binner.bin_indices() == i_bin)
        print("%6.3f - %6.3f" % binner.bin_d_range(i_bin), end=' ')
        print("%.4f" % ccano(f_obs.select(sel), f_model.select(sel)))

    print("overall: %.4f" % ccano(f_model, f_obs))



# x. x.info().label_string(), arrays)



if __name__ == "__main__":
    import sys
    run(sys.argv[1])
