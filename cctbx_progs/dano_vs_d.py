"""
Usage:
 phenix.python dano_vs_d.py your.sca 20
"""
from __future__ import print_function
from __future__ import unicode_literals

import iotbx.file_reader
from cctbx.array_family import flex

def run(hklin, n_bins):
    for array in iotbx.file_reader.any_file(hklin).file_server.miller_arrays:
        # skip if not anomalous intensity data
        if not (array.is_xray_intensity_array() and array.anomalous_flag()):
            print("skipping", array.info())
            continue

        # We assume that data is already merged
        assert array.is_unique_set_under_symmetry()

        # take anomalous differences
        dano = array.anomalous_differences()

        # process with binning
        dano.setup_binner(n_bins=n_bins)
        binner = dano.binner()

        print("Array:", array.info())
        print("   dmax    dmin  nrefs  dano")
        for i_bin in binner.range_used():
            # selection for this bin. sel is flex.bool object (list of True of False)
            sel = binner.selection(i_bin)

            # take mean of absolute value of anomalous differences in a bin
            bin_mean = flex.mean(flex.abs(dano.select(sel).data()))
            d_max, d_min = binner.bin_d_range(i_bin)
            print("%7.2f %7.2f %6d %.2f" % (d_max, d_min, binner.count(i_bin), bin_mean))
# run()

if __name__ == "__main__":
    import sys
    hklin = sys.argv[1]
    n_bins = int(sys.argv[2])
    run(hklin, n_bins)
