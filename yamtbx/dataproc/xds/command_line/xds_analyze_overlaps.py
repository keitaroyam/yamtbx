#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

"""
Analyze overlapped spots due to non-merohedral twinning or multiple crystals.
Give two XDS_ASCII.HKL files to this script.


"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from cctbx import miller
from cctbx.array_family import flex
from iotbx.merging_statistics import dataset_statistics, filter_intensities_by_sigma
from yamtbx.dataproc.xds import xds_ascii
from scipy import spatial
import numpy
import sys
import os

def run(files):
    assert len(files) == 2

    hkl1 = xds_ascii.XDS_ASCII(files[0], sys.stdout)
    hkl2 = xds_ascii.XDS_ASCII(files[1], sys.stdout)

    hkl1_points = numpy.column_stack((hkl1.xd, hkl1.yd, hkl1.zd))
    tree1 = spatial.cKDTree(hkl1_points)

    n_ovl, n_nonovl = 0, 0

    novl_indices, novl_i, novl_sigma = flex.miller_index(), flex.double(), flex.double()

    for i in range(len(hkl2.indices)):
        x, y, z = hkl2.xd[i], hkl2.yd[i], hkl2.zd[i]
        #if z > 180:
        #    continue
        dists, idxs = tree1.query((x,y,z), k=3, p=1)

        overlaps = []
        for dist, idx in zip(dists, idxs):
            idx = int(idx)
            xo, yo, zo = hkl1.xd[idx], hkl1.yd[idx], hkl1.zd[idx]

            if abs(z-zo) < 2.5 and (xo-x)**2+(yo-y)**2 < 15**2: # FIXME MAGIC NUMBER!
                overlaps.append((dist,idx))

        if len(overlaps) == 0:
            novl_indices.append(hkl2.indices[i])
            novl_i.append(hkl2.iobs[i])
            novl_sigma.append(hkl2.sigma_iobs[i])
            n_nonovl += 1
        else:
            print(hkl2.indices[i], x, y, z)
            for dist, idx in overlaps:
                xo, yo, zo = hkl1.xd[idx], hkl1.yd[idx], hkl1.zd[idx]
                print(hkl1.indices[idx], xo, yo, zo)
                print(dist, idx)
            print()

    print()
    n_ref = len(hkl2.indices)
    print("%.2f%% overlap!" % (100.*(n_ref-n_nonovl)/n_ref))

    novl_array = miller.array(miller_set=miller.set(crystal_symmetry=hkl2.symm, indices=novl_indices),
                              data=novl_i, sigmas=novl_sigma)

    stats = dataset_statistics(novl_array, anomalous=False, sigma_filtering="xds")
    stats.show(out=sys.stdout)

    novl_array = novl_array.customized_copy(anomalous_flag=False).map_to_asu()
    novl_array = novl_array.eliminate_sys_absent()
    novl_array = novl_array.select(novl_array.sigmas() >= 0)

    filtr = filter_intensities_by_sigma(novl_array, "xds")
    hklout = os.path.splitext(os.path.basename(files[1]))[0] + "_novl.mtz"
    filtr.array_merged.set_observation_type_xray_intensity().as_mtz_dataset(column_root_label="IMEAN").mtz_object().write(hklout)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: %s XDS_ASCII.HKL.1 XDS_ASCII.HKL.2" % sys.argv[0])
        quit()

    run(sys.argv[1:])
