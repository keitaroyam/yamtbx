"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from yamtbx.dataproc import aimless
from yamtbx import util
import collections

import os
from functools import reduce

class AimlessCycles(object):
    def __init__(self, workdir, anomalous_flag, d_min, d_max,
                 reject_method, cc_cutoff, delta_cchalf_bin,
                 mtzin, batch_info, out, nproc=1, nproc_each=None, batchjobs=None):
        self._counter = 0
        self.workdir_org = workdir
        self.anomalous_flag = anomalous_flag
        self.d_min = d_min
        self.d_max = d_max
        self.cc_cutoff = cc_cutoff
        self.reject_method = reject_method
        self.mtzin = mtzin # by pointless
        self.batch_info = batch_info # xdsfile, batch_range
        self.out = out
        self.nproc = nproc
        self.nproc_each = nproc_each
        self.batchjobs = batchjobs

        if delta_cchalf_bin == "total-then-outer":
            self.delta_cchalf_bin = "total"
            self.next_delta_cchalf_bin = ["outer"]
        else:
            self.delta_cchalf_bin = delta_cchalf_bin
            self.next_delta_cchalf_bin = []
            
        self.workdir = self.request_next_workdir()

        if self.d_min is not None and self.d_max is None:
            self.d_max = 100

    # __init__()

    def request_next_workdir(self):
        self._counter += 1
        new_wd = os.path.join(self.workdir_org, "run_%.2d" % self._counter)
        os.mkdir(new_wd)
        print("\nIn %s" % new_wd, file=self.out)

        return new_wd
    # request_next_workdir()

    def get_last_cycle_number(self): return self._counter

    def average_cells(self, files):
        cells = []
        sg = None
        for f in files:
            if not os.path.isfile(f):
                continue
            for l in open(f):
                if l.startswith("!SPACE_GROUP_NUMBER="):
                    sg = l[l.index("=")+1:].strip()
                if l.startswith("!UNIT_CELL_CONSTANTS="):
                    cell = list(map(float, l[l.index("=")+1:].split()))
                    assert len(cell) == 6
                    cells.append(cell)
                    break

        cell_sum = reduce(lambda x,y:[x[a]+y[a] for a in range(6)], cells)
        return sg, " ".join(["%.3f"%(x/float(len(cells))) for x in cell_sum])
    # average_cells()

    def run_cycles(self, xds_files):
        self.removed_files = []
        self.removed_reason = {}
        print("********************* START FUNCTION ***********************", file=self.out)

        # Remove unconnected files here.
        remove_idxes = find_unconnected_xds_files(xds_files, min_ios=3, d_min=self.d_min)
        print("DEBUG:: Need to remove unconnected %d files" % len(remove_idxes), file=self.out)
        for i in sorted(remove_idxes): 
            print(" %.3d %s" % (i+1, xds_files[i]), file=self.out)
            self.removed_files.append(xds_files[i])
            self.removed_reason[xds_files[i]] = "no_common_refls"

        keep_idxes = [x for x in range(len(xds_files)) if x not in remove_idxes]
        self.run_cycle([xds_files[i] for i in keep_idxes])

        # Final cycle; average cells and merge again.
        used_files = set(xds_files).difference(set(self.removed_files))
        sg, cell = self.average_cells(used_files)
        print("Final scaling with averaged cell.", file=self.out)
        print("Averaged cell= %s" % (cell), file=self.out)

        # call mtzutils to change cell!
        mtz_finalcell = "pointless_finalcell.mtz"
        util.call(cmd="mtzutils hklin %s hklout %s" % (os.path.relpath(self.mtzin, self.workdir),
                                                       mtz_finalcell),
                  wdir=self.workdir,
                  expects_in=[os.path.relpath(self.mtzin, self.workdir)],
                  expects_out=[mtz_finalcell],
                  stdin="cell %s\nend\n"%cell,
                  stdout=open(os.path.join(self.workdir, "mtzutils.log"), "w")
                  )

        inp_str = ""
        for i, f in enumerate(used_files):
            brange = self.batch_info[f]
            inp_str += "RUN %3d BATCH %4d to %4d\n" % (i+1, brange[0], brange[1])

        aimless.run_aimless(mtzin=mtz_finalcell,
                            wdir=self.workdir,
                            anomalous=self.anomalous_flag, d_min=self.d_min, prefix="aimless",
                            add_stdin=inp_str)

        ctruncate_arg = "-hklin aimless.mtz -hklout ctruncate.mtz -colin '/*/*/[IMEAN,SIGIMEAN]'"
        if self.anomalous_flag: ctruncate_arg += " -colano '/*/*/[I(+),SIGI(+),I(-),SIGI(-)]'"

        util.call(cmd="ctruncate %s" % ctruncate_arg,
                  wdir=self.workdir,
                  expects_in=["aimless.mtz"],
                  expects_out=["ctruncate.mtz"],
                  stdout=open(os.path.join(self.workdir, "ctruncate.log"), "w")
                  )

        return self.removed_files, self.removed_reason
    # run_cycles()
   
    def current_working_dir(self): return self.workdir

    def run_cycle(self, xds_files, do_rejection=True):
        if len(xds_files) == 0:
            print("Error: no files given.", file=self.out)
            return

        inp_str = ""
        for i, f in enumerate(xds_files):
            brange = self.batch_info[f]
            inp_str += "RUN %3d BATCH %4d to %4d\n" % (i+1, brange[0], brange[1])

        print("DEBUG:: running aimless with %3d files.." % len(xds_files), file=self.out)
        aimless.run_aimless(mtzin=os.path.relpath(self.mtzin, self.workdir),
                            wdir=self.workdir,
                            anomalous=self.anomalous_flag, d_min=self.d_min, prefix="aimless",
                            add_stdin=inp_str)
        aimless_log = os.path.join(self.workdir, "aimless.log")

        # XXX Aimless error handling here.

        if not do_rejection:
            return

        # Remove bad data
        remove_idxes = []

        if self.reject_method == "delta_cc1/2":
            print("Rejection based on delta_CC1/2 in %s shell" % self.delta_cchalf_bin, file=self.out)
            table = aimless.read_summary(aimless_log)
            i_stat = 0 if self.delta_cchalf_bin == "total" else 2
            prev_cchalf = table["cc_half"][i_stat]
            prev_nuniq = table["nuniq"][i_stat]
            # file_name->idx table
            remaining_files = collections.OrderedDict([x[::-1] for x in enumerate(xds_files)])

            for i in range(len(xds_files)-1): # if only one file, cannot proceed.
                tmpdir = os.path.join(self.workdir, "reject_test_%.3d" % i)

                cchalf_list = aimless.calc_cchalf_by_removing(wdir=tmpdir,
                                                              mtzin=self.mtzin,
                                                              batch_info=self.batch_info,
                                                              inpfiles=list(remaining_files.keys()),
                                                              anomalous_flag=self.anomalous_flag,
                                                              d_min=self.d_min,
                                                              stat_bin=self.delta_cchalf_bin,
                                                              nproc=self.nproc,
                                                              nproc_each=self.nproc_each,
                                                              batchjobs=self.batchjobs)

                rem_idx, cc_i, nuniq_i = cchalf_list[0] # First (largest) is worst one to remove.
                rem_idx_in_org = remaining_files[list(remaining_files.keys())[rem_idx]]
                
                # Decision making by CC1/2
                print("DEBUG:: remove %3d if %.4f*%d > %.4f*%d" % (rem_idx_in_org, 
                                                                               cc_i, nuniq_i,
                                                                               prev_cchalf, prev_nuniq), file=self.out)
                if cc_i*nuniq_i <= prev_cchalf*prev_nuniq: break
                print("Removing idx= %3d gains CC1/2 by %.4f" % (rem_idx_in_org, cc_i-prev_cchalf), file=self.out)

                prev_cchalf, prev_nuniq = cc_i, nuniq_i
                remove_idxes.append(rem_idx_in_org)
                del remaining_files[list(remaining_files.keys())[rem_idx]] # remove file from table
        else:
            print("ERROR:: Unsupported reject_method (%s)" % reject_method, file=self.out)

        if len(remove_idxes) > 0:
            print("DEBUG:: Need to remove %d files" % len(remove_idxes), file=self.out)
            for i in sorted(remove_idxes): 
                print(" %.3d %s" % (i+1, xds_files[i]), file=self.out)
                self.removed_files.append(xds_files[i])
                self.removed_reason[xds_files[i]] = "badcc"

        if self.next_delta_cchalf_bin != []:
            self.delta_cchalf_bin = self.next_delta_cchalf_bin.pop(0)
            do_rejection = True
        else:
            do_rejection = False

        if do_rejection or len(remove_idxes) > 0:
            keep_idxes = [x for x in range(len(xds_files)) if x not in remove_idxes]
            self.workdir = self.request_next_workdir()
            self.run_cycle([xds_files[i] for i in keep_idxes], do_rejection=do_rejection)
            # XXX second rejection trials are (sometimes?) waste.
    # run_cycle()
# class AimlessCycles

def find_unconnected_xds_files(xds_files, min_ios=3, d_min=None):
    import networkx as nx
    from yamtbx.dataproc.xds import xds_ascii
    from yamtbx.dataproc.xds import integrate_hkl_as_flex

    if len(xds_files) < 2:
        return []

    G = nx.Graph()
    
    arrays = []
    for f in xds_files:
        if xds_ascii.is_xds_ascii(f):
            a = xds_ascii.XDS_ASCII(f, i_only=True).i_obs().resolution_filter(d_min=d_min)
        elif integrate_hkl_as_flex.is_integrate_hkl(f):
            a = integrate_hkl_as_flex.reader(f, ["IOBS","SIGMA"]).i_obs().resolution_filter(d_min=d_min)
        else:
            raise "Never reaches here"
            
        a = a.select(a.sigmas()>0)
        a = a.merge_equivalents(use_internal_variance=False).array()
        arrays.append(a.select(a.data()/a.sigmas()>=min_ios))

    for i in range(len(arrays)-1):
        for j in range(i+1, len(arrays)):
            matchs = arrays[i].match_indices(other=arrays[j], assert_is_similar_symmetry=False)
            if matchs.pairs().size() >= 10:
                G.add_edge(i, j)
                print("edge", i, j)
    
    ccomps = [x for x in nx.connected_components(G)]
    print("DEBUG:: Connected components=", ccomps)
    keep_idxes = ccomps[0]
    remove_idxes = [x for x in range(len(xds_files)) if x not in keep_idxes]
    return remove_idxes
# find_unconnected_xds_files()
