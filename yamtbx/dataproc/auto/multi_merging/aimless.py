"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc import aimless
from yamtbx import util
import collections

import os

class AimlessCycles:
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
        print >>self.out, "\nIn %s" % new_wd

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
                    cell = map(float, l[l.index("=")+1:].split())
                    assert len(cell) == 6
                    cells.append(cell)
                    break

        cell_sum = reduce(lambda x,y:map(lambda a:x[a]+y[a], xrange(6)), cells)
        return sg, " ".join(map(lambda x:"%.3f"%(x/float(len(cells))), cell_sum))
    # average_cells()

    def run_cycles(self, xds_files):
        self.removed_files = []
        self.removed_reason = {}
        print >>self.out, "********************* START FUNCTION ***********************"

        # Remove unconnected files here.
        remove_idxes = find_unconnected_xds_files(xds_files, min_ios=3, d_min=self.d_min)
        print >>self.out, "DEBUG:: Need to remove unconnected %d files" % len(remove_idxes)
        for i in sorted(remove_idxes): 
            print >>self.out, " %.3d %s" % (i+1, xds_files[i])
            self.removed_files.append(xds_files[i])
            self.removed_reason[xds_files[i]] = "no_common_refls"

        keep_idxes = filter(lambda x: x not in remove_idxes, xrange(len(xds_files)))
        self.run_cycle(map(lambda i: xds_files[i], keep_idxes))

        # Final cycle; average cells and merge again.
        used_files = set(xds_files).difference(set(self.removed_files))
        sg, cell = self.average_cells(used_files)
        print >>self.out, "Final scaling with averaged cell."
        print >>self.out, "Averaged cell= %s" % (cell)

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
            print >>self.out, "Error: no files given."
            return

        inp_str = ""
        for i, f in enumerate(xds_files):
            brange = self.batch_info[f]
            inp_str += "RUN %3d BATCH %4d to %4d\n" % (i+1, brange[0], brange[1])

        print >>self.out, "DEBUG:: running aimless with %3d files.." % len(xds_files)
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
            print >>self.out, "Rejection based on delta_CC1/2 in %s shell" % self.delta_cchalf_bin
            table = aimless.read_summary(aimless_log)
            i_stat = 0 if self.delta_cchalf_bin == "total" else 2
            prev_cchalf = table["cc_half"][i_stat]
            prev_nuniq = table["nuniq"][i_stat]
            # file_name->idx table
            remaining_files = collections.OrderedDict(map(lambda x: x[::-1], enumerate(xds_files)))

            for i in xrange(len(xds_files)-1): # if only one file, cannot proceed.
                tmpdir = os.path.join(self.workdir, "reject_test_%.3d" % i)

                cchalf_list = aimless.calc_cchalf_by_removing(wdir=tmpdir,
                                                              mtzin=self.mtzin,
                                                              batch_info=self.batch_info,
                                                              inpfiles=remaining_files.keys(),
                                                              anomalous_flag=self.anomalous_flag,
                                                              d_min=self.d_min,
                                                              stat_bin=self.delta_cchalf_bin,
                                                              nproc=self.nproc,
                                                              nproc_each=self.nproc_each,
                                                              batchjobs=self.batchjobs)

                rem_idx, cc_i, nuniq_i = cchalf_list[0] # First (largest) is worst one to remove.
                rem_idx_in_org = remaining_files[remaining_files.keys()[rem_idx]]
                
                # Decision making by CC1/2
                print >>self.out, "DEBUG:: remove %3d if %.4f*%d > %.4f*%d" % (rem_idx_in_org, 
                                                                               cc_i, nuniq_i,
                                                                               prev_cchalf, prev_nuniq)
                if cc_i*nuniq_i <= prev_cchalf*prev_nuniq: break
                print >>self.out, "Removing idx= %3d gains CC1/2 by %.4f" % (rem_idx_in_org, cc_i-prev_cchalf)

                prev_cchalf, prev_nuniq = cc_i, nuniq_i
                remove_idxes.append(rem_idx_in_org)
                del remaining_files[remaining_files.keys()[rem_idx]] # remove file from table
        else:
            print >>self.out, "ERROR:: Unsupported reject_method (%s)" % reject_method

        if len(remove_idxes) > 0:
            print >>self.out, "DEBUG:: Need to remove %d files" % len(remove_idxes)
            for i in sorted(remove_idxes): 
                print >>self.out, " %.3d %s" % (i+1, xds_files[i])
                self.removed_files.append(xds_files[i])
                self.removed_reason[xds_files[i]] = "badcc"

        if self.next_delta_cchalf_bin != []:
            self.delta_cchalf_bin = self.next_delta_cchalf_bin.pop(0)
            do_rejection = True
        else:
            do_rejection = False

        if do_rejection or len(remove_idxes) > 0:
            keep_idxes = filter(lambda x: x not in remove_idxes, xrange(len(xds_files)))
            self.workdir = self.request_next_workdir()
            self.run_cycle(map(lambda i: xds_files[i], keep_idxes), do_rejection=do_rejection)
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

    for i in xrange(len(arrays)-1):
        for j in xrange(i+1, len(arrays)):
            matchs = arrays[i].match_indices(other=arrays[j], assert_is_similar_symmetry=False)
            if matchs.pairs().size() >= 10:
                G.add_edge(i, j)
                print "edge", i, j
    
    ccomps = map(lambda x:x, nx.connected_components(G))
    print "DEBUG:: Connected components=", ccomps
    keep_idxes = ccomps[0]
    remove_idxes = filter(lambda x: x not in keep_idxes, xrange(len(xds_files)))
    return remove_idxes
# find_unconnected_xds_files()
