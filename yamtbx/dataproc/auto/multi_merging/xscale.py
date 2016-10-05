"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.xds import xscale
from yamtbx.dataproc.xds import xscalelp
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.dataproc.xds.command_line import xds2mtz
from yamtbx.dataproc import blend_lcv
from yamtbx import util
import collections

import os
import glob
import traceback
import networkx as nx
import numpy

xscale_comm = "xscale_par"

class XscaleCycles:
    def __init__(self, workdir, anomalous_flag, d_min, d_max,
                 reject_method, reject_params, xscale_params,
                 reference_file, out, nproc=1, nproc_each=None, batchjobs=None):
        self.reference_file = None
        self._counter = 0
        self.workdir_org = workdir
        self.anomalous_flag = anomalous_flag
        self.d_min = d_min
        self.d_max = d_max
        self.reject_params = reject_params
        self.xscale_params = xscale_params
        self.reference_choice = xscale_params.reference if not reference_file else None
        self.out = out
        self.nproc = nproc
        self.nproc_each = nproc_each
        self.batchjobs = batchjobs
        self.all_data_root = None # the root directory for all data
        self.altfile = {} # Modified files
        self.cell_info_at_cycles = {}

        if reject_params.delta_cchalf.bin == "total-then-outer":
            self.delta_cchalf_bin = "total"
            self.next_delta_cchalf_bin = ["outer"]
        else:
            self.delta_cchalf_bin = reject_params.delta_cchalf.bin
            self.next_delta_cchalf_bin = []

        # Define order of rejections
        self.reject_method = []
        if "framecc" in reject_method: self.reject_method.append("framecc")
        if "lpstats" in reject_method: self.reject_method.append("lpstats")
        if "delta_cc1/2" in reject_method: self.reject_method.append("delta_cc1/2")

        self.workdir = self.request_next_workdir()

        if self.d_min is not None and self.d_max is None:
            self.d_max = 100

        self.xscale_inp_head = """\
MINIMUM_I/SIGMA= %.2f
OUTPUT_FILE= xscale.hkl
""" % xscale_params.min_i_over_sigma
        if self.anomalous_flag is not None:
            self.xscale_inp_head += "FRIEDEL'S_LAW= %s\n" % ("FALSE" if self.anomalous_flag else "TRUE")

        if self.d_min is not None:
            nbins = 9
            step = ( 1./(self.d_min**2) - 1./(self.d_max**2) ) / float(nbins)
            start = 1./(self.d_max**2)
            rshells = " ".join(map(lambda i: "%.2f" % (start + i * step)**(-1./2), xrange(1, nbins+1)))
            self.xscale_inp_head += "RESOLUTION_SHELLS= %s\n\n" % rshells

        if reference_file:
            self.reference_file = os.path.join(workdir, "reference.hkl")
            os.symlink(os.path.relpath(reference_file, workdir), self.reference_file)
    # __init__()

    def request_next_workdir(self):
        self._counter += 1
        new_wd = os.path.join(self.workdir_org, "run_%.2d" % self._counter)
        os.mkdir(new_wd)
        print >>self.out, "\nIn %s" % new_wd

        return new_wd
    # request_next_workdir()

    def get_last_cycle_number(self): return self._counter

    def request_file_modify(self, filename):
        relp = os.path.relpath(filename, self.all_data_root)
        newpath = os.path.join(self.workdir_org, "files", relp)
        newpathd = os.path.dirname(newpath)
        if not os.path.exists(newpathd): os.makedirs(newpathd)
        self.altfile[filename] = newpath
        return newpath
    # request_file_modify

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

        cells = numpy.array(cells)
        mean_cell = map(lambda i: cells[:,i].mean(), xrange(6))
        cell_std = map(lambda i: numpy.std(cells[:,i]), xrange(6))
        lcv, alcv = blend_lcv.calc_lcv(cells)

        mean_cell_str = " ".join(map(lambda x:"%.3f"%x, mean_cell))
        cell_std_str = " ".join(map(lambda x:"%.1e"%x, cell_std))

        print >>self.out, "Averaged cell= %s (%d)" % (mean_cell_str, len(cells))
        print >>self.out, "Cell Std dev.= %s" % cell_std_str
        print >>self.out, "LCV, aLCV= %.3f%%, %.3f A" % (lcv, alcv)

        return sg, mean_cell_str, lcv, alcv
    # average_cells()

    def run_cycles(self, xds_ascii_files):
        self.all_data_root = os.path.dirname(os.path.commonprefix(xds_ascii_files))
        self.removed_files = []
        self.removed_reason = {}
        print >>self.out, "********************* START FUNCTION ***********************"
        if self.reference_file:
            self.run_cycle([self.reference_file,]+xds_ascii_files)
        else:
            self.run_cycle(xds_ascii_files)

        for i in xrange(1, self.get_last_cycle_number()+1):
            wd = os.path.join(self.workdir_org, "run_%.2d"%i)
            try:
                xds2mtz.xds2mtz(xds_file=os.path.abspath(os.path.join(wd, "xscale.hkl")),
                                dir_name=os.path.join(wd, "ccp4"),
                                run_xtriage=True, run_ctruncate=True)
            except:
                # Don't want to stop the program.
                print >>self.out, traceback.format_exc()

        return self.removed_files, self.removed_reason
    # run_cycles()

    def check_remove_list(self, remove_idxes):
        new_list = []
        skip_num = 0
        if self.reference_file: skip_num += 1 # because first is reference file
        #if self._counter > 1: skip_num += 1 # because first (or second) is the previously merged file
        remove_idxes = filter(lambda x: x >= skip_num, remove_idxes)
        return remove_idxes
    # check_remove_list()
    
    def current_working_dir(self): return self.workdir

    def run_cycle(self, xds_ascii_files, reference_idx=None):
        if len(xds_ascii_files) == 0:
            print >>self.out, "Error: no files given."
            return

        xscale_inp = os.path.join(self.workdir, "XSCALE.INP")
        xscale_lp = os.path.join(self.workdir, "XSCALE.LP")

        # Get averaged cell for scaling
        sg, cell, lcv, alcv = self.average_cells(xds_ascii_files)
        self.cell_info_at_cycles[self.get_last_cycle_number()] = (cell, lcv, alcv)
        
        # Choose directory containing XDS_ASCII.HKL and set space group (but how??)
        inp_out = open(xscale_inp, "w")
        inp_out.write("MAXIMUM_NUMBER_OF_PROCESSORS= %d\n" % self.nproc)
        inp_out.write("SPACE_GROUP_NUMBER= %s\nUNIT_CELL_CONSTANTS= %s\n\n" % (sg, cell))
        inp_out.write(self.xscale_inp_head)

        for i, xds_ascii in enumerate(xds_ascii_files):
            f = self.altfile.get(xds_ascii, xds_ascii)
            tmp = min(os.path.relpath(f, self.workdir), f, key=lambda x:len(x))
            refstr = "*" if i==reference_idx else " "
            inp_out.write(" INPUT_FILE=%s%s\n" % (refstr,tmp))
            if len(self.xscale_params.corrections) != 3:
                inp_out.write("  CORRECTIONS= %s\n" % " ".join(self.xscale_params.corrections))
            if self.xscale_params.frames_per_batch is not None:
                frame_range = XDS_ASCII(f, read_data=False).get_frame_range()
                nframes = frame_range[1] - frame_range[0]
                nbatch = int(numpy.ceil(nframes / self.xscale_params.frames_per_batch))
                print >>self.out, "frame range of %s is %d,%d setting NBATCH= %d" % (f, frame_range[0], frame_range[1], nbatch)
                inp_out.write("  NBATCH= %d\n" % nbatch)

        inp_out.close()

        print >>self.out, "DEBUG:: running xscale with %3d files.." % len(xds_ascii_files)
        try:
            xscale.run_xscale(xscale_inp, cbf_to_dat=True,
                              use_tmpdir_if_available=self.xscale_params.use_tmpdir_if_available)
        except:
            print >>self.out, traceback.format_exc()

        xscale_log = open(xscale_lp).read()
        if "!!! ERROR !!! INSUFFICIENT NUMBER OF COMMON STRONG REFLECTIONS." in xscale_log:
            print >>self.out, "DEBUG:: Need to choose files."

            # From XDS ver. March 1, 2015, it kindly informs which dataset has no common reflections.
            # ..but does not print the table. Sometimes only one dataset is left. Should we make table by ourselves?
            # Older versions just print correlation table and stop.
            if "CORRELATIONS BETWEEN INPUT DATA SETS AFTER CORRECTIONS" in xscale_log:
                G = xscalelp.construct_data_graph(xscale_lp, min_common_refs=10)
                #nx.write_dot(G, os.path.join(self.workdir, "common_set_graph.dot"))
                cliques = [c for c in nx.find_cliques(G)]
                cliques.sort(key=lambda x:len(x))
                if self._counter == 1:
                    max_clique = cliques[-1]
                else:
                    idx_prevfile = 1 if self.reference_file else 0
                    max_clique = filter(lambda x: idx_prevfile in x, cliques)[-1] # xscale.hkl must be included!

                if self.reference_file:
                    max_clique = [0,] + filter(lambda x: x!=0, max_clique)

                for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(self.workdir, f))

                try_later = map(lambda i: xds_ascii_files[i], filter(lambda x: x not in max_clique, G.nodes()))

                print >>self.out, "DEBUG:: %d files can be merged. %d files will be merged later." % (len(max_clique),
                                                                                                      len(try_later))
                print >>self.out, "DEBUG:: %d files are of no use." % (len(xds_ascii_files)-len(G.nodes()))
                for i in filter(lambda j: j not in G.nodes(), xrange(len(xds_ascii_files))):
                    self.removed_files.append(xds_ascii_files[i])
                    self.removed_reason[xds_ascii_files[i]] = "no_common_refls"

                self.run_cycle(map(lambda i: xds_ascii_files[i], max_clique))

                assert len(try_later) <= 0 # Never be the case with newer xscale!! (if the case, check_remove_list() should be modified to skip_num+=1
                if len(try_later) > 0:
                    print >>self.out, "Trying to merge %d remaining files.." % len(try_later)
                    next_files = [os.path.join(self.workdir, "xscale.hkl")] + try_later
                    if self.reference_file: next_files = [self.reference_file,] + next_files
                    self.workdir = self.request_next_workdir()
                    self.run_cycle(next_files)
                    return
            else:
                bad_idxes = xscalelp.read_no_common_ref_datasets(xscale_lp)
                print >>self.out, "DEBUG:: %d files are of no use." % (len(bad_idxes))

                for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(self.workdir, f))

                # XXX Actually, not all datasets need to be thrown.. some of them are useful..
                for i in bad_idxes:
                    self.removed_files.append(xds_ascii_files[i])
                    self.removed_reason[xds_ascii_files[i]] = "no_common_refls"

                self.run_cycle(map(lambda i: xds_ascii_files[i], 
                                   filter(lambda j: j not in bad_idxes, xrange(len(xds_ascii_files)))))

            return
        elif "!!! ERROR !!! USELESS DATA ON INPUT REFLECTION FILE" in xscale_log:
            print >>self.out, "DEBUG:: Need to discard useless data."
            unuseful_data = [xscalelp.get_read_data(xscale_lp)[-1]] #filter(lambda x: x[2]==0, xscalelp.get_read_data(xscale_lp))
            if len(unuseful_data) == 0:
                print >>self.out, "I don't know how to fix it.."
                return
            remove_idxes = map(lambda x: x[0]-1, unuseful_data)
            remove_idxes = self.check_remove_list(remove_idxes)
            keep_idxes = filter(lambda x: x not in remove_idxes, xrange(len(xds_ascii_files)))
            for i in remove_idxes:
                self.removed_files.append(xds_ascii_files[i])
                self.removed_reason[xds_ascii_files[i]] = "useless"

            for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(self.workdir, f))
            self.run_cycle(map(lambda i: xds_ascii_files[i], keep_idxes))
            return
        elif "INACCURATE SCALING FACTORS." in xscale_log:
            # Actually I don't know how to fix this.. (bug?) but worth proceeding (discarding bad data may solve problem).
            print >>self.out, "'INACCURATE SCALING FACTORS' happened.. but ignored."
        elif "!!! ERROR !!!" in xscale_log:
            print >>self.out, "Unknown error! please check the XSCALE.LP and fix the program."
            return

        # Re-scale by changing reference
        rescale_for = None
        if len(self.reject_method) == 0:
            rescale_for = self.reference_choice # may be None
        elif reference_idx is None:
            rescale_for = "bmed"
        
        if rescale_for is not None and len(xds_ascii_files) > 1:
            ref_num = xscale.decide_scaling_reference_based_on_bfactor(xscale_lp, rescale_for, return_as="index")
            if reference_idx != ref_num:
                print >>self.out, "Rescaling with %s" % rescale_for
                for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(self.workdir, f))
                self.run_cycle(xds_ascii_files, reference_idx=ref_num)

        if len(self.reject_method) == 0:
            return

        # Remove bad data
        remove_idxes = []
        remove_reasons = {}

        if self.reject_method[0] == "framecc":
            print >>self.out, "Rejections based on frame CC"
            from yamtbx.dataproc.xds.command_line import xscale_cc_against_merged

            # list of [frame, n_all, n_common, cc] in the same order
            framecc = xscale_cc_against_merged.run(hklin=os.path.join(self.workdir, "xscale.hkl"),
                                                   output_dir=self.workdir,
                                                   nproc=self.nproc).values()
            if self.reject_params.framecc.method == "tukey":
                ccs = numpy.array(map(lambda x: x[3], reduce(lambda x,y:x+y,framecc)))
                ccs = ccs[ccs==ccs] # Remove nan
                q25, q75 = numpy.percentile(ccs, [25, 75])
                cc_cutoff  = q25 - self.reject_params.framecc.iqr_coeff * (q75 - q25)
                print >>self.out, " frameCC cutoff = %.4f (%.2f*IQR)" % (cc_cutoff, self.reject_params.framecc.iqr_coeff)
            else:
                cc_cutoff = self.reject_params.framecc.abs_cutoff
                print >>self.out, " frameCC cutoff = %.4f (value specified)" % cc_cutoff

            for i, cclist in enumerate(framecc):
                useframes = map(lambda x: x[0], filter(lambda x: x[3] > cc_cutoff, cclist))
                if len(useframes) == 0:
                    remove_idxes.append(i)
                    remove_reasons.setdefault(i, []).append("allbadframe")
                    continue

                f = xds_ascii_files[i]
                xac = XDS_ASCII(f)
                if set(useframes).issuperset(set(range(min(xac.iframe), max(xac.iframe)))):
                    continue # All useful frames.

                sel = xac.iframe == useframes[0]
                for x in useframes[1:]: sel |= xac.iframe == x
                if sum(sel) < 10: # XXX care I/sigma
                    remove_idxes.append(i)
                    remove_reasons.setdefault(i, []).append("allbadframe")
                    continue

                print >>self.out, "Extracting frames %s out of %d-%d in %s" % (",".join(map(str,useframes)),
                                                                               min(xac.iframe), max(xac.iframe),
                                                                               f)

                newf = self.request_file_modify(f)
                xac.write_selected(sel, newf)

            self.reject_method.pop(0) # Perform only once

        elif self.reject_method[0] == "lpstats":
            if "bfactor" in self.reject_params.lpstats.stats:
                iqrc = self.reject_params.lpstats.iqr_coeff
                print >>self.out, "Rejections based on B-factor outliers (%.2f*IQR)" % iqrc
                Bs = numpy.array(map(lambda x:x[1], xscalelp.get_k_b(xscale_lp)))
                q25, q75 = numpy.percentile(Bs, [25, 75])
                iqr = q75 - q25
                lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
                count = 0
                for i, b in enumerate(Bs):
                    if b < lowlim or b > highlim:
                        remove_idxes.append(i)
                        remove_reasons.setdefault(i, []).append("bad_B")
                        count += 1

                print >>self.out, " %4d B-factor outliers (<%.2f, >%.2f) removed"% (count, lowlim, highlim)

            if "em.b" in self.reject_params.lpstats.stats:
                iqrc = self.reject_params.lpstats.iqr_coeff
                print >>self.out, "Rejections based on error model b outliers (%.2f*IQR)" % iqrc
                bs = numpy.array(map(lambda x:x[1], xscalelp.get_ISa(xscale_lp)))
                q25, q75 = numpy.percentile(bs, [25, 75])
                iqr = q75 - q25
                lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
                count = 0
                for i, b in enumerate(bs):
                    if b < lowlim or b > highlim:
                        remove_idxes.append(i)
                        remove_reasons.setdefault(i, []).append("bad_em.b")
                        count += 1

                print >>self.out, " %4d error model b outliers (<%.2f, >%.2f) removed"% (count, lowlim, highlim)

            if "em.ab" in self.reject_params.lpstats.stats:
                iqrc = self.reject_params.lpstats.iqr_coeff
                print >>self.out, "Rejections based on error model a*b outliers (%.2f*IQR)" % iqrc
                vals = numpy.array(map(lambda x:x[0]*x[1], xscalelp.get_ISa(xscale_lp)))
                q25, q75 = numpy.percentile(vals, [25, 75])
                iqr = q75 - q25
                lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
                count = 0
                for i, ab in enumerate(vals):
                    if ab < lowlim or ab > highlim:
                        remove_idxes.append(i)
                        remove_reasons.setdefault(i, []).append("bad_em.ab")
                        count += 1

                print >>self.out, " %4d error model a*b outliers (<%.2f, >%.2f) removed"% (count, lowlim, highlim)

            if "rfactor" in self.reject_params.lpstats.stats:
                iqrc = self.reject_params.lpstats.iqr_coeff
                print >>self.out, "Rejections based on R-factor outliers (%.2f*IQR)" % iqrc
                rstats = xscalelp.get_rfactors_for_each(xscale_lp)
                vals = numpy.array(map(lambda x:rstats[x][-1][1], rstats)) # Read total R-factor
                q25, q75 = numpy.percentile(vals, [25, 75])
                iqr = q75 - q25
                lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
                count = 0
                for i, v in enumerate(vals):
                    if v < lowlim or v > highlim:
                        remove_idxes.append(i)
                        remove_reasons.setdefault(i, []).append("bad_R")
                        count += 1

                print >>self.out, " %4d R-factor outliers (<%.2f, >%.2f) removed"% (count, lowlim, highlim)

            if "pairwise_cc" in self.reject_params.lpstats.stats:
                corrs = xscalelp.get_pairwise_correlations(xscale_lp)
                if self.reject_params.lpstats.pwcc.method == "tukey":
                    q25, q75 = numpy.percentile(map(lambda x: x[3], corrs), [25, 75])
                    iqr = q75 - q25
                    lowlim = q25 - self.reject_params.lpstats.pwcc.iqr_coeff * iqr
                    print >>self.out, "Rejections based on pairwise_cc < %.4f (IQR=%.2f)" % (lowlim, iqr)
                else:
                    lowlim = self.reject_params.lpstats.pwcc.abs_cutoff
                    print >>self.out, "Rejections based on pairwise_cc < %.4f" % lowlim

                bad_corrs = filter(lambda x: x[3] < lowlim, corrs)
                idx_bad = {}
                for i, j, common_refs, corr, ratio, bfac in bad_corrs:
                    idx_bad[i] = idx_bad.get(i, 0) + 1
                    idx_bad[j] = idx_bad.get(j, 0) + 1

                idx_bad = idx_bad.items()
                idx_bad.sort(key=lambda x:x[1])
                count = 0
                for idx, badcount in reversed(idx_bad):
                    remove_idxes.append(idx-1)
                    remove_reasons.setdefault(idx-1, []).append("bad_pwcc")
                    bad_corrs = filter(lambda x: idx not in x[:2], bad_corrs)
                    if len(bad_corrs) == 0: break
                    fun_key = lambda x: x[3]
                    print >>self.out, " Removing idx=%d (CC %.3f..%.3f) remaining %d bad pairs" % (idx, 
                                                                                                   min(bad_corrs,key=fun_key)[3],
                                                                                                   max(bad_corrs,key=fun_key)[3],
                                                                                                   len(bad_corrs))
                    count += 1
                print >>self.out, " %4d pairwise CC outliers removed" % count

            self.reject_method.pop(0) # Perform only once
        elif self.reject_method[0] == "delta_cc1/2":
            print >>self.out, "Rejection based on delta_CC1/2 in %s shell" % self.delta_cchalf_bin
            table = xscalelp.read_stats_table(xscale_lp)
            i_stat = -1 if self.delta_cchalf_bin == "total" else -2
            prev_cchalf = table["cc_half"][i_stat]
            prev_nuniq = table["nuniq"][i_stat]
            # file_name->idx table
            remaining_files = collections.OrderedDict(map(lambda x: x[::-1], enumerate(xds_ascii_files)))

            # For consistent resolution limit
            inp_head = self.xscale_inp_head + "SPACE_GROUP_NUMBER= %s\nUNIT_CELL_CONSTANTS= %s\n\n" % (sg, cell)
            count = 0
            for i in xrange(len(xds_ascii_files)-1): # if only one file, cannot proceed.
                tmpdir = os.path.join(self.workdir, "reject_test_%.3d" % i)

                cchalf_list = xscale.calc_cchalf_by_removing(wdir=tmpdir, inp_head=inp_head,
                                                             inpfiles=remaining_files.keys(),
                                                             stat_bin=self.delta_cchalf_bin,
                                                             nproc=self.nproc,
                                                             nproc_each=self.nproc_each,
                                                             batchjobs=self.batchjobs)

                rem_idx, cc_i, nuniq_i = cchalf_list[0] # First (largest) is worst one to remove.
                rem_idx_in_org = remaining_files[remaining_files.keys()[rem_idx]]
                
                # Decision making by CC1/2
                print >>self.out, "DEBUG:: cycle %.3d remove %3d if %.2f*%d > %.2f*%d" % (i, rem_idx_in_org, 
                                                                                          cc_i, nuniq_i,
                                                                                          prev_cchalf, prev_nuniq)
                if cc_i*nuniq_i <= prev_cchalf*prev_nuniq: break
                print >>self.out, "Removing idx= %3d gained CC1/2 by %.2f" % (rem_idx_in_org, cc_i-prev_cchalf)

                prev_cchalf, prev_nuniq = cc_i, nuniq_i
                remove_idxes.append(rem_idx_in_org)
                remove_reasons.setdefault(rem_idx_in_org, []).append("bad_cchalf")
                del remaining_files[remaining_files.keys()[rem_idx]] # remove file from table
                count += 1

            print >>self.out, " %4d removed by DeltaCC1/2 method" % count

            if self.next_delta_cchalf_bin != []:
                self.delta_cchalf_bin = self.next_delta_cchalf_bin.pop(0)
            else:
                self.reject_method.pop(0)
        else:
            print >>self.out, "ERROR:: Unsupported reject_method (%s)" % reject_method

        # Remove duplicates
        remove_idxes = list(set(remove_idxes))
        remove_idxes = self.check_remove_list(remove_idxes)
        if len(remove_idxes) > 0:
            print >>self.out, "DEBUG:: Need to remove %d files" % len(remove_idxes)
            for i in sorted(remove_idxes): 
                print >>self.out, " %.3d %s" % (i, xds_ascii_files[i])
                self.removed_files.append(xds_ascii_files[i])
                self.removed_reason[xds_ascii_files[i]] = ",".join(remove_reasons[i])

        # Next run
        keep_idxes = filter(lambda x: x not in remove_idxes, xrange(len(xds_ascii_files)))
        if len(self.reject_method) > 0 or len(remove_idxes) > 0:
            self.workdir = self.request_next_workdir()
            self.run_cycle(map(lambda i: xds_ascii_files[i], keep_idxes))
        elif self.reference_choice is not None and len(keep_idxes) > 1:
            # Just re-scale with B reference
            ref_num = xscale.decide_scaling_reference_based_on_bfactor(xscale_lp, self.reference_choice, return_as="index")
            if reference_idx != ref_num:
                print >>self.out, "Rescaling2 with %s" % self.reference_choice
                for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(self.workdir, f))
                self.run_cycle(map(lambda i: xds_ascii_files[i], keep_idxes), reference_idx=ref_num)

    # run_cycle()
# class XscaleCycles
