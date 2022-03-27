"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from yamtbx.dataproc.xds import xscale
from yamtbx.dataproc.xds import xscalelp
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.dataproc.xds.command_line import xds2mtz
from yamtbx.dataproc.xds import modify_xdsinp, check_xds_version
from yamtbx.dataproc.pointless import Pointless
from yamtbx.dataproc import blend_lcv
from yamtbx.dataproc.auto.resolution_cutoff import estimate_resolution_based_on_cc_half, initial_estimate_byfit_cchalf
from yamtbx import util
from yamtbx.util import batchjob

import collections
import shutil
import os
import glob
import traceback
import networkx as nx
import numpy
from functools import reduce

xscale_comm = "xscale_par"

def make_bin_str(d_min, d_max, nbins=9):
    start = 1./(d_max**2) if d_max is not None else 0.
    step = ( 1./(d_min**2) - start ) / float(nbins)
    rshells = " ".join(["%.2f" % (start + i * step)**(-1/2) for i in range(1, nbins+1)])
    return "RESOLUTION_SHELLS= %s\n" % rshells
# make_bin_str()

class XscaleCycles(object):
    def __init__(self, workdir, anomalous_flag, d_min, d_max,
                 reject_method, reject_params, xscale_params, res_params,
                 reference_file, space_group, ref_mtz, out, batch_params, nproc=1):
        self.reference_file = None
        self._counter = 0
        self.workdir_org = workdir
        self.anomalous_flag = anomalous_flag
        self.d_min = d_min # XXX when None..
        self.d_max = d_max
        self.reject_params = reject_params
        self.xscale_params = xscale_params
        self.res_params = res_params
        self.reference_choice = xscale_params.reference if not reference_file else None
        self.space_group = space_group # sgtbx.space_group object. None when user not supplied.
        self.ref_mtz = ref_mtz
        self.out = out
        self.nproc = nproc
        self.nproc_each = batch_params.nproc_each
        if batch_params.engine == "sge": self.batchjobs = batchjob.SGE(pe_name=batch_params.sge_pe_name)
        elif batch_params.engine == "sh": self.batchjobs = batchjob.ExecLocal(max_parallel=batch_params.sh_max_jobs)
        self.all_data_root = None # the root directory for all data
        self.altfile = {} # Modified files
        self.cell_info_at_cycles = {}
        self.dmin_est_at_cycles = {}

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
        self.xscale_inp_head = ""
        if xscale_params.min_i_over_sigma is not None:
            xds_built = check_xds_version()[1]
            snrc_kwd = "MINIMUM_I/SIGMA" if xds_built is not None and xds_built < "20191015" else "SNRC"
            self.xscale_inp_head += "%s= %.2f\n" % (snrc_kwd, xscale_params.min_i_over_sigma)

        self.xscale_inp_head += "OUTPUT_FILE= xscale.hkl\n"

        if self.anomalous_flag is not None:
            self.xscale_inp_head += "FRIEDEL'S_LAW= %s\n" % ("FALSE" if self.anomalous_flag else "TRUE")

        if self.d_min is not None:
            self.xscale_inp_head += make_bin_str(self.d_min, self.d_max) + "\n"

        if reference_file:
            self.reference_file = os.path.join(workdir, "reference.hkl")
            os.symlink(os.path.relpath(reference_file, workdir), self.reference_file)
    # __init__()

    def request_next_workdir(self):
        self._counter += 1
        new_wd = os.path.join(self.workdir_org, "run_%.2d" % self._counter)
        os.mkdir(new_wd)
        print("\nIn %s" % new_wd, file=self.out)

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
                    cell = list(map(float, l[l.index("=")+1:].split()))
                    assert len(cell) == 6
                    cells.append(cell)
                    break

        if self.space_group is not None:
            sg = self.space_group.type().number()

        cells = numpy.array(cells)
        mean_cell = [cells[:,i].mean() for i in range(6)]
        cell_std = [numpy.std(cells[:,i]) for i in range(6)]
        lcv, alcv = blend_lcv.calc_lcv(cells)

        mean_cell_str = " ".join(["%.3f"%x for x in mean_cell])
        cell_std_str = " ".join(["%.1e"%x for x in cell_std])

        print("Averaged cell= %s (%d)" % (mean_cell_str, len(cells)), file=self.out)
        print("Cell Std dev.= %s" % cell_std_str, file=self.out)
        print("LCV, aLCV= %.3f%%, %.3f A" % (lcv, alcv), file=self.out)

        return sg, mean_cell_str, lcv, alcv
    # average_cells()

    def cut_resolution(self, cycle_number):
        def est_resol(xscale_hkl, res_params, plt_out):
            iobs = XDS_ASCII(xscale_hkl, i_only=True).i_obs()
            est = estimate_resolution_based_on_cc_half(iobs, res_params.cc_one_half_min,
                                                       res_params.cc_half_tol,
                                                       res_params.n_bins, log_out=self.out)
            est.show_plot(False, plt_out)
            if None not in (est.d_min, est.cc_at_d_min):
                self.out.write("Best resolution cutoff= %.2f A @CC1/2= %.4f\n" % (est.d_min, est.cc_at_d_min))
            else:
                self.out.write("Can't decide resolution cutoff. No reflections??\n")
            return est.d_min
        # est_resol()

        print("**** Determining resolution cutoff in run_%.2d ****" % cycle_number, file=self.out)
        last_wd = os.path.join(self.workdir_org, "run_%.2d"%cycle_number)
        xscale_hkl = os.path.abspath(os.path.join(last_wd, "xscale.hkl"))

        tmpwd = os.path.join(self.workdir_org, "run_%.2d_tmp"%cycle_number)
        os.mkdir(tmpwd)

        for i, cc_cut in enumerate((self.res_params.cc_one_half_min*.7, self.res_params.cc_one_half_min)):
            self.res_params.cc_one_half_min = cc_cut
            d_min = est_resol(xscale_hkl, self.res_params,
                              os.path.join(tmpwd, "ccfit_%d.pdf"%(i+1)))
            if d_min is not None and d_min > self.d_min + 0.001:
                for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(tmpwd, f))
                inp_new = os.path.join(tmpwd, "XSCALE.INP")
                shutil.copyfile(os.path.join(last_wd, "XSCALE.INP"), inp_new)
                modify_xdsinp(inp_new, [make_bin_str(d_min, self.d_max).split("= ")])

                try:
                    xscale.run_xscale(inp_new, cbf_to_dat=True, aniso_analysis=True,
                                      use_tmpdir_if_available=self.xscale_params.use_tmpdir_if_available)
                except:
                    print(traceback.format_exc(), file=self.out)

                xscale_hkl = os.path.abspath(os.path.join(tmpwd, "xscale.hkl"))

        if not os.path.isfile(os.path.join(tmpwd, "XSCALE.INP")):
            for f in "XSCALE.INP", "XSCALE.LP", "xscale.hkl", "pointless.log", "ccp4":
                os.symlink(os.path.relpath(os.path.join(last_wd, f), tmpwd),
                           os.path.join(tmpwd, f))

        if d_min is not None:
            self.dmin_est_at_cycles[cycle_number] = d_min
            os.rename(tmpwd, os.path.join(self.workdir_org, "run_%.2d_%.2fA"%(cycle_number, d_min)))
    # cut_resolution()

    def estimate_resolution(self, cycle_number):
        print("**** Determining resolution cutoff in run_%.2d ****" % cycle_number, file=self.out)
        last_wd = os.path.join(self.workdir_org, "run_%.2d"%cycle_number)
        xscale_hkl = os.path.abspath(os.path.join(last_wd, "xscale.hkl"))

        i_obs = XDS_ASCII(xscale_hkl, i_only=True).i_obs()
        d_min_est, _ = initial_estimate_byfit_cchalf(i_obs, cc_half_min=self.res_params.cc_one_half_min,
                                                 anomalous_flag=False, log_out=self.out)
                
        self.out.write("Estimated resolution cutoff= %.2f A @CC1/2= %.4f\n" % (d_min_est, self.res_params.cc_one_half_min))
        
        self.dmin_est_at_cycles[cycle_number] = d_min_est
    # estimate_resolution()

    def run_cycles(self, xds_ascii_files):
        self.all_data_root = os.path.dirname(os.path.commonprefix(xds_ascii_files))
        self.removed_files = []
        self.removed_reason = {}
        print("********************* START FUNCTION ***********************", file=self.out)
        if self.reference_file:
            self.run_cycle([self.reference_file,]+xds_ascii_files)
        else:
            self.run_cycle(xds_ascii_files)

        if self.res_params.estimate:
            #self.cut_resolution(self.get_last_cycle_number())
            for run_i in range(1, self.get_last_cycle_number()+1):
                try: self.estimate_resolution(run_i)
                except: print(traceback.format_exc(), file=self.out) # Don't want to stop the program.

        for wd in glob.glob(os.path.join(self.workdir_org, "run_*")):
            if os.path.exists(os.path.join(wd, "ccp4")): continue
            xscale_hkl = os.path.abspath(os.path.join(wd, "xscale.hkl"))
            sg = None # Use user-specified one. Otherwise follow pointless.
            try:
                sg = XDS_ASCII(xscale_hkl, read_data=False).symm.space_group()
                laue_symm_str = str(sg.build_derived_reflection_intensity_group(False).info())
                worker = Pointless()
                result = worker.run_for_symm(xdsin=xscale_hkl,
                                             logout=os.path.join(wd, "pointless.log"),
                                             choose_laue=laue_symm_str,
                                             xdsin_to_p1=True)
                
                if "symm" in result:
                    print("Pointless suggestion (forcing %s symmetry):" % laue_symm_str, file=self.out)
                    result["symm"].show_summary(self.out, " ")
                    sg = str(result["symm"].space_group_info())
                else:
                    print("Pointless failed.", file=self.out)
            except:
                # Don't want to stop the program.
                print(traceback.format_exc(), file=self.out)

            if self.space_group is not None:
                sg = str(self.space_group.info())

            try:
                xds2mtz.xds2mtz(xds_file=xscale_hkl,
                                dir_name=os.path.join(wd, "ccp4"),
                                run_xtriage=True, run_ctruncate=True,
                                with_multiplicity=True,
                                space_group=sg,
                                flag_source=self.ref_mtz)
            except:
                # Don't want to stop the program.
                print(traceback.format_exc(), file=self.out)

        return self.removed_files, self.removed_reason
    # run_cycles()

    def check_remove_list(self, remove_idxes):
        new_list = []
        skip_num = 0
        if self.reference_file: skip_num += 1 # because first is reference file
        #if self._counter > 1: skip_num += 1 # because first (or second) is the previously merged file
        remove_idxes = [x for x in remove_idxes if x >= skip_num]
        return remove_idxes
    # check_remove_list()
    
    def current_working_dir(self): return self.workdir

    def run_cycle(self, xds_ascii_files, reference_idx=None):
        if len(xds_ascii_files) == 0:
            print("Error: no files given.", file=self.out)
            return

        xscale_inp = os.path.join(self.workdir, "XSCALE.INP")
        xscale_lp = os.path.join(self.workdir, "XSCALE.LP")

        # Get averaged cell for scaling
        sg, cell, lcv, alcv = self.average_cells(xds_ascii_files)
        self.cell_info_at_cycles[self.get_last_cycle_number()] = (cell, lcv, alcv)
        
        # Choose directory containing XDS_ASCII.HKL and set space group (but how??)
        inp_out = open(xscale_inp, "w")
        inp_out.write("! This XSCALE.INP is generated by kamo.multi_merge.\n")
        inp_out.write("! You may want to use yamtbx.run_xscale to re-run xscale by yourself\n")
        inp_out.write("! because number of characters in line may exceed the limit of xscale.\n")
        inp_out.write("MAXIMUM_NUMBER_OF_PROCESSORS= %d\n" % self.nproc)
        inp_out.write("SPACE_GROUP_NUMBER= %s\nUNIT_CELL_CONSTANTS= %s\n\n" % (sg, cell))
        inp_out.write(self.xscale_inp_head)

        for i, xds_ascii in enumerate(xds_ascii_files):
            f = self.altfile.get(xds_ascii, xds_ascii)
            tmp = min(os.path.relpath(f, self.workdir), f, key=lambda x:len(x))
            refstr = "*" if i==reference_idx else " "
            inp_out.write(" INPUT_FILE=%s%s\n" % (refstr,tmp))
            if self.d_max is not None:
                d_range = (float("inf") if self.d_max is None else self.d_max,
                           0.           if self.d_min is None else self.d_min)
                inp_out.write("  INCLUDE_RESOLUTION_RANGE= %.4f %.4f\n" % d_range)
            if len(self.xscale_params.corrections) != 3:
                inp_out.write("  CORRECTIONS= %s\n" % " ".join(self.xscale_params.corrections))
            if (self.xscale_params.frames_per_batch, self.xscale_params.degrees_per_batch).count(None) < 2:
                xactmp = XDS_ASCII(f, read_data=False)
                frame_range = xactmp.get_frame_range()
                osc_range = xactmp.osc_range
                nframes = frame_range[1] - frame_range[0] + 1
                if self.xscale_params.frames_per_batch is not None:
                    nbatch = int(numpy.ceil(nframes / self.xscale_params.frames_per_batch))
                else:
                    nbatch = int(numpy.ceil(nframes / self.xscale_params.degrees_per_batch * osc_range))
                print("frame range of %s is %d,%d setting NBATCH= %d" % (f, frame_range[0], frame_range[1], nbatch), file=self.out)
                inp_out.write("  NBATCH= %d\n" % nbatch)

        inp_out.close()

        print("DEBUG:: running xscale with %3d files.." % len(xds_ascii_files), file=self.out)
        try:
            xscale.run_xscale(xscale_inp, cbf_to_dat=True, aniso_analysis=True,
                              use_tmpdir_if_available=self.xscale_params.use_tmpdir_if_available)
        except:
            print(traceback.format_exc(), file=self.out)

        xscale_log = open(xscale_lp).read()
        if "!!! ERROR !!! INSUFFICIENT NUMBER OF COMMON STRONG REFLECTIONS." in xscale_log:
            print("DEBUG:: Need to choose files.", file=self.out)

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
                    max_clique = [x for x in cliques if idx_prevfile in x][-1] # xscale.hkl must be included!

                if self.reference_file:
                    max_clique = [0,] + [x for x in max_clique if x!=0]

                for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(self.workdir, f))

                try_later = [xds_ascii_files[i] for i in [x for x in G.nodes() if x not in max_clique]]

                print("DEBUG:: %d files can be merged. %d files will be merged later." % (len(max_clique),
                                                                                                      len(try_later)), file=self.out)
                print("DEBUG:: %d files are of no use." % (len(xds_ascii_files)-len(G.nodes())), file=self.out)
                for i in [j for j in range(len(xds_ascii_files)) if j not in G.nodes()]:
                    self.removed_files.append(xds_ascii_files[i])
                    self.removed_reason[xds_ascii_files[i]] = "no_common_refls"

                self.run_cycle([xds_ascii_files[i] for i in max_clique])

                assert len(try_later) <= 0 # Never be the case with newer xscale!! (if the case, check_remove_list() should be modified to skip_num+=1
                if len(try_later) > 0:
                    print("Trying to merge %d remaining files.." % len(try_later), file=self.out)
                    next_files = [os.path.join(self.workdir, "xscale.hkl")] + try_later
                    if self.reference_file: next_files = [self.reference_file,] + next_files
                    self.workdir = self.request_next_workdir()
                    self.run_cycle(next_files)
                    return
            else:
                bad_idxes = xscalelp.read_no_common_ref_datasets(xscale_lp)
                print("DEBUG:: %d files are of no use." % (len(bad_idxes)), file=self.out)

                for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(self.workdir, f))

                # XXX Actually, not all datasets need to be thrown.. some of them are useful..
                for i in bad_idxes:
                    self.removed_files.append(xds_ascii_files[i])
                    self.removed_reason[xds_ascii_files[i]] = "no_common_refls"

                self.run_cycle([xds_ascii_files[i] for i in [j for j in range(len(xds_ascii_files)) if j not in bad_idxes]])

            return
        elif "!!! ERROR !!! USELESS DATA ON INPUT REFLECTION FILE" in xscale_log:
            print("DEBUG:: Need to discard useless data.", file=self.out)
            unuseful_data = [xscalelp.get_read_data(xscale_lp)[-1]] #filter(lambda x: x[2]==0, xscalelp.get_read_data(xscale_lp))
            if len(unuseful_data) == 0:
                print("I don't know how to fix it..", file=self.out)
                return
            remove_idxes = [x[0]-1 for x in unuseful_data]
            remove_idxes = self.check_remove_list(remove_idxes)
            keep_idxes = [x for x in range(len(xds_ascii_files)) if x not in remove_idxes]
            for i in remove_idxes:
                self.removed_files.append(xds_ascii_files[i])
                self.removed_reason[xds_ascii_files[i]] = "useless"

            for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(self.workdir, f))
            self.run_cycle([xds_ascii_files[i] for i in keep_idxes])
            return
        elif "INACCURATE SCALING FACTORS." in xscale_log:
            # Actually I don't know how to fix this.. (bug?) but worth proceeding (discarding bad data may solve problem).
            print("'INACCURATE SCALING FACTORS' happened.. but ignored.", file=self.out)
        elif "!!! ERROR !!!" in xscale_log:
            print("Unknown error! please check the XSCALE.LP and fix the program.", file=self.out)
            return

        # Re-scale by changing reference
        rescale_for = None
        if len(self.reject_method) == 0:
            rescale_for = self.reference_choice # may be None
        elif reference_idx is None:
            rescale_for = "bmed"
        
        if rescale_for is not None and len(xds_ascii_files) > 1:
            ref_num = xscale.decide_scaling_reference_based_on_bfactor(xscale_lp, rescale_for, return_as="index")
            if reference_idx is None:
                print("Rescaling with %s" % rescale_for, file=self.out)
                for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(self.workdir, f))
                self.run_cycle(xds_ascii_files, reference_idx=ref_num)

        if len(self.reject_method) == 0:
            return

        # Remove bad data
        remove_idxes = []
        remove_reasons = {}

        if self.reject_method[0] == "framecc":
            print("Rejections based on frame CC", file=self.out)
            from yamtbx.dataproc.xds.command_line import xscale_cc_against_merged

            # list of [frame, n_all, n_common, cc] in the same order
            framecc = list(xscale_cc_against_merged.run(hklin=os.path.join(self.workdir, "xscale.hkl"),
                                                   output_dir=self.workdir).values())
            if self.reject_params.framecc.method == "tukey":
                ccs = numpy.array([x[3] for x in reduce(lambda x,y:x+y,framecc)])
                ccs = ccs[ccs==ccs] # Remove nan
                q25, q75 = numpy.percentile(ccs, [25, 75])
                cc_cutoff  = q25 - self.reject_params.framecc.iqr_coeff * (q75 - q25)
                print(" frameCC cutoff = %.4f (%.2f*IQR)" % (cc_cutoff, self.reject_params.framecc.iqr_coeff), file=self.out)
            else:
                cc_cutoff = self.reject_params.framecc.abs_cutoff
                print(" frameCC cutoff = %.4f (value specified)" % cc_cutoff, file=self.out)

            for i, cclist in enumerate(framecc):
                useframes = [x[0] for x in [x for x in cclist if x[3] > cc_cutoff]]
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

                print("Extracting frames %s out of %d-%d in %s" % (",".join(map(str,useframes)),
                                                                               min(xac.iframe), max(xac.iframe),
                                                                               f), file=self.out)

                newf = self.request_file_modify(f)
                xac.write_selected(sel, newf)

            self.reject_method.pop(0) # Perform only once

        elif self.reject_method[0] == "lpstats":
            if "bfactor" in self.reject_params.lpstats.stats:
                iqrc = self.reject_params.lpstats.iqr_coeff
                print("Rejections based on B-factor outliers (%.2f*IQR)" % iqrc, file=self.out)
                Bs = numpy.array([x[1] for x in xscalelp.get_k_b(xscale_lp)])
                if len(Bs)>1: # If one data, K & B table is not available.
                    q25, q75 = numpy.percentile(Bs, [25, 75])
                    iqr = q75 - q25
                    lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
                    count = 0
                    for i, b in enumerate(Bs):
                        if b < lowlim or b > highlim:
                            remove_idxes.append(i)
                            remove_reasons.setdefault(i, []).append("bad_B")
                            count += 1

                    print(" %4d B-factor outliers (<%.2f, >%.2f) removed"% (count, lowlim, highlim), file=self.out)
                else:
                    print(" B-factor outlier rejection is not available.", file=self.out)

            if "em.b" in self.reject_params.lpstats.stats:
                iqrc = self.reject_params.lpstats.iqr_coeff
                print("Rejections based on error model b outliers (%.2f*IQR)" % iqrc, file=self.out)
                bs = numpy.array([x[1] for x in xscalelp.get_ISa(xscale_lp)])
                q25, q75 = numpy.percentile(bs, [25, 75])
                iqr = q75 - q25
                lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
                count = 0
                for i, b in enumerate(bs):
                    if b < lowlim or b > highlim:
                        remove_idxes.append(i)
                        remove_reasons.setdefault(i, []).append("bad_em.b")
                        count += 1

                print(" %4d error model b outliers (<%.2f, >%.2f) removed"% (count, lowlim, highlim), file=self.out)

            if "em.ab" in self.reject_params.lpstats.stats:
                iqrc = self.reject_params.lpstats.iqr_coeff
                print("Rejections based on error model a*b outliers (%.2f*IQR)" % iqrc, file=self.out)
                vals = numpy.array([x[0]*x[1] for x in xscalelp.get_ISa(xscale_lp)])
                q25, q75 = numpy.percentile(vals, [25, 75])
                iqr = q75 - q25
                lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
                count = 0
                for i, ab in enumerate(vals):
                    if ab < lowlim or ab > highlim:
                        remove_idxes.append(i)
                        remove_reasons.setdefault(i, []).append("bad_em.ab")
                        count += 1

                print(" %4d error model a*b outliers (<%.2f, >%.2f) removed"% (count, lowlim, highlim), file=self.out)

            if "rfactor" in self.reject_params.lpstats.stats:
                iqrc = self.reject_params.lpstats.iqr_coeff
                print("Rejections based on R-factor outliers (%.2f*IQR)" % iqrc, file=self.out)
                rstats = xscalelp.get_rfactors_for_each(xscale_lp)
                vals = numpy.array([rstats[x][-1][1] for x in rstats]) # Read total R-factor
                q25, q75 = numpy.percentile(vals, [25, 75])
                iqr = q75 - q25
                lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
                count = 0
                for i, v in enumerate(vals):
                    if v < lowlim or v > highlim:
                        remove_idxes.append(i)
                        remove_reasons.setdefault(i, []).append("bad_R")
                        count += 1

                print(" %4d R-factor outliers (<%.2f, >%.2f) removed"% (count, lowlim, highlim), file=self.out)

            if "pairwise_cc" in self.reject_params.lpstats.stats:
                corrs = xscalelp.get_pairwise_correlations(xscale_lp)
                if self.reject_params.lpstats.pwcc.method == "tukey":
                    q25, q75 = numpy.percentile([x[3] for x in corrs], [25, 75])
                    iqr = q75 - q25
                    lowlim = q25 - self.reject_params.lpstats.pwcc.iqr_coeff * iqr
                    print("Rejections based on pairwise_cc < %.4f (IQR=%.2f)" % (lowlim, iqr), file=self.out)
                else:
                    lowlim = self.reject_params.lpstats.pwcc.abs_cutoff
                    print("Rejections based on pairwise_cc < %.4f" % lowlim, file=self.out)

                bad_corrs = [x for x in corrs if x[3] < lowlim]
                idx_bad = {}
                for i, j, common_refs, corr, ratio, bfac in bad_corrs:
                    idx_bad[i] = idx_bad.get(i, 0) + 1
                    idx_bad[j] = idx_bad.get(j, 0) + 1

                idx_bad = list(idx_bad.items())
                idx_bad.sort(key=lambda x:x[1])
                count = 0
                for idx, badcount in reversed(idx_bad):
                    remove_idxes.append(idx-1)
                    remove_reasons.setdefault(idx-1, []).append("bad_pwcc")
                    bad_corrs = [x for x in bad_corrs if idx not in x[:2]]
                    if len(bad_corrs) == 0: break
                    fun_key = lambda x: x[3]
                    print(" Removing idx=%d (CC %.3f..%.3f) remaining %d bad pairs" % (idx, 
                                                                                                   min(bad_corrs,key=fun_key)[3],
                                                                                                   max(bad_corrs,key=fun_key)[3],
                                                                                                   len(bad_corrs)), file=self.out)
                    count += 1
                print(" %4d pairwise CC outliers removed" % count, file=self.out)

            self.reject_method.pop(0) # Perform only once
        elif self.reject_method[0] == "delta_cc1/2":
            print("Rejection based on delta_CC1/2 in %s shell" % self.delta_cchalf_bin, file=self.out)
            table = xscalelp.read_stats_table(xscale_lp)
            i_stat = -1 if self.delta_cchalf_bin == "total" else -2
            prev_cchalf = table["cc_half"][i_stat]
            prev_nuniq = table["nuniq"][i_stat]
            # file_name->idx table
            remaining_files = collections.OrderedDict([x[::-1] for x in enumerate(xds_ascii_files)])

            # For consistent resolution limit
            inp_head = self.xscale_inp_head + "SPACE_GROUP_NUMBER= %s\nUNIT_CELL_CONSTANTS= %s\n\n" % (sg, cell)
            count = 0
            for i in range(len(xds_ascii_files)-1): # if only one file, cannot proceed.
                tmpdir = os.path.join(self.workdir, "reject_test_%.3d" % i)

                cchalf_list = xscale.calc_cchalf_by_removing(wdir=tmpdir, inp_head=inp_head,
                                                             inpfiles=list(remaining_files.keys()),
                                                             stat_bin=self.delta_cchalf_bin,
                                                             nproc=self.nproc,
                                                             nproc_each=self.nproc_each,
                                                             batchjobs=self.batchjobs)

                rem_idx, cc_i, nuniq_i = cchalf_list[0] # First (largest) is worst one to remove.
                rem_idx_in_org = remaining_files[list(remaining_files.keys())[rem_idx]]
                
                # Decision making by CC1/2
                print("DEBUG:: cycle %.3d remove %3d if %.2f*%d > %.2f*%d" % (i, rem_idx_in_org, 
                                                                                          cc_i, nuniq_i,
                                                                                          prev_cchalf, prev_nuniq), file=self.out)
                if cc_i*nuniq_i <= prev_cchalf*prev_nuniq: break
                print("Removing idx= %3d gained CC1/2 by %.2f" % (rem_idx_in_org, cc_i-prev_cchalf), file=self.out)

                prev_cchalf, prev_nuniq = cc_i, nuniq_i
                remove_idxes.append(rem_idx_in_org)
                remove_reasons.setdefault(rem_idx_in_org, []).append("bad_cchalf")
                del remaining_files[list(remaining_files.keys())[rem_idx]] # remove file from table
                count += 1

            print(" %4d removed by DeltaCC1/2 method" % count, file=self.out)

            if self.next_delta_cchalf_bin != []:
                self.delta_cchalf_bin = self.next_delta_cchalf_bin.pop(0)
            else:
                self.reject_method.pop(0)
        else:
            print("ERROR:: Unsupported reject_method (%s)" % reject_method, file=self.out)

        # Remove duplicates
        remove_idxes = list(set(remove_idxes))
        remove_idxes = self.check_remove_list(remove_idxes)
        if len(remove_idxes) > 0:
            print("DEBUG:: Need to remove %d files" % len(remove_idxes), file=self.out)
            for i in sorted(remove_idxes): 
                print(" %.3d %s" % (i, xds_ascii_files[i]), file=self.out)
                self.removed_files.append(xds_ascii_files[i])
                self.removed_reason[xds_ascii_files[i]] = ",".join(remove_reasons[i])

        # Next run
        keep_idxes = [x for x in range(len(xds_ascii_files)) if x not in remove_idxes]
        if len(self.reject_method) > 0 or len(remove_idxes) > 0:
            self.workdir = self.request_next_workdir()
            self.run_cycle([xds_ascii_files[i] for i in keep_idxes])
        elif self.reference_choice is not None and len(keep_idxes) > 1:
            # Just re-scale with B reference
            ref_num = xscale.decide_scaling_reference_based_on_bfactor(xscale_lp, self.reference_choice, return_as="index")
            if reference_idx is None: # this would not happen? why do we need this?
                print("Rescaling2 with %s" % self.reference_choice, file=self.out)
                for f in "XSCALE.INP", "XSCALE.LP": util.rotate_file(os.path.join(self.workdir, f))
                self.run_cycle([xds_ascii_files[i] for i in keep_idxes], reference_idx=ref_num)

    # run_cycle()
# class XscaleCycles
