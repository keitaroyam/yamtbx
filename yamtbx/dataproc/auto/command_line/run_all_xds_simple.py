"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import sys
import os
import shutil
import traceback
import subprocess
import re
import pickle
import time
import glob
import numpy

from yamtbx.dataproc.xds import get_xdsinp_keyword, modify_xdsinp, optimal_delphi_by_nproc, make_backup, revert_files, remove_backups
from yamtbx.dataproc.xds import idxreflp
from yamtbx.dataproc.xds import correctlp
from yamtbx.dataproc.xds.command_line import xds_plot_integrate
from yamtbx.dataproc.xds import files as xds_files
from yamtbx.dataproc.xds.xparm import XPARM
from yamtbx.dataproc.auto import resolution_cutoff
from yamtbx.dataproc.auto import html_report
from yamtbx.dataproc.pointless import Pointless
from yamtbx import util
from yamtbx.util import xtal

import iotbx.phil
from libtbx.utils import Sorry
from libtbx import easy_mp
from libtbx.utils import multi_out
from cctbx import crystal
from cctbx.crystal import reindex

master_params_str = """
topdir = None
 .type = path
dont_overwrite = False
 .type = bool
fast_delphi = True
 .type = bool
use_pointless = True
 .type = bool
show_progress = True
 .type = bool
nproc = None
 .type = int
 .help = number of processors for single xds job OR number of parallel jobs
multiproc = False
 .type = bool
 .help = Parallel processing of multiple xds jobs
parmethod = *multiprocessing sge
 .type = choice(multi=False)
 .help = multiprocessing method
resume = False
 .type = bool
 .help = Skip processed dir
tryhard = False
 .type = bool
 .help = Do my best.
mode = *initial recycle
 .type = choice(multi=False)
 .help = initial= initial runs. recyle=re-integrate with post-refined parameter
cut_resolution = True
 .type = bool
 .help = automatically cut resolution
no_scaling = False
 .type = bool
 .help = Run CORRECT, but no scaling is applied (for merging small wedge)
make_report = True
 .type = bool
 .help = Create html report
use_tmpdir_if_available = true
 .type = bool
 .help = Use ramdisk or tempdir if sufficient size is available
auto_frame_exclude_spot_based = false
 .type = bool
 .help = automatic frame exclusion from integration based on spot search result.
cell_prior {
 method = *not_use_first use_first symm_constraint_only correct_only
  .type = choice(multi=False)
  .help = "not_use_first: Try indexing without prior information first, and if failed, use prior."
          "use_first: Try indexing with prior information."
          "symm_constraint_only: Try indexing without prior information, and apply symmetry constraints for determined unit cell."
          "correct_only: Use given symmetry in CORRECT. May be useful in recycling."
 check = true
  .type = bool
 cell = None
  .type = floats(size=6)
 sgnum = 0
  .type = int
 tol_length = 0.05
  .type = float
  .help = relative_length_tolerance
 tol_angle = 5
  .type = float
  .help = absolute_angle_tolerance in degree
}
"""

re_running_job = re.compile("\*\*\*\*\* ([^ ]*) \*\*\*\*\*")
re_running_integrate = re.compile("PROCESSING OF IMAGES *([0-9]*) *\.\.\. *([0-9]*)")

def find_mosaicity_for_image(line):
    if line[6:10] == "   0":
        sp = line.split()
        if len(sp)== 10:
            return float(sp[-1])
# find_mosaicity_for_image()

def calc_merging_stats(xac_file, cut_resolution=True):
    import iotbx.merging_statistics
    from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII

    wdir = os.path.dirname(xac_file)
    pklout = os.path.join(wdir, "merging_stats.pkl")
    logout = open(os.path.join(wdir, "merging_stats.log"), "w")

    print >>logout, xac_file
    print >>logout, ""
    print >>logout, "Estimate cutoff"
    print >>logout, "================"

    obj = XDS_ASCII(xac_file, i_only=True)
    i_obs = obj.i_obs()
    d_min = None
    if i_obs.size() < 10: return

    try:
        cutoffs = resolution_cutoff.estimate_crude_resolution_cutoffs(i_obs=i_obs)
        cutoffs.show(out=logout)

        if cutoffs.cc_one_half_cut != float("inf") and cut_resolution:
            d_min = cutoffs.cc_one_half_cut
    except Sorry, e:
        print >>logout, e.message

    print >>logout, ""
    print >>logout, "Merging statistics"
    print >>logout, "==================="

    try:
        stats = iotbx.merging_statistics.dataset_statistics(i_obs=i_obs,
                                                            crystal_symmetry=obj.symm,
                                                            d_min=d_min,
                                                            d_max=None,
                                                            n_bins=10,
                                                            anomalous=obj.anomalous,
                                                            sigma_filtering="xds",
                                                            log=logout)
        stats.show(out=logout)
    except (Sorry, RuntimeError) as e:
        print >>logout, e.message
        return

    ret = dict(cutoff=d_min, cutoffs=cutoffs, stats=stats)
    pickle.dump(ret, open(pklout, "w"))
    return d_min, cutoffs, stats
# calc_merging_stats()

def run_xds(wdir, comm="xds_par", show_progress=True):
    env = None

    if "SGE_O_PATH" in os.environ:
        env = os.environ.copy()
        env["PATH"] = env["SGE_O_PATH"] + ":" + env["PATH"]

    if show_progress:
        p = subprocess.Popen(comm, cwd=wdir, stdout=subprocess.PIPE, env=env)
        read_mosaicity_flag = False
        mosaicity_sofar = []
        for l in iter(p.stdout.readline,''):
            #print l
            r = re_running_job.search(l)
            if r:
                sys.stdout.write("\r\x1b[K Running %s " % r.group(1))
                if r.group(1) == "INTEGRATE":
                    read_mosaicity_flag = True

            if read_mosaicity_flag:
                m = find_mosaicity_for_image(l)
                if m: mosaicity_sofar.append(m)

            r = re_running_integrate.search(l)
            if r:
                sys.stdout.write("\r\x1b[K Running INTEGRATE  %5d ...%5d" % tuple(map(int, r.groups())))
                if len(mosaicity_sofar) > 0:
                    mean_mosaicity = sum(mosaicity_sofar)/len(mosaicity_sofar)
                    sys.stdout.write(" (mosaicity so far: mean=%.2f)" % mean_mosaicity)
    else:
        p = subprocess.Popen(comm, cwd=wdir, stdout=subprocess.PIPE, env=env)
        p.wait()
# run_xds()

def run_xdsstat(wdir):
    comm = "xdsstat"
    out = open(os.path.join(wdir, "XDSSTAT.LP"), "w")
    sys.stdout.write("\r\x1b[K Running XDSSTAT")
    util.call(comm, wdir=wdir,
              stdin="\n", stdout=out)
# run_xdsstat()

def try_indexing_hard(wdir, show_progress, decilog, 
                      known_sgnum=None, known_cell=None, tol_length=None, tol_angle=None):
    idxref_lp = os.path.join(wdir, "IDXREF.LP")
    xdsinp = os.path.join(wdir, "XDS.INP")

    lp_org = idxreflp.IdxrefLp(idxref_lp)

    if lp_org.is_cell_maybe_half():
        backup_needed = ("XDS.INP",) + xds_files.generated_by_IDXREF

        print >>decilog, " !! Cell may be halved. Trying doubled cell."
        bk_prefix = make_backup(backup_needed, wdir=wdir, quiet=True)

        cell = lp_org.deduce_correct_cell_based_on_integerness()
        cell = " ".join(map(lambda x:"%.2f"%x, cell.parameters()))
        modify_xdsinp(xdsinp, inp_params=[("JOB", "IDXREF"),
                                          ("SPACE_GROUP_NUMBER", "1"),
                                          ("UNIT_CELL_CONSTANTS", cell)
                                          ])

        for f in xds_files.generated_by_IDXREF: util.rotate_file(os.path.join(wdir, f), copy=(f=="SPOT.XDS"))
        run_xds(wdir=wdir, show_progress=show_progress)

        if idxreflp.IdxrefLp(idxref_lp).is_cell_maybe_half():
            revert_files(backup_needed, bk_prefix, wdir=wdir, quiet=True)

            print >>decilog, "  .. not solved. Next, try decreasing SEPMIN= and CLUSTER_RADIUS=."
            bk_prefix = make_backup(backup_needed, wdir=wdir, quiet=True)

            modify_xdsinp(xdsinp, inp_params=[("JOB", "IDXREF"),
                                              ("SEPMIN", "4"),
                                              ("CLUSTER_RADIUS", "2")
                                              ])
            for f in xds_files.generated_by_IDXREF: util.rotate_file(os.path.join(wdir, f), copy=(f=="SPOT.XDS"))
            run_xds(wdir=wdir, show_progress=show_progress)

            if idxreflp.IdxrefLp(idxref_lp).is_cell_maybe_half():
                print >>decilog, "  .. not solved. Give up."
                revert_files(backup_needed, bk_prefix, wdir=wdir, quiet=True)
        else:
            print >>decilog, "  Now OK."
            remove_backups(backup_needed, bk_prefix, wdir=wdir)
            modify_xdsinp(xdsinp, inp_params=[("SPACE_GROUP_NUMBER", "0"),
                                              ])

    # If Cell hint exists, try to use it..
    if known_sgnum > 0:
        flag_try_cell_hint = False
        xparm = os.path.join(wdir, "XPARM.XDS")
        if not os.path.isfile(xparm):
            flag_try_cell_hint = True
        else:
            xsxds = XPARM(xparm).crystal_symmetry()

            xsref = crystal.symmetry(known_cell, known_sgnum)
            cosets = reindex.reindexing_operators(xsref, xsxds,
                                                  tol_length, tol_angle)

            if cosets.double_cosets is None: flag_try_cell_hint = True

        if flag_try_cell_hint:
            print >>decilog, " Worth trying to use prior cell for indexing."
            modify_xdsinp(xdsinp, inp_params=[("JOB", "IDXREF"),
                                              ("UNIT_CELL_CONSTANTS",
                                               " ".join(map(lambda x: "%.3f"%x, known_cell))),
                                              ("SPACE_GROUP_NUMBER", "%d"%known_sgnum),
                                              ])
            for f in xds_files.generated_by_IDXREF: util.rotate_file(os.path.join(wdir, f),
                                                                     copy=(f=="SPOT.XDS"))
            run_xds(wdir=wdir, show_progress=False)
            modify_xdsinp(xdsinp, inp_params=[("SPACE_GROUP_NUMBER", "0"),
                                              ])           
# try_indexing_hard()

def xds_sequence(root, params):
    print
    print os.path.relpath(root, params.topdir)

    xparm = os.path.join(root, "XPARM.XDS")
    gxparm = os.path.join(root, "GXPARM.XDS")
    defpix_lp = os.path.join(root, "DEFPIX.LP")
    correct_lp = os.path.join(root, "CORRECT.LP")
    integrate_hkl = os.path.join(root, "INTEGRATE.HKL")
    xac_hkl = os.path.join(root, "XDS_ASCII.HKL")
    integrate_lp = os.path.join(root, "INTEGRATE.LP")
    spot_xds = os.path.join(root, "SPOT.XDS")
    xdsinp = os.path.join(root, "XDS.INP")

    assert os.path.isfile(xdsinp)

    xdsinp_dict = dict(get_xdsinp_keyword(xdsinp))

    decilog = multi_out()
    decilog.register("log", open(os.path.join(root, "decision.log"), "a"), atexit_send_to=None)

    print >>decilog, "xds_sequence started at %s in %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), root)
    
    if params.show_progress:
        decilog.register("stdout", sys.stdout)

    if params.mode=="initial" and params.resume and os.path.isfile(correct_lp):
        print " Already processed."
        return

    if params.mode == "recycle" and not os.path.isfile(gxparm):
        print "GXPARM.XDS not found. Cannot do recycle."
        return

    if params.fast_delphi and (params.nproc is None or params.nproc > 1):
        delphi = optimal_delphi_by_nproc(xdsinp=xdsinp, nproc=params.nproc)
        print " Setting delphi to ", delphi
        modify_xdsinp(xdsinp, inp_params=[("DELPHI", str(delphi)),
                                          ])

    if params.nproc is not None and params.nproc > 1:
        modify_xdsinp(xdsinp, inp_params=[("MAXIMUM_NUMBER_OF_PROCESSORS", str(params.nproc)),
                                          ])

    if params.mode == "initial":
        # Peak search
        modify_xdsinp(xdsinp, inp_params=[("JOB", "XYCORR INIT COLSPOT")])
        run_xds(wdir=root, show_progress=params.show_progress)
        if params.auto_frame_exclude_spot_based:
            sx = idxreflp.SpotXds(spot_xds)
            sx.set_xdsinp(xdsinp)
            spots = filter(lambda x: 5 < x[-1] < 30, sx.collected_spots()) # low-res (5 A)
            frame_numbers = numpy.array(map(lambda x: int(x[2])+1, spots))
            data_range = map(int, xdsinp_dict["DATA_RANGE"].split())
            # XXX this assumes SPOT_RANGE equals to DATA_RANGE. Is this guaranteed?
            h = numpy.histogram(frame_numbers,
                                bins=numpy.arange(data_range[0], data_range[1]+2, step=1))
            q14 = numpy.percentile(h[0], [25,75])
            iqr = q14[1] - q14[0]
            cutoff = max(h[0][h[0]<=iqr*1.5+q14[1]]) / 5 # magic number
            print "DEBUG:: IQR= %.2f, Q1/4= %s, cutoff= %.2f" % (iqr,q14, cutoff)
            cut_frames = h[1][h[0]<cutoff]
            keep_frames = h[1][h[0]>=cutoff]
            print "DEBUG:: keep_frames=", keep_frames
            print "DEBUG::  cut_frames=", cut_frames

            if len(cut_frames) > 0:
                cut_ranges = [[cut_frames[0], cut_frames[0]], ]
                for fn in cut_frames:
                    if fn - cut_ranges[-1][1] <= 1: cut_ranges[-1][1] = fn
                    else: cut_ranges.append([fn, fn])

                # Edit XDS.INP
                cut_inp_str = "".join(map(lambda x: "EXCLUDE_DATA_RANGE= %6d %6d\n"%tuple(x), cut_ranges))
                open(xdsinp, "a").write("\n"+cut_inp_str)

                # Edit SPOT.XDS
                shutil.copyfile(spot_xds, spot_xds+".org")
                sx.write(open(spot_xds, "w"), frame_selection=set(keep_frames))

        # Indexing
        if params.cell_prior.method == "use_first":
            modify_xdsinp(xdsinp, inp_params=[("JOB", "IDXREF"),
                                              ("UNIT_CELL_CONSTANTS",
                                               " ".join(map(lambda x: "%.3f"%x, params.cell_prior.cell))),
                                              ("SPACE_GROUP_NUMBER", "%d"%params.cell_prior.sgnum),
                                              ])
        else:
            modify_xdsinp(xdsinp, inp_params=[("JOB", "IDXREF")])

        run_xds(wdir=root, show_progress=params.show_progress)
        print # TODO indexing stats like indexed percentage here.

        if params.tryhard:
            try_indexing_hard(root, params.show_progress, decilog,
                              known_sgnum=params.cell_prior.sgnum,
                              known_cell=params.cell_prior.cell,
                              tol_length=params.cell_prior.tol_length,
                              tol_angle=params.cell_prior.tol_angle)

        if not os.path.isfile(xparm):
            print >>decilog, " Indexing failed."
            return

        if params.cell_prior.sgnum > 0:
            # Check anyway
            xsxds = XPARM(xparm).crystal_symmetry()
            xsref = crystal.symmetry(params.cell_prior.cell, params.cell_prior.sgnum)
            cosets = reindex.reindexing_operators(xsref, xsxds,
                                                  params.cell_prior.tol_length, params.cell_prior.tol_angle)
            if cosets.double_cosets is None:
                if params.cell_prior.check:
                    print >>decilog, " Incompatible cell. Indexing failed."
                    return
                else:
                    print >>decilog, " Warning: Incompatible cell."

            elif params.cell_prior.method == "symm_constraint_only":
                cell = xsxds.unit_cell().change_basis(cosets.combined_cb_ops()[0])
                print >>decilog, " Trying symmetry-constrained cell parameter:", cell
                modify_xdsinp(xdsinp, inp_params=[("JOB", "IDXREF"),
                                                  ("UNIT_CELL_CONSTANTS",
                                                   " ".join(map(lambda x: "%.3f"%x, cell.parameters()))),
                                                  ("SPACE_GROUP_NUMBER", "%d"%params.cell_prior.sgnum),
                                                  ])
                for f in xds_files.generated_by_IDXREF:
                    util.rotate_file(os.path.join(root, f), copy=(f=="SPOT.XDS"))
                
                run_xds(wdir=root, show_progress=params.show_progress)
                
                if not os.path.isfile(xparm):
                    print >>decilog, " Indexing failed."
                    return

                # Check again
                xsxds = XPARM(xparm).crystal_symmetry()
                if not xsxds.unit_cell().is_similar_to(xsref.unit_cell(),
                                                       params.cell_prior.tol_length, params.cell_prior.tol_angle):
                    print >>decilog, "  Resulted in different cell. Indexing failed."
                    return

    elif params.mode == "recycle":
        print " Start recycle. original ISa= %.2f" % correctlp.get_ISa(correct_lp, check_valid=True)
        for f in xds_files.generated_after_DEFPIX + ("XPARM.XDS", "plot_integrate.log"):
            util.rotate_file(os.path.join(root, f), copy=True)
        shutil.copyfile(gxparm+".1", xparm)
    else:
        raise "Unknown mode (%s)" % params.mode

    # To Integration
    modify_xdsinp(xdsinp, inp_params=[("JOB", "DEFPIX INTEGRATE"),
                                      ("INCLUDE_RESOLUTION_RANGE", "50 0")])
    run_xds(wdir=root, show_progress=params.show_progress)
    if os.path.isfile(integrate_lp):
        xds_plot_integrate.run(integrate_lp, os.path.join(root, "plot_integrate.log"))
    if not os.path.isfile(integrate_hkl):
        print >>decilog, " Integration failed."
        return


    # Make _noscale.HKL if needed
    if params.no_scaling:
        bk_prefix = make_backup(("XDS.INP",), wdir=root, quiet=True)
        xparm_obj = XPARM(xparm)
        modify_xdsinp(xdsinp, inp_params=[("JOB", "CORRECT"),
                                          ("CORRECTIONS", ""),
                                          ("NBATCH", "1"),
                                          ("MINIMUM_I/SIGMA", "50"),
                                          ("REFINE(CORRECT)", ""),
                                          ("UNIT_CELL_CONSTANTS", " ".join(map(lambda x:"%.3f"%x, xparm_obj.unit_cell))),
                                          ("SPACE_GROUP_NUMBER", "%d"%xparm_obj.spacegroup),])
        print >>decilog, " running CORRECT without empirical scaling"
        run_xds(wdir=root, show_progress=params.show_progress)
        for f in xds_files.generated_by_CORRECT + ("XDS.INP",):
            ff = os.path.join(root, f)
            if not os.path.isfile(ff): continue
            if ff.endswith(".cbf"):
                os.remove(ff)
            else:
                os.rename(ff, ff+"_noscale")

        revert_files(("XDS.INP",), bk_prefix, wdir=root, quiet=True)

    # Run pointless
    symm_by_integrate = None
    if params.use_pointless:
        worker = Pointless()
        result = worker.run_for_symm(xdsin=integrate_hkl, 
                                     logout=os.path.join(root, "pointless_integrate.log"))
        if "symm" in result:
            symm = result["symm"]
            print >>decilog, " pointless using INTEGRATE.HKL suggested", symm.space_group_info()
            sgnum = symm.space_group_info().type().number()
            cell = " ".join(map(lambda x:"%.2f"%x, symm.unit_cell().parameters()))
            modify_xdsinp(xdsinp, inp_params=[("SPACE_GROUP_NUMBER", "%d"%sgnum),
                                              ("UNIT_CELL_CONSTANTS", cell)])
            symm_by_integrate = symm
        else:
            print >>decilog, " pointless failed."

    flag_do_not_change_symm = False
    if params.cell_prior.method == "correct_only":
        xsxds = XPARM(xparm).crystal_symmetry()
        xsref = crystal.symmetry(params.cell_prior.cell, params.cell_prior.sgnum)
        cosets = reindex.reindexing_operators(xsref, xsxds,
                                              params.cell_prior.tol_length, params.cell_prior.tol_angle)
        if cosets.double_cosets is not None:
            cell = xsxds.unit_cell().change_basis(cosets.combined_cb_ops()[0])
            print >>decilog, " Using given symmetry in CORRECT with symmetry constraints:", cell
            modify_xdsinp(xdsinp, inp_params=[("UNIT_CELL_CONSTANTS",
                                               " ".join(map(lambda x: "%.3f"%x, cell.parameters()))),
                                              ("SPACE_GROUP_NUMBER", "%d"%params.cell_prior.sgnum),
            ])
            flag_do_not_change_symm = True
        else:
            print >>decilog, " Tried to use given symmetry in CORRECT, but cell in integration is incompatible."
            
            
    # Do Scaling
    modify_xdsinp(xdsinp, inp_params=[("JOB", "CORRECT"),])

    run_xds(wdir=root, show_progress=params.show_progress)

    if not os.path.isfile(gxparm):
        print >>decilog, " Scaling failed."
        return

    print >>decilog, " OK. ISa= %.2f" % correctlp.get_ISa(correct_lp, check_valid=True)

    ret = calc_merging_stats(os.path.join(root, "XDS_ASCII.HKL"))
    if params.cut_resolution:
        if ret is not None and ret[0] is not None:
            d_min = ret[0]
            modify_xdsinp(xdsinp, inp_params=[("JOB", "CORRECT"),
                                              ("INCLUDE_RESOLUTION_RANGE", "50 %.2f"%d_min)])
            print >>decilog, " Re-scale at %.2f A" % d_min
            os.rename(os.path.join(root, "CORRECT.LP"), os.path.join(root, "CORRECT_fullres.LP"))
            os.rename(os.path.join(root, "XDS_ASCII.HKL"), os.path.join(root, "XDS_ASCII_fullres.HKL"))
            run_xds(wdir=root, show_progress=params.show_progress)
            print >>decilog, " OK. ISa= %.2f" % correctlp.get_ISa(correct_lp, check_valid=True)
            print >>decilog, " (Original files are saved as *_fullres.*)"
        else:
            print >>decilog, "error: Can't decide resolution."

    last_ISa = correctlp.get_ISa(correct_lp, check_valid=True)

    # Run pointless and (if result is different from INTEGRATE) re-scale.
    if params.use_pointless:
        worker = Pointless()
        result = worker.run_for_symm(xdsin=xac_hkl,
                                     logout=os.path.join(root, "pointless_correct.log"))
        if "symm" in result:
            symm = result["symm"]
            need_rescale = False

            if symm_by_integrate is not None:
                if not xtal.is_same_laue_symmetry(symm_by_integrate.space_group(), symm.space_group()):
                    print >>decilog, "pointless suggested %s, which is different Laue symmetry from INTEGRATE.HKL (%s)" % (symm.space_group_info(), symm_by_integrate.space_group_info())
                    need_rescale = True
            else:
                print >>decilog, "pointless using XDS_ASCII.HKL suggested %s" % symm.space_group_info()
                need_rescale = True

            if need_rescale and not flag_do_not_change_symm:
                # make backup, and do correct and compare ISa
                # if ISa got worse, revert the result.
                backup_needed = ("XDS.INP", "XDS_ASCII_fullres.HKL","CORRECT_fullres.LP",
                                 "merging_stats.pkl","merging_stats.log")
                backup_needed += xds_files.generated_by_CORRECT
                bk_prefix = make_backup(backup_needed, wdir=root, quiet=True)

                sgnum = symm.space_group_info().type().number()
                cell = " ".join(map(lambda x:"%.2f"%x, symm.unit_cell().parameters()))
                modify_xdsinp(xdsinp, inp_params=[("JOB", "CORRECT"),
                                                  ("SPACE_GROUP_NUMBER", "%d"%sgnum),
                                                  ("UNIT_CELL_CONSTANTS", cell),
                                                  ("INCLUDE_RESOLUTION_RANGE", "50 0")])

                run_xds(wdir=root, show_progress=params.show_progress)

                ret = calc_merging_stats(os.path.join(root, "XDS_ASCII.HKL"))
                
                if params.cut_resolution:
                    if ret is not None and ret[0] is not None:
                        d_min = ret[0]
                        modify_xdsinp(xdsinp, inp_params=[("JOB", "CORRECT"),
                                                          ("INCLUDE_RESOLUTION_RANGE", "50 %.2f"%d_min)])
                        print >>decilog, " Re-scale at %.2f A" % d_min
                        os.rename(os.path.join(root, "CORRECT.LP"), os.path.join(root, "CORRECT_fullres.LP"))
                        os.rename(os.path.join(root, "XDS_ASCII.HKL"), os.path.join(root, "XDS_ASCII_fullres.HKL"))
                        run_xds(wdir=root, show_progress=params.show_progress)
                        print >>decilog, " OK. ISa= %.2f" % correctlp.get_ISa(correct_lp, check_valid=True)
                        print >>decilog, " (Original files are saved as *_fullres.*)"
                    else:
                        print >>decilog, "error: Can't decide resolution."
                        for f in ("CORRECT_fullres.LP", "XDS_ASCII_fullres.HKL"):
                            if os.path.isfile(os.path.join(root, f)):
                                print >>decilog, "removing", f
                                os.remove(os.path.join(root, f))

                ISa = correctlp.get_ISa(correct_lp, check_valid=True)

                if ISa >= last_ISa or last_ISa!=last_ISa: # if improved or last_ISa is nan
                    print >>decilog, "ISa improved= %.2f" % ISa
                    remove_backups(backup_needed, bk_prefix, wdir=root)
                else:
                    print >>decilog, "ISa got worse= %.2f" % ISa
                    for f in backup_needed:
                        if os.path.isfile(os.path.join(root, f)): os.remove(os.path.join(root, f))

                    revert_files(backup_needed, bk_prefix, wdir=root, quiet=True)

    run_xdsstat(wdir=root)
    print
    if params.make_report: html_report.make_individual_report(root, root)
    print >>decilog, "xds_sequence finished at %s\n" % time.strftime("%Y-%m-%d %H:%M:%S")
    decilog.close()
# xds_sequence()

def run_xds_sequence(root, params):
    tmpdir = None
    
    if params.use_tmpdir_if_available:
        tmpdir = util.get_temp_local_dir("xdskamo", min_gb=2) # TODO guess required tempdir size
        if tmpdir is None:
            print "Can't get temp dir with sufficient size."

    # If tmpdir is not used
    if tmpdir is None:
        return xds_sequence(root, params)

    print "Using %s as temp dir.." % tmpdir

    # If tempdir is used
    for f in glob.glob(os.path.join(root, "*")): shutil.copy2(f, tmpdir)
    xdsinp = os.path.join(tmpdir, "XDS.INP")
    xdsinp_dict = dict(get_xdsinp_keyword(xdsinp))

    # Make a link to data dir
    org_data_template = xdsinp_dict["NAME_TEMPLATE_OF_DATA_FRAMES"]
    ord_data_dir = os.path.dirname(org_data_template)
    if not os.path.isabs(ord_data_dir): ord_data_dir = os.path.join(root, ord_data_dir)
    datadir_lns = os.path.join(tmpdir, "data_loc")
    os.symlink(ord_data_dir, datadir_lns)

    # Modify XDS.INP
    modify_xdsinp(xdsinp, inp_params=[("NAME_TEMPLATE_OF_DATA_FRAMES",
                                       os.path.join("data_loc", os.path.basename(org_data_template)))])

    try:
        ret = xds_sequence(tmpdir, params)
    finally:
        # Revert XDS.INP
        modify_xdsinp(xdsinp, inp_params=[("NAME_TEMPLATE_OF_DATA_FRAMES", org_data_template)])

        # Remove link
        os.remove(datadir_lns)
        
        # Move to original directory
        for f in glob.glob(os.path.join(tmpdir, "*")):
            f_dest = os.path.join(root, os.path.relpath(f, tmpdir))
            if os.path.isfile(f_dest):
                # Copy only if newer
                if os.stat(f).st_mtime > os.stat(f_dest).st_mtime:
                    shutil.copy2(f, root)
                else:
                    print "%s already exists and not modified. skip." % f

                os.remove(f)
            else:
                shutil.move(f, root)
        
        # Remove tmpdir
        shutil.rmtree(tmpdir)

    return ret
# run_xds_sequence()

class xds_runmanager(object):
    def __init__(self, params):
        self.params = params

    def __call__(self, arg):
        try:
            return run_xds_sequence(arg, self.params)
        except:
            print traceback.format_exc()

# xds_runmanager()

def run(params):
    params.topdir = os.path.abspath(params.topdir)

    xds_dirs = []
    print "Found xds directories:"
    for root, dirnames, filenames in os.walk(params.topdir, followlinks=True):
        if "XDS.INP" in filenames:
            if "decision.log" in filenames and params.dont_overwrite:
                print "Already done - skip:", os.path.relpath(root, params.topdir)
                continue
            print "", os.path.relpath(root, params.topdir)
            xds_dirs.append(root)

    print
    print "Start running.."

    import functools
    from yamtbx.dataproc.auto.command_line import run_all_xds_simple

    if params.multiproc:
        npar = util.get_number_of_processors() if params.nproc is None else params.nproc

        # Override nproc
        if len(xds_dirs) < npar: params.nproc  = npar // len(xds_dirs)
        else: params.nproc = 1
        print "nproc=", params.nproc

        if params.parmethod == "sge": npar = len(xds_dirs)

        fun = run_all_xds_simple.xds_runmanager(params)
        easy_mp.parallel_map(func=fun,
                             iterable=map(lambda x: os.path.abspath(x), xds_dirs),
                             processes=npar,
                             method=params.parmethod,
                             preserve_exception_message=True)

        """
        fun_local = lambda x: xds_sequence(x, params)
        easy_mp.pool_map(fixed_func=fun_local,
                         args=xds_dirs,
                         processes=npar)
        """
    else:
        for root in xds_dirs:
            run_xds_sequence(root, params)
                
# run()

def run_from_args(argv):
    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    for arg in args:
        if os.path.isdir(arg) and params.topdir is None:
            params.topdir = arg

    if params.topdir is None:
        params.topdir = os.getcwd()

    run(params)
# run_from_args()

if __name__ == "__main__":
    import sys

    if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
        print "Parameter syntax:\n"
        iotbx.phil.parse(master_params_str).show(prefix="  ")
        print
        print "Usage: run_all_xds_simple.py [top-dir/] [tryhard=true] [multiproc=true]"
        quit()

    run_from_args(sys.argv[1:])
