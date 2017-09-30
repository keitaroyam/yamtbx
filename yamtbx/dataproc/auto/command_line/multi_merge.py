"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.xds import xscalelp
from yamtbx.dataproc.xds import correctlp
from yamtbx.dataproc.xds.xparm import XPARM
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.dataproc.xds.command_line import xds_aniso_analysis
from yamtbx.dataproc.pointless import Pointless
from yamtbx.dataproc import aimless
from yamtbx.dataproc import xtriage
from yamtbx.dataproc.auto.command_line.run_all_xds_simple import calc_merging_stats
from yamtbx.dataproc.auto import blend
from yamtbx.dataproc.auto import cc_clustering
from yamtbx.dataproc.auto import multi_merging
from yamtbx import util
from yamtbx.util import batchjob

import iotbx.phil
import libtbx.phil
import iotbx.file_reader
import iotbx.mtz
from libtbx.utils import multi_out
from cctbx import sgtbx
from cctbx import crystal
from cctbx import miller

import os
import sys
import time
import collections
import StringIO
import networkx as nx
import traceback
import numpy
import cPickle as pickle

master_params_str = """
lstin = None
 .type = path
 .help = list of XDS_ASCII.HKL
workdir = None
 .type = path
 .help = working directory for merging data
anomalous = None
 .type = bool
d_min = None
 .type = float
d_max = None
 .type = float
clustering = *no blend cc
 .type = choice(multi=False)
program = *xscale aimless
 .type = choice(multi=False)
reject_method = delta_cc1/2 *framecc *lpstats
 .type = choice(multi=True)
 .help = Dataset rejection method
reference_file = None
 .type = path
 .help = reference (for example, low resolution but complete) XDS_ASCII.HKL.
space_group = None
 .type = str
 .help = Space group for merging
add_test_flag = False
 .type = bool
 .help = "Add test flag (FreeR_flag) in MTZ"
nproc = 1
 .type = int
 .help = number of processors that can be used.

xscale { 
 min_i_over_sigma = 3
  .type = float
  .help = MINIMUM_I/SIGMA= parameter passed to XSCALE. (criterion for calculating cc etc)
 reference = bmax bmed *bmin
  .type = choice(multi=False)
  .help = Reference selection for final scaling
 frames_per_batch = None
  .type = int(value_min=1)
  .help = affects NBATCH=. When 1, NBATCH= matches the number of frames in the dataset.
 degrees_per_batch = None
  .type = float(value_min=0.0001)
  .help = "angle version of frames_per_batch"
 corrections = *MODPIX *ABSORP *DECAY
  .type = choice(multi=True)
  .help = controls CORRECTIONS=. Use lower case to specify
 use_tmpdir_if_available = True
  .type = bool
}

rejection {
 framecc {
  method = abs *tukey
   .type = choice(multi=False)
  iqr_coeff = 1.5
   .type = float
  abs_cutoff = 0.9
   .type = float
 }
 lpstats {
  stats = *em.b em.ab pairwise_cc bfactor rfactor
   .type = choice(multi=True)
  iqr_coeff = 1.5
   .type = float
  pwcc {
   method = *abs tukey
    .type = choice(multi=False)
   iqr_coeff = 1.5
    .type = float
   abs_cutoff = 0.8
    .type = float
    .help = Minimum acceptable CC for each dataset in XSCALE (when stats=pairwise_cc)
  }
 }
 delta_cchalf {
  bin = *total outer total-then-outer
  .type = choice(multi=False)
  .help = choice of resolution bin of CC1/2 (when reject_method=delta_cc1/2)
 }
}

max_clusters = None
 .type = int
 .help = "Upper limit of clusters for merging"

blend {
 use_old_result = None
  .type = path
  .help = Directory where BLEND ran
 min_cmpl = 90
  .type = float
  .help = minimum completeness of cluster for merging
 min_acmpl = None
  .type = float
  .help = minimum anomalous completeness of cluster for merging
 min_redun = 2
  .type = float
  .help = minimum redundancy of cluster for merging
 min_aredun = None
  .type = float
  .help = minimum anomalous redundancy of cluster for merging
 max_LCV = None
  .type = float
  .help = minimum LCV of cluster for merging
 max_aLCV = None
  .type = float
  .help = minimum aLCV of cluster for merging
}

cc_clustering {
# use_old_result = None
#  .type = path
#  .help = Directory where pre-calculated results exist
 d_min = None
  .type = float
  .help = d_min for CC calculation
 b_scale = false
  .type = bool
  .help = Scale data by B-factor before CC calculation
 use_normalized = false
  .type = bool
  .help = Use normalized structure factors in CC calculation
 min_ios = None
  .type = float
  .help = minimum I/sigma for CC calculation
 min_cmpl = 90
  .type = float
  .help = minimum completeness of cluster for merging
 min_acmpl = None
  .type = float
  .help = minimum anomalous completeness of cluster for merging
 min_redun = 2
  .type = float
  .help = minimum redundancy of cluster for merging
 min_aredun = None
  .type = float
  .help = minimum anomalous redundancy of cluster for merging
 max_clheight = None
  .type = float
  .help = maximum cluster height for merging
 nproc = 1
  .type = int
}

reference {
 data = None
  .type = path
  .help = "Reference hkl data for space group and test flag transfer (optional). Reindexing must be done beforehand."
 copy_test_flag = True
  .type = bool
  .help = Copy test flag when available
}

resolution {
 estimate = False
  .type = bool
  .help = Automatically estimate resolution cutoff
 n_bins = 9
  .type = int
  .help = "Number of shells"
 cc_one_half_min = 0.5
  .type = float
  .help = "Minimum CC1/2 value at highest resolution shell"
 cc_half_tol = 0.03
  .type = float
  .help = "CC1/2 value tolerance"
}

batch {
 par_run = *deltacchalf merging
  .type = choice(multi=True)
  .help = What to run in parallel
 engine = sge *sh
  .type = choice(multi=False)
 sge_pe_name = par
  .type = str
  .help = pe name (put after -pe option)
 nproc_each = 4
  .type = int
  .help = maximum number of cores used for single data processing
 sh_max_jobs = 1
  .type = int
  .help = maximum number of concurrent jobs when engine=sh
}
"""

def merge_datasets(params, workdir, xds_files, cells, space_group):
    if not os.path.exists(workdir): os.makedirs(workdir)
    out = open(os.path.join(workdir, "merge.log"), "w")

    if params.program == "xscale":
        cycles = multi_merging.xscale.XscaleCycles(workdir, 
                                                   anomalous_flag=params.anomalous,
                                                   d_min=params.d_min, d_max=params.d_max, 
                                                   reject_method=params.reject_method,
                                                   reject_params=params.rejection,
                                                   xscale_params=params.xscale,
                                                   res_params=params.resolution,
                                                   reference_file=params.reference_file,
                                                   space_group=space_group,
                                                   ref_mtz=params.reference.data if params.reference.copy_test_flag else None,
                                                   out=out, nproc=params.nproc,
                                                   batch_params=params.batch)

        unused_files, reasons = cycles.run_cycles(xds_files)
        used_files = set(xds_files).difference(set(unused_files))

        print >>out
        print >>out, " SUMMARY "
        print >>out, "========================"
        for i, files in enumerate((used_files, unused_files)):
            print >>out, "\n%6s %4d files:\n" % (("Used", "Unused")[i], len(files))
            if len(files) == 0:
                continue

            maxlen_f = max(map(lambda f: len(os.path.relpath(f, params.workdir)), files))

            for f in files:
                cell = cells[f]
                merge_log = os.path.join(os.path.dirname(f), "merging_stats.log")
                try:
                    lines = open(merge_log).readlines()
                    resn = float(filter(lambda x:x.startswith("Resolution:"), lines)[0].split()[-1])
                    cmpl = float(filter(lambda x:x.startswith("Completeness:"), lines)[0].split()[-1].replace("%",""))
                except:
                    resn = float("nan")
                    cmpl = float("nan")

                if i == 1: # print reason
                    print >>out, "%-15s"%reasons.get(f, "unknown"),
                print >>out, ("%-"+str(maxlen_f)+"s")%os.path.relpath(f, params.workdir), cell,
                #print >>out, "ISa=%5.1f" % correctlp.get_ISa(os.path.join(os.path.dirname(f), "CORRECT.LP")),
                print >>out, "Cmpl=%3.0f%%, Resn= %.1f" % (cmpl, resn)

        ret = []
        tkvals = lambda x: (x[-1], x[0], x[-2]) # overall, inner, outer

        for i in xrange(1, cycles.get_last_cycle_number()+1):
            wd = os.path.join(workdir, "run_%.2d"%i)
            xscale_lp = os.path.join(wd, "XSCALE.LP")
            table = xscalelp.read_stats_table(xscale_lp)
            num_files = len(xscalelp.get_read_data(xscale_lp))
            xtriage_logfile = os.path.join(wd, "ccp4", "logfile.log")
            aniso =  xds_aniso_analysis.parse_logfile(os.path.join(wd, "aniso.log"))
            cellinfo = cycles.cell_info_at_cycles[i]
            ret.append([i, wd, num_files, 
                        dict(cmpl=tkvals(table["cmpl"]),
                             redundancy=tkvals(table["redundancy"]),
                             i_over_sigma=tkvals(table["i_over_sigma"]),
                             r_meas=tkvals(table["r_meas"]),
                             cc_half=tkvals(table["cc_half"]),
                             sig_ano=tkvals(table["sig_ano"]),
                             cc_ano=tkvals(table["cc_ano"]),
                             drange=tkvals(table["d_range"]),
                             lp=xscale_lp,
                             xtriage_log=xtriage.XtriageLogfile(xtriage_logfile),
                             aniso=aniso,
                             lcv=cellinfo[1],
                             alcv=cellinfo[2],
                             dmin_est=cycles.dmin_est_at_cycles.get(i, float("nan")))
                        ])

        xscale_lp = os.path.join(cycles.current_working_dir(), "XSCALE.LP")
        print >>out, "\nFinal statistics:\n"
        print >>out, xscalelp.snip_stats_table(xscale_lp)

        return ret

    elif params.program == "aimless":
        worker = Pointless()
        print >>out, "\nRunning pointless"
        runinfo = worker.run_copy(hklout="pointless.mtz", wdir=workdir,
                                  xdsin=xds_files,
                                  logout=os.path.join(workdir, "pointless.log"),
                                  tolerance=30)

        # Table of file name -> Batch range
        assert len(xds_files) == len(runinfo)
        batch_info = collections.OrderedDict(map(lambda x: (x[0], (x[1][1:3])), zip(xds_files, runinfo)))

        cycles = multi_merging.aimless.AimlessCycles(workdir, 
                                                     anomalous_flag=params.anomalous,
                                                     d_min=params.d_min, d_max=params.d_max, 
                                                     reject_method=params.reject_method,
                                                     cc_cutoff=params.rejection.lpstats.pwcc.abs_cutoff,
                                                     delta_cchalf_bin=params.rejection.delta_cchalf.bin,
                                                     mtzin=os.path.join(workdir, "pointless.mtz"),
                                                     batch_info=batch_info,
                                                     out=out, nproc=params.nproc,
                                                     nproc_each=params.batch.nproc_each,
                                                     batchjobs=None) # FIXME batchjobs
        unused_files, reasons = cycles.run_cycles(xds_files)
        used_files = set(xds_files).difference(set(unused_files))

        print >>out
        print >>out, " SUMMARY "
        print >>out, "========================"
        for i, files in enumerate((used_files, unused_files)):
            print >>out, "\n%6s %4d files:\n" % (("Used", "Unused")[i], len(files))
            if len(files) == 0:
                continue

            maxlen_f = max(map(lambda f: len(os.path.relpath(f, params.workdir)), files))

            for f in files:
                cell = cells[f]
                merge_log = os.path.join(os.path.dirname(f), "merging_stats.log")
                try:
                    lines = open(merge_log).readlines()
                    resn = float(filter(lambda x:x.startswith("Resolution:"), lines)[0].split()[-1])
                    cmpl = float(filter(lambda x:x.startswith("Completeness:"), lines)[0].split()[-1].replace("%",""))
                except:
                    resn = float("nan")
                    cmpl = float("nan")

                if i == 1: # print reason
                    print >>out, "%-15s"%reasons.get(f, "unknown"),
                print >>out, ("%-"+str(maxlen_f)+"s")%os.path.relpath(f, params.workdir), cell,
                print >>out, "ISa=%5.1f" % correctlp.get_ISa(os.path.join(os.path.dirname(f), "CORRECT.LP")),
                print >>out, "Cmpl=%3.0f%%, Resn= %.1f" % (cmpl, resn)

        aimless_log = os.path.join(cycles.current_working_dir(), "aimless.log")
        print >>out, "\nFinal statistics:\n"
        print >>out, aimless.snip_summary(aimless_log)

        # Write summary
        table = aimless.read_summary(aimless_log)

        tkvals = lambda x: (x[0], x[1], x[2]) # overall, inner, outer
        return [[cycles.get_last_cycle_number(), cycles.current_working_dir(), len(used_files),
                dict(cmpl=tkvals(table["cmpl"]),
                     redundancy=tkvals(table["redundancy"]),
                     i_over_sigma=tkvals(table["i_over_sigma"]),
                     r_meas=tkvals(table["r_meas"]),
                     cc_half=tkvals(table["cc_half"]),
                     sig_ano=(float("nan"),)*3,
                     cc_ano=tkvals(table["cc_ano"]))], ]

        #print >>out, "\nRunning aimless"
        #aimless.run_aimless(mtzin="pointless.mtz",
        #                    wdir=workdir,
        #                    anomalous=params.anomalous, d_min=params.d_min, prefix=None)

    else:
        print >>out, "Unknown program:", params.program
        return []
# merge_datasets()

def run(params):
    if os.path.isdir(params.workdir) and os.listdir(params.workdir):
        print "Directory already exists and not empty:", params.workdir
        return

    # Check parameters
    if params.program == "xscale":
        if (params.xscale.frames_per_batch, params.xscale.degrees_per_batch).count(None) == 0:
            print "ERROR! You can't specify both of xscale.frames_per_batch and xscale.degrees_per_batch"
            return
            
    
    if params.reference_file is not None and params.program != "xscale":
        print "WARNING - reference file is not used unless program=xscale."

    if not os.path.isdir(params.workdir):
        os.makedirs(params.workdir)

    if params.batch.engine == "sge":
        batchjobs = batchjob.SGE(pe_name=params.batch.sge_pe_name)
    elif params.batch.engine == "sh":
        batchjobs = batchjob.ExecLocal(max_parallel=params.batch.sh_max_jobs)
    else:
        raise "Unknown batch engine: %s" % params.batch.engine

    out = multi_out()
    out.register("log", open(os.path.join(params.workdir, "multi_merge.log"), "w"), atexit_send_to=None)
    out.register("stdout", sys.stdout)
    out.write("kamo.multi_merge started at %s\n\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
    time_started = time.time()

    print >>out, "Paramters:"
    libtbx.phil.parse(master_params_str).format(params).show(out=out, prefix=" ")
    print >>out, ""

    # XXX Not works when clustering is used..
    html_report = multi_merging.html_report.HtmlReportMulti(os.path.abspath(params.workdir))
    try: html_report.add_params(params, master_params_str)
    except: print >>out, traceback.format_exc()

    xds_ascii_files = util.read_path_list(params.lstin, only_exists=True, as_abspath=True, err_out=out)

    if not xds_ascii_files:
        print >>out, "ERROR! Cannot find (existing) files in %s." % params.lstin
        return

    if len(xds_ascii_files) < 2:
        print >>out, "ERROR! Only one file in %s." % params.lstin
        print >>out, "       Give at least two files for merging."
        return        

    cells = collections.OrderedDict()
    laues = {} # for check
    for xac in xds_ascii_files:
        try:
            symm = XDS_ASCII(xac, read_data=False).symm
        except:
            print >>out, "Error in reading %s" % xac
            print >>out, traceback.format_exc()
            return
        cells[xac] = symm.unit_cell().parameters()
        laue = symm.space_group().build_derived_reflection_intensity_group(False).info()
        laues.setdefault(str(laue),{}).setdefault(symm.space_group_info().type().number(), []).append(xac)

    if len(laues) > 1:
        print >>out, "ERROR! more than one space group included."
        for laue in laues:
            print "Laue symmetry", laue
            for sg in laues[laue]:
                print >>out, " SPACE_GROUP_NUMBER= %d (%d data)" % (sg, len(laues[laue][sg]))
                for f in laues[laue][sg]: print >>out, "  %s" % f
                print >>out, ""
        return

    space_group = None
    if params.space_group is not None:
        space_group = sgtbx.space_group_info(params.space_group).group()
        laue_given = str(space_group.build_derived_reflection_intensity_group(False).info())
        if laue_given != laues.keys()[0]:
            print >>out, "ERROR! user-specified space group (space_group=%s) is not compatible with input files (%s)" % (params.space_group, laues.keys()[0])
            return
    else:
        tmp = sgtbx.space_group_info(laues.values()[0].keys()[0]).group().build_derived_reflection_intensity_group(True)
        print >>out, "Space group for merging:", tmp.info()

    test_flag_will_be_transferred = False

    if params.reference.data is not None:
        params.reference.data = os.path.abspath(params.reference.data)
        print >>out, "Reading reference data file: %s" % params.reference.data

        tmp = iotbx.file_reader.any_file(params.reference.data, force_type="hkl", raise_sorry_if_errors=True)
        if params.reference.copy_test_flag:
            from yamtbx.dataproc.command_line import copy_free_R_flag
            if None in copy_free_R_flag.get_flag_array(tmp.file_server.miller_arrays, log_out=out):
                print >>out, " Warning: no test flag found in reference file (%s)" % params.reference.data
            else:
                test_flag_will_be_transferred = True
                print >>out, " test flag will be transferred"

        if space_group is not None:
            if space_group != tmp.file_server.miller_arrays[0].space_group():
                print >>out, " ERROR! space_group=(%s) and that of reference.data (%s) do not match." % (space_group.info(), tmp.file_server.miller_arrays[0].space_group_info())
                return
        else:
            space_group = tmp.file_server.miller_arrays[0].space_group()
            print >>out, " space group for merging: %s" % space_group.info()

    if params.add_test_flag:
        if test_flag_will_be_transferred:
            print >>out, "Warning: add_test_flag=True was set, but the flag will be transferred from the reference file given."
        else:
            from cctbx import r_free_utils

            med_cell = numpy.median(cells.values(), axis=0)
            d_min = max(params.d_min-0.2, 1.0) if params.d_min is not None else 1.5 # to prevent infinite set
            sg = space_group
            if not sg: sg = sgtbx.space_group_info(laues.values()[0].keys()[0]).group().build_derived_reflection_intensity_group(True)
            tmp = miller.build_set(crystal.symmetry(tuple(med_cell), space_group=sg), False,
                                   d_min=d_min, d_max=None)
            print >>out, "Generating test set using the reference symmetry:"
            crystal.symmetry.show_summary(tmp, out, " ")
            tmp = tmp.generate_r_free_flags(fraction=0.05,
                                            max_free=None,
                                            lattice_symmetry_max_delta=5.0,
                                            use_lattice_symmetry=True,
                                            n_shells=20)
            tmp.show_r_free_flags_info(out=out, prefix=" ")
            tmp = tmp.customized_copy(data=r_free_utils.export_r_free_flags_for_ccp4(flags=tmp.data(),
                                                                                     test_flag_value=True))

            mtz_object = tmp.as_mtz_dataset(column_root_label="FreeR_flag").mtz_object()
            test_flag_mtz = os.path.abspath(os.path.join(params.workdir, "test_flag.mtz"))
            mtz_object.write(file_name=test_flag_mtz)

            # Override the parameters
            params.reference.copy_test_flag = True
            params.reference.data = test_flag_mtz
            
    try: html_report.add_cells_and_files(cells, laues.keys()[0])
    except: print >>out, traceback.format_exc()

    data_for_merge = []
    if params.clustering == "blend":
        if params.blend.use_old_result is None:
            blend_wdir = os.path.join(params.workdir, "blend")
            os.mkdir(blend_wdir)
            blend.run_blend0R(blend_wdir, xds_ascii_files)
            print >>out, "\nRunning BLEND with analysis mode"
        else:
            blend_wdir = params.blend.use_old_result
            print >>out, "\nUsing precalculated BLEND result in %s" % params.blend.use_old_result

        blend_clusters = blend.BlendClusters(workdir=blend_wdir, d_min=params.d_min)
        summary_out = os.path.join(blend_wdir, "blend_cluster_summary.dat")
        clusters = blend_clusters.show_cluster_summary(out=open(summary_out, "w"))
        print >>out, "Clusters found by BLEND were summarized in %s" % summary_out

        if params.blend.min_cmpl is not None:
            clusters = filter(lambda x: x[3] >= params.blend.min_cmpl, clusters)
        if params.blend.min_acmpl is not None:
            clusters = filter(lambda x: x[5] >= params.blend.min_acmpl, clusters)            
        if params.blend.min_redun is not None:
            clusters = filter(lambda x: x[4] >= params.blend.min_redun, clusters)
        if params.blend.min_aredun is not None:
            clusters = filter(lambda x: x[6] >= params.blend.min_aredun, clusters)            
        if params.blend.max_LCV is not None:
            clusters = filter(lambda x: x[7] <= params.blend.max_LCV, clusters)
        if params.blend.max_aLCV is not None:
            clusters = filter(lambda x: x[8] <= params.blend.max_aLCV, clusters)

        if params.max_clusters is not None and len(clusters) > params.max_clusters:
            print >>out, "Only first %d (/%d) clusters will be merged (as specified by max_clusters=)" % (params.max_clusters, len(clusters))
            clusters = clusters[:params.max_clusters]

        if clusters:
            print >>out, "With specified conditions, following %d clusters will be merged:" % len(clusters)
        else:
            print >>out, "\nERROR: No clusters satisfied the specified conditions for merging!"
            print >>out, "Please change criteria of completeness or redundancy"
            print >>out, "Here is the table of completeness and redundancy for each cluster:\n"
            print >>out, open(summary_out).read()

        for clno, IDs, clh, cmpl, redun, acmpl, aredun, LCV, aLCV in clusters: # process largest first
            print >>out, " Cluster_%.4d NumDS= %4d CLh= %5.1f Cmpl= %6.2f Redun= %4.1f ACmpl=%6.2f ARedun=%4.1f LCV= %5.1f aLCV=%5.1f" % (clno, len(IDs), clh, cmpl, redun, acmpl, aredun, LCV, aLCV)
            data_for_merge.append((os.path.join(params.workdir, "cluster_%.4d"%clno),
                                   map(lambda x: blend_clusters.files[x-1], IDs),
                                   LCV, aLCV,clh))
        print >>out
        try: html_report.add_clutering_result(clusters, "blend")
        except: print >>out, traceback.format_exc()

    elif params.clustering == "cc":
        ccc_wdir = os.path.join(params.workdir, "cc_clustering")
        os.mkdir(ccc_wdir)
        cc_clusters = cc_clustering.CCClustering(ccc_wdir, xds_ascii_files,
                                                 d_min=params.cc_clustering.d_min if params.cc_clustering.d_min is not None else params.d_min,
                                                 min_ios=params.cc_clustering.min_ios)
        print >>out, "\nRunning CC-based clustering"

        cc_clusters.do_clustering(nproc=params.cc_clustering.nproc,
                                  b_scale=params.cc_clustering.b_scale,
                                  use_normalized=params.cc_clustering.use_normalized,
                                  html_maker=html_report)
        summary_out = os.path.join(ccc_wdir, "cc_cluster_summary.dat")
        clusters = cc_clusters.show_cluster_summary(d_min=params.d_min, out=open(summary_out, "w"))
        print >>out, "Clusters were summarized in %s" % summary_out

        if params.cc_clustering.min_cmpl is not None:
            clusters = filter(lambda x: x[3] >= params.cc_clustering.min_cmpl, clusters)
        if params.cc_clustering.min_acmpl is not None:
            clusters = filter(lambda x: x[5] >= params.cc_clustering.min_acmpl, clusters)            
        if params.cc_clustering.min_redun is not None:
            clusters = filter(lambda x: x[4] >= params.cc_clustering.min_redun, clusters)
        if params.cc_clustering.min_aredun is not None:
            clusters = filter(lambda x: x[6] >= params.cc_clustering.min_aredun, clusters)            
        if params.cc_clustering.max_clheight is not None:
            clusters = filter(lambda x: x[2] <= params.cc_clustering.max_clheight, clusters)

        if params.max_clusters is not None and len(clusters) > params.max_clusters:
            print >>out, "Only first %d (/%d) clusters will be merged (as specified by max_clusters=)" % (params.max_clusters, len(clusters))
            clusters = clusters[:params.max_clusters]

        if clusters:
            print >>out, "With specified conditions, following %d clusters will be merged:" % len(clusters)
        else:
            print >>out, "\nERROR: No clusters satisfied the specified conditions for merging!"
            print >>out, "Please change criteria of completeness or redundancy"
            print >>out, "Here is the table of completeness and redundancy for each cluster:\n"
            print >>out, open(summary_out).read()

        for clno, IDs, clh, cmpl, redun, acmpl, aredun, ccmean, ccmin in clusters: # process largest first
            print >>out, " Cluster_%.4d NumDS= %4d CLh= %5.1f Cmpl= %6.2f Redun= %4.1f ACmpl=%6.2f ARedun=%4.1f CCmean=% .4f CCmin=% .4f" % (clno, len(IDs), clh, cmpl, redun, acmpl, aredun, ccmean, ccmin)
            data_for_merge.append((os.path.join(params.workdir, "cluster_%.4d"%clno),
                                   map(lambda x: xds_ascii_files[x-1], IDs),
                                   float("nan"),float("nan"),clh))
        print >>out

        try: html_report.add_clutering_result(clusters, "cc_clustering")
        except: print >>out, traceback.format_exc()
        
    else:
        data_for_merge.append((os.path.join(params.workdir, "all_data"),
                               xds_ascii_files, float("nan"), float("nan"), 0))

    ofs_summary = open(os.path.join(params.workdir, "cluster_summary.dat"), "w")
    ofs_summary.write("# d_min= %.3f A\n" % (params.d_min if params.d_min is not None else float("nan")))
    ofs_summary.write("# LCV and aLCV are values of all data\n")
    ofs_summary.write("     cluster    ClH  LCV aLCV run ds.all ds.used  Cmpl Redun I/sigI Rmeas CC1/2 Cmpl.ou Red.ou I/sig.ou Rmeas.ou CC1/2.ou Cmpl.in Red.in I/sig.in Rmeas.in CC1/2.in SigAno.in CCano.in WilsonB Aniso.bst Aniso.wst dmin.est\n")

    out.flush()

    def write_ofs_summary(workdir, cycle, clh, LCV, aLCV, xds_files, num_files, stats):
        tmps = "%12s %6.2f %4.1f %4.1f %3d %6d %7d %5.1f %5.1f %6.2f %5.1f %5.1f %7.1f %6.1f % 8.2f % 8.1f %8.1f %7.1f %6.1f % 8.2f % 8.1f %8.1f %9.1f %8.1f %7.2f %9.2f %9.2f %.2f\n"
        ofs_summary.write(tmps % (os.path.relpath(workdir, params.workdir), clh, LCV, aLCV, cycle,
                                  len(xds_files), num_files,
                                  stats["cmpl"][0],
                                  stats["redundancy"][0],
                                  stats["i_over_sigma"][0],
                                  stats["r_meas"][0],
                                  stats["cc_half"][0],
                                  stats["cmpl"][2],
                                  stats["redundancy"][2],
                                  stats["i_over_sigma"][2],
                                  stats["r_meas"][2],
                                  stats["cc_half"][2],
                                  stats["cmpl"][1],
                                  stats["redundancy"][1],
                                  stats["i_over_sigma"][1],
                                  stats["r_meas"][1],
                                  stats["cc_half"][1],
                                  stats["sig_ano"][1],
                                  stats["cc_ano"][1],
                                  stats["xtriage_log"].wilson_b,
                                  #stats["xtriage_log"].anisotropy,
                                  stats["aniso"]["d_min_best"],
                                  stats["aniso"]["d_min_worst"],
                                  stats["dmin_est"],
                                  ))
        ofs_summary.flush()
    # write_ofs_summary()

    if "merging" in params.batch.par_run:
        params.nproc = params.batch.nproc_each
        jobs = []
        for workdir, xds_files, LCV, aLCV, clh in data_for_merge:
            if not os.path.exists(workdir): os.makedirs(workdir)
            shname = "merge_%s.sh" % os.path.relpath(workdir, params.workdir)
            pickle.dump((params, os.path.abspath(workdir), xds_files, cells, space_group), open(os.path.join(workdir, "args.pkl"), "w"), -1)
            job = batchjob.Job(workdir, shname, nproc=params.batch.nproc_each)
            job.write_script("""\
cd "%s" || exit 1
"%s" -c '\
import pickle; \
from yamtbx.dataproc.auto.command_line.multi_merge import merge_datasets; \
args = pickle.load(open("args.pkl")); \
ret = merge_datasets(*args); \
pickle.dump(ret, open("result.pkl","w")); \
'
""" % (os.path.abspath(workdir), sys.executable))
            batchjobs.submit(job)
            jobs.append(job)

        batchjobs.wait_all(jobs)
        for workdir, xds_files, LCV, aLCV, clh in data_for_merge:
            try:
                results = pickle.load(open(os.path.join(workdir, "result.pkl")))
            except:
                print >>out, "Error in unpickling result in %s" % workdir
                print >>out, traceback.format_exc()
                results = []

            if len(results) == 0:
                ofs_summary.write("#%s failed\n" % os.path.relpath(workdir, params.workdir))

            lcv, alcv = float("nan"), float("nan")
            for cycle, wd, num_files, stats in results:
                lcv, alcv = stats.get("lcv", LCV), stats.get("alcv", aLCV)
                write_ofs_summary(workdir, cycle, clh, lcv, alcv, xds_files, num_files, stats)

            # Last lcv & alcv
            try: html_report.add_merge_result(workdir, clh, lcv, alcv, xds_files, results[-1][2], results[-1][3])
            except: print >>out, traceback.format_exc()
    else:
        for workdir, xds_files, LCV, aLCV, clh in data_for_merge:
            print >>out, "Merging %s..." % os.path.relpath(workdir, params.workdir)
            out.flush()
            results = merge_datasets(params, workdir, xds_files, cells, space_group)
            
            if len(results) == 0:
                ofs_summary.write("#%s failed\n" % os.path.relpath(workdir, params.workdir))

            for cycle, wd, num_files, stats in results:
                lcv, alcv = stats.get("lcv", LCV), stats.get("alcv", aLCV)
                write_ofs_summary(workdir, cycle, clh, lcv, alcv, xds_files, num_files, stats)

            try: html_report.add_merge_result(workdir, clh, lcv, alcv, xds_files, results[-1][2], results[-1][3])
            except: print >>out, traceback.format_exc()

    try: html_report.write_html()
    except: print >>out, traceback.format_exc()

    print "firefox %s" % os.path.join(html_report.root, "report.html")

    out.write("\nNormal exit at %s\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
    out.write("Total wall-clock time: %.2f sec.\n" % (time.time()-time_started))

    return
# run()

def run_from_args(argv):
    if "-h" in argv or "--help" in argv:
        print """
kamo.multi_merge is a program for merging multiple (small wedge) datasets.
This is an alpha-version. If you found something wrong, please let staff know! We would appreciate your feedback.

* Use cases (options) *

1) Use BLEND clustering + default outlier rejections (using XSCALE)

  workdir=blend_3A_framecc_b lstin=formerge.lst d_min=3.0 anomalous=false clustering=blend blend.min_cmpl=90

2) Use CC-based clustering instead

  workdir=ccc_3A_framecc_b lstin=formerge.lst d_min=3.0 anomalous=false clustering=cc cc_clustering.d_min=3.5 cc_clustering.min_cmpl=90

All parameters:
"""
        iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
        return

    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args
    flag_unrecognized = False
    
    for arg in args:
        if os.path.isdir(arg) and params.topdir is None:
            params.topdir = arg
        elif os.path.isfile(arg) and params.lstin is None:
            params.lstin = arg
        else:
            print "ERROR: unrecognized arg:", arg
            flag_unrecognized = True

    if flag_unrecognized:
        return

    if params.lstin is None:
        print "Give lstin="
        print "Use -h option to see help"
        return

    params.lstin = os.path.abspath(params.lstin)

    if params.workdir is None:
        params.workdir = time.strftime("merge_%y%m%d-%H%M%S")

    run(params)
# run_from_args()

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
