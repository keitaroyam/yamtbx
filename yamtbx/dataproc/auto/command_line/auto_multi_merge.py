"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import libtbx.phil
import iotbx.phil
import iotbx.file_reader
from cctbx.crystal import reindex
from cctbx.array_family import flex
from cctbx import miller
from iotbx import crystal_symmetry_from_any
from libtbx.utils import multi_out
from yamtbx.util import batchjob
from yamtbx.util import replace_forbidden_chars
from yamtbx.util.xtal import format_unit_cell
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.dataproc.auto.command_line import multi_check_cell_consistency
from yamtbx.dataproc.auto.command_line import multi_merge
from yamtbx.dataproc.auto.command_line.multi_prep_merging import rescale_with_specified_symm, reindex_with_specified_symm
from yamtbx.dataproc.auto.multi_merging import resolve_reindex
from yamtbx.dataproc.auto.resolution_cutoff import estimate_resolution_based_on_cc_half
import os
import glob
import pickle
import sys
import csv
import copy
import collections
import time
import subprocess
import traceback

gui_phil_str = """\
csv = None
 .type = path
 .help = CSV file to define sample information
workdir = None
 .type = path
 .help = top directory where merging directories are created. current directory by default.
datadir = None
 .type = path
 .help = the root directory where the processed results exist.
prefix = merge_
 .type = str
 .help = Prefix of directory names
postrefine = False
 .type = bool
 .help = The method to determine the unit cell for each data. If true, do post-refinement with the constraints by symmetry. Otherwise just use transformed P1 cell.
reference = None
 .type = path
 .help = Reference for symmetry and reindexing

merge {
  d_min_start = 1.8
   .type = float
   .help = Starting value for merging
  %s
}

rescut {
 auto = true
  .type = bool
 n_bins = 9
  .type = int
 cc_one_half_min = 0.5
  .type = float
 cc_half_tol = 0.03
  .type = float
}

batch {
 engine = sge sh *no
  .type = choice(multi=False)
 sge_pe_name = par
  .type = str
  .help = pe name (put after -pe option)
 nproc_each = 1
  .type = int
  .help = maximum number of cores used for single data processing
}
""" % multi_merge.master_params_str

def read_sample_info(csvin, datadir=None):
    reader = csv.reader(open(csvin, "rU"))
    header = map(lambda x: x.strip(), next(reader))
    hidxes = dict(map(lambda x:(x[1],x[0]), enumerate(header)))
    puck_flag = set(["uname","puck","pin","name"]).issubset(header)

    ret = collections.OrderedDict()

    if not puck_flag and not set(["name","topdir"]).issubset(header):
        print "Error! CSV file is invalid (header is not good)"
        return ret

    for vals in reader:
        if not vals: continue

        if puck_flag:
            if "root_dir" in hidxes: rdir = vals[hidxes["root_dir"]].strip()
            else: rdir = datadir

            if datadir is None:
                raise Exception("Provide datadir= parameter or give root_dir column in csv file!")

            uname = vals[hidxes["uname"]].strip()
            puck = vals[hidxes["puck"]].strip()
            pin = int(vals[hidxes["pin"]])
            ddirs = glob.glob(os.path.join(rdir, "%s-%s-%.2d" % (uname, puck, pin)))
        else:
            tmp = os.path.expanduser(vals[hidxes["topdir"]].strip())
            if "*" in tmp: ddirs = glob.glob(tmp)
            else: ddirs = [tmp]

        sample_name = vals[hidxes["name"]].strip()

        if sample_name not in ret: ret[sample_name] = [[], {}]

        ret[sample_name][0].extend(ddirs)

        # Custom parameters
        if "anomalous" in hidxes:
            tmp = vals[hidxes["anomalous"]].strip().lower()
            if tmp != "":
                assert tmp in ("yes", "no")
                ret[sample_name][1]["anomalous"] = (tmp=="yes")
        if "reference" in hidxes:
            tmp = vals[hidxes["reference"]].strip()
            if tmp != "":
                ret[sample_name][1]["reference"] = tmp
            
    return ret
# read_sample_info()

def choose_best_result(summarydat, log_out):
    wdir = os.path.dirname(summarydat)

    lines = filter(lambda x: not x.startswith("#"), open(summarydat))
    
    header = lines[0].split()
    i_cchalf = header.index("CC1/2")
    i_cchalf_ou = header.index("CC1/2.ou")
    i_redun = header.index("Redun")
    i_cls = header.index("cluster")
    i_run = header.index("run")

    results = []
    cls_runs = {}
    
    for l in lines[1:]:
        if l.startswith("#"): continue
        sp = l.split()
        run = int(sp[i_run])
        hklfile = os.path.join(wdir, sp[i_cls], "run_%.2d" % run, "xscale.hkl")
        results.append((hklfile, sp[i_cls], run, # 0,1,2
                        float(sp[i_cchalf]), float(sp[i_cchalf_ou]), float(sp[i_redun]) # 3,4,5
                        ))
        cls_runs.setdefault(sp[i_cls], []).append(run)

    if not results: return None

    # Remove non-final runs
    results = filter(lambda x: x[2]==max(cls_runs[x[1]]), results)

    results.sort(key=lambda x:x[5], reverse=True)
    results = results[:len(results)//2] # First half of top redundancy

    results.sort(key=lambda x:(x[3], x[4]), reverse=True)
    best_result = results[0][0]
    #open(os.path.join(os.path.dirname(best_result), "THIS_MAY_BE_THE_BEST"), "w")
    log_out.write("The best result chosen from %s is %s (CC1/2=%s CC1/2.ou=%s Red=%s)" % (summarydat, best_result,
                                                                                          results[0][3], results[0][4], results[0][5]))
    return best_result
# choose_best_result()

def decide_resolution(summarydat, params, log_out):
    best = choose_best_result(summarydat, log_out)
    if best is None:
        log_out.write("No data for deciding resolution cutoff.\n")
        return None

    log_out.write("Using %s for deciding resolution cutoff.\n" % best)
    iobs = XDS_ASCII(best, i_only=True).i_obs() # Result with max CC1/2

    est = estimate_resolution_based_on_cc_half(iobs, params.cc_one_half_min, params.cc_half_tol, params.n_bins, log_out=log_out)
    if None not in (est.d_min, est.cc_at_d_min):
        log_out.write("Best resolution cutoff= %.2f A @CC1/2= %.4f\n" % (est.d_min, est.cc_at_d_min))
    else:
        log_out.write("Can't decide resolution cutoff. No reflections??\n")
    return est.d_min
# decide_resolution()

def auto_merge(workdir, topdirs, do_postrefine, ref_array, merge_params, rescut_params, log_out_all=None):
    log_out = multi_out()
    log_out.register("log", open(os.path.join(workdir, "multi_merge.log"), "w"))
    if log_out_all: log_out.register("original", log_out_all)

    log_out.write("Inspecting data for %s\n" % workdir)
    cell_params = libtbx.phil.parse(input_string=multi_check_cell_consistency.master_params_str).extract()
    cell_params.xdsdir = []
    for topdir in topdirs:
        for root, dirnames, filenames in os.walk(topdir):
            if "XDS.INP" in filenames: cell_params.xdsdir.append(root)

    cell_params.xdsdir.sort()
    log_out.write("NOTE: %6d xds directories detected.\n" % len(cell_params.xdsdir))

    #if config.params.merging.cell_grouping.tol_length is not None:
    #    cell_params.tol_length = config.params.merging.cell_grouping.tol_length
    #if config.params.merging.cell_grouping.tol_angle is not None:
    #    cell_params.tol_angle = config.params.merging.cell_grouping.tol_angle 

    cm = multi_check_cell_consistency.run(cell_params, out=log_out)
    if len(cm.groups) == 0:
        log_out.write("Error: No data!\n")
        return

    if ref_array:
        group = 1 # TODO select group id by comparing to reference cell
        reference_symm = cm.get_symmetry_reference_matched(group-1, ref_array)
    else:
        group = 1
        reference_symm = cm.get_most_frequent_symmetry(group-1)
        if reference_symm is None:
            raise "Failed to get reference symmetry"

    log_out.write("\n\nNOTE: group_choice= %d, symmetry= %s (%s)\n" % (group, reference_symm.space_group_info(),
                                                                 format_unit_cell(reference_symm.unit_cell())))

    symms = map(lambda i: cm.symms[i], cm.groups[group-1])
    dirs = map(lambda i: cm.dirs[i], cm.groups[group-1])
    topdir = os.path.dirname(os.path.commonprefix(dirs))

    if do_postrefine:
        cell_and_files, reference_symm = rescale_with_specified_symm(topdir, dirs, symms, reference_symm=reference_symm, out=log_out)
    else:
        cell_and_files = reindex_with_specified_symm(topdir, reference_symm, dirs, out=log_out)

    log_out.flush()

    log_out.write("NOTE: %6d datasets for merging.\n" % len(cell_and_files))

    # Make list for merging
    lstname = os.path.join(workdir, "formerge.lst")
    ofs_lst = open(lstname, "w")
    ofs_dat = open(os.path.join(workdir, "cells.dat"), "w")
    ofs_dat.write("file a b c al be ga\n")
    for wd in cell_and_files:
        cell, xas = cell_and_files[wd]
        ofs_lst.write(xas+"\n")
        ofs_dat.write(xas+" "+" ".join(map(lambda x:"%7.3f"%x, cell))+"\n")

    ofs_dat.close()
    ofs_lst.close()

    cosets = reindex.reindexing_operators(reference_symm, reference_symm, max_delta=5)
    reidx_ops = cosets.combined_cb_ops()
    if len(reidx_ops) > 1:
        log_out.write("!! ATTENTION !! Reindex operators found.\n")
        for rop in reidx_ops:
            log_out.write(" operator: %-16s Cell: (%s)\n" % (rop.as_hkl(),
                                                             format_unit_cell(reference_symm.unit_cell().change_basis(rop))))

        if ref_array and ref_array.size() > 0:
            rb = resolve_reindex.ReferenceBased(map(lambda x: x[1], cell_and_files.values()),
                                                ref_array, max_delta=5,
                                                d_min=3, min_ios=None,
                                                nproc=1, log_out=log_out)
            rb.assign_operators()
        elif len(cell_and_files) > 1: # no need if only 1 data
            rb = resolve_reindex.KabschSelectiveBreeding(map(lambda x: x[1], cell_and_files.values()), max_delta=5,
                                                         d_min=3, min_ios=None,
                                                         nproc=1, log_out=log_out)
            rb.assign_operators()#max_cycle=params.max_cycles)

        lstname = os.path.join(workdir, "formerge_reindexed.lst")
        new_files = rb.modify_xds_ascii_files(cells_dat_out=open(os.path.join(workdir, "cells_reindexed.dat"), "w"))
        open(lstname, "w").write("\n".join(new_files)+"\n")

        # TODO Put citation!!

    if len(cell_and_files) < 2:
        log_out.write("Too few datasets. giveup.\n")
        return

    # Start merging
    merge_params.lstin = lstname
    merge_params.d_min = merge_params.d_min_start

    merge_params.workdir = os.path.join(workdir, "%s_%.2fA"%(merge_params.clustering, merge_params.d_min))
    multi_merge.run(merge_params)

    if rescut_params.auto:
        for cc_cut in (rescut_params.cc_one_half_min*.7, rescut_params.cc_one_half_min):
            rescut_params.cc_one_half_min = cc_cut
            rescut = decide_resolution(os.path.join(merge_params.workdir, "cluster_summary.dat"), rescut_params, log_out)
            if rescut:
                tmp = os.path.join(workdir, "%s_%.2fA"%(merge_params.clustering, rescut))
                if os.path.exists(tmp): break

                merge_params.workdir = tmp
                merge_params.d_min = rescut
                multi_merge.run(merge_params)
                choose_best_result(os.path.join(merge_params.workdir, "cluster_summary.dat"), log_out)
            else:
                break

    os.rename(merge_params.workdir, merge_params.workdir+"_final")
    log_out.flush()
# auto_merge()

def read_reference_data(filename, log_out):
    log_out.write("Reading reference file: %s\n" % filename)
    f = iotbx.file_reader.any_file(filename)
    if f.file_type == "hkl":
        array = f.file_server.get_xray_data(file_name=None,
                                            labels=None,
                                            ignore_all_zeros=True,
                                            parameter_scope="",
                                            prefer_anomalous=False,
                                            prefer_amplitudes=False)
        log_out.write(" Read data: %s\n\n" % array.info())
        return array.as_intensity_array()
    elif f.file_type == "pdb":
        import mmtbx.utils
        xrs = f.file_content.xray_structure_simple()
        fmodel_params = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params.extract()
        fmodel_params.fmodel.k_sol = 0.35
        fmodel_params.fmodel.b_sol = 50
        fmodel_params.high_resolution = 3.0
        log_out.write(" Constructed reference from model:\n")
        xrs.show_summary(log_out, "  ")
        log_out.write("\n")
        return mmtbx.utils.fmodel_from_xray_structure(xray_structure=xrs, params=fmodel_params).f_model.as_intensity_array()
    else:
        xs = crystal_symmetry_from_any.extract_from(filename)
        if xs is None:
            raise "Not useful file: %s" % filename
        log_out.write(" Read symmetry reference: \n")
        xs.show_summary(log_out, "  ")
        log_out.write("\n")
        return miller.array(miller.set(xs, flex.miller_index()), flex.double()) # dummy array
# read_reference_data()

def run(params):
    if not os.path.exists(params.workdir): os.makedirs(params.workdir)
    log_out = multi_out()
    log_out.register("log", open(os.path.join(params.workdir, time.strftime("automerge_%y%m%d-%H%M%S.log")), "w"), atexit_send_to=None)
    log_out.register("stdout", sys.stdout)

    log_out.write("Paramters:\n")
    libtbx.phil.parse(gui_phil_str).format(params).show(out=log_out, prefix=" ")
    log_out.write("\n")


    ref_arrays = {}

    if params.reference:
        ref_arrays[params.reference] = read_reference_data(params.reference, log_out)

    samples = read_sample_info(params.csv, params.datadir)

    log_out.write("Loaded from %s\n"%params.csv)
    for k in samples:
        log_out.write(" %s\n"%k)
        log_out.write("  # Overridden parameters: %s\n"%samples[k][1])
        for d in samples[k][0]:
            log_out.write("  %s\n"%d)
        log_out.write("\n")

    log_out.flush()

    # Load custom reference data
    for k in samples:
        if not "reference" in samples[k][1]: continue
        ref_filename = samples[k][1]["reference"]
        if ref_filename in ref_arrays: continue
        ref_arrays[ref_filename] = read_reference_data(ref_filename, log_out)

    if params.batch.engine == "sge":
        batchjobs = batchjob.SGE(pe_name=params.batch.sge_pe_name)
    elif params.batch.engine == "sh":
        batchjobs = batchjob.ExecLocal()
    else:
        batchjobs = None

    jobs = []

    for k in samples:
        params2 = copy.deepcopy(params)
        ref_array = ref_arrays.get(params.reference, None)
        params2.merge.reference.data = params.reference
            

        # Reflect custom parameters
        if "anomalous" in samples[k][1]:
            params2.merge.anomalous = samples[k][1]["anomalous"]
        if "reference" in samples[k][1]:
            ref_array = ref_arrays[samples[k][1]["reference"]]
            params2.merge.reference.data = samples[k][1]["reference"]

        if params2.merge.reference.data:
            params2.merge.reference.data = os.path.abspath(params2.merge.reference.data)

        log_out.write("\n\n")
        workdir = os.path.join(params2.workdir, "%s%s" % (params2.prefix,
                                                         replace_forbidden_chars(k).replace(" ","_")))
        os.mkdir(workdir)
        if batchjobs:
            shname = "multimerge.sh"
            pickle.dump(dict(workdir=workdir, topdirs=samples[k][0],
                             do_postrefine=params2.postrefine, ref_array=ref_array,
                             merge_params=params2.merge, rescut_params=params2.rescut),
                        open(os.path.join(workdir, "kwargs.pkl"), "w"), -1)
            job = batchjob.Job(workdir, shname, nproc=params2.batch.nproc_each)
            job.write_script("""\
cd "%s" || exit 1
"%s" -c '\
import pickle; \
from yamtbx.dataproc.auto.command_line.auto_multi_merge import auto_merge; \
kwargs = pickle.load(open("kwargs.pkl")); \
auto_merge(**kwargs); \
'
""" % (os.path.abspath(workdir), sys.executable))
            batchjobs.submit(job)
            jobs.append(job)
        else:
            try:
                auto_merge(workdir=workdir, topdirs=samples[k][0],
                           do_postrefine=params2.postrefine, ref_array=ref_array,
                           merge_params=params2.merge, rescut_params=params2.rescut,
                           log_out_all=log_out)
            except:
                log_out.write("Error occurred in %s\n%s\n"%(workdir, traceback.format_exc()))

        log_out.flush()

    if batchjobs:
        batchjobs.wait_all(jobs)

# run()
def run_from_args(argv):
    if "-h" in argv or "--help" in argv or not argv:
        print "All parameters:\n"
        iotbx.phil.parse(gui_phil_str).show(prefix="  ", attributes_level=1)
        print
        print
        print r"""Typical usage:

kamo.auto_multi_merge.dev \
  csv=automerge.csv \
  workdir=$PWD \
  prefix=merge_ \
  postrefine=False \
  merge.max_clusters=15 \
  merge.d_min_start=1.4 \
  merge.clustering=cc \
  merge.cc_clustering.min_acmpl=90 \
  merge.cc_clustering.min_aredun=2 \
  batch.engine=sge \
  merge.batch.engine=sge \
  merge.batch.par_run=merging \

CSV file should look like (first line must be a header line):

topdir,name,anomalous,reference
/home/hoge/_kamoproc/01, sample1, yes, xxxx.mtz
/home/hoge/_kamoproc/02, sample1, yes, xxxx.mtz
...

where anomalous,reference columns are optional.

Alternatively, it can look like (for Zoo mode):

root_dir,uname,puck,pin,name,anomalous,reference
/home/hoge/zoo/_kamoproc/, user1, CPS-0000, 1, sample1, yes, xxxx.mtz
/home/hoge/zoo/_kamoproc/, user1, CPS-0000, 2, sample1, yes, xxxx.mtz
...
"""

        return

    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=gui_phil_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    if params.workdir is None:
        params.workdir = os.getcwd()

    for arg in args:
        if not os.path.exists(arg):
            print "Error: Given path does not exist: %s" % arg
            return
        if params.csv is None and arg.endswith(".csv"):
            params.csv = arg

    run(params)
# run_from_args()

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
