"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
import libtbx.phil
import iotbx.phil
import iotbx.file_reader
from cctbx import crystal
from cctbx.crystal import reindex
from cctbx.array_family import flex
from cctbx import miller
from iotbx import crystal_symmetry_from_any
from libtbx.utils import multi_out
from yamtbx.util import batchjob
from yamtbx.util import replace_forbidden_chars
from yamtbx.util.xtal import format_unit_cell
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.dataproc.auto.command_line.multi_check_cell_consistency import CellGraph
from yamtbx.dataproc.auto.command_line import multi_merge
from yamtbx.dataproc.auto.command_line.multi_prep_merging import PrepMerging
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
import numpy

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
cell_method = *reindex refine
 .type = choice(multi=False)
 .help = "The method to determine the unit cell for each data. If refine, do refinement with the constraints by symmetry. Otherwise just use transformed P1 cell."
reference = None
 .type = path
 .help = Reference for symmetry and reindexing

space_group = None
 .type = str
 .help = Preferred space group (need unit_cell also)
unit_cell = None
 .type = floats(6)
 .help = Preferred unit cell parameters (need space_group also)

merge {
  d_min_start = 1.8
   .type = float
   .help = Starting value for merging
  %s
}

filtering {
  choice = cell
   .type = choice(multi=True)
  cell_iqr_scale = 2.5
   .type = float
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
 sh_max_jobs = 1
  .type = int
  .help = maximum number of concurrent jobs when engine=sh
}
""" % multi_merge.master_params_str

def read_sample_info(csvin, datadir=None):
    reader = csv.reader(open(csvin, "rU"))
    header = [x.strip().lower() for x in next(reader)]
    hidxes = dict([(x[1],x[0]) for x in enumerate(header)])
    puck_flag = set(["uname","puck","pin","name"]).issubset(header)

    ret = collections.OrderedDict()

    if not puck_flag and not set(["name","topdir"]).issubset(header):
        print("Error! CSV file is invalid (header is not good)")
        return ret

    if len(set(["space_group","unit_cell"]).intersection(header)) == 1:
        print("Error! Both space_group and unit_cell should be specified if you have preferred cell.")
        return ret

    for vals in reader:
        if not vals: continue

        if puck_flag:
            if "root_dir" in hidxes: rdir = vals[hidxes["root_dir"]].strip()
            else: rdir = datadir

            if rdir is None:
                raise Exception("Provide datadir= parameter or give root_dir column in csv file!")

            uname = vals[hidxes["uname"]].strip()
            puck = vals[hidxes["puck"]].strip()
            pin = int(vals[hidxes["pin"]])
            ddirs = glob.glob(os.path.join(rdir, "%s-%s-%.2d" % (uname, puck, pin)))
        else:
            tmp = os.path.expanduser(vals[hidxes["topdir"]].strip())
            ddirs = glob.glob(tmp)

        sample_name = vals[hidxes["name"]].strip()

        if sample_name not in ret: ret[sample_name] = [[], {}]

        ret[sample_name][0].extend(ddirs)

        # Custom parameters
        # XXX if specified more than once for a sample?
        if "anomalous" in hidxes:
            tmp = vals[hidxes["anomalous"]].strip().lower()
            if tmp != "":
                assert tmp in ("yes", "no")
                ret[sample_name][1]["anomalous"] = (tmp=="yes")
        if "reference" in hidxes:
            tmp = vals[hidxes["reference"]].strip()
            if tmp != "":
                ret[sample_name][1]["reference"] = tmp

        if "unit_cell" in hidxes and "space_group" in hidxes:
            tmp1 = vals[hidxes["unit_cell"]].strip()
            tmp2 = vals[hidxes["space_group"]].strip()
            if tmp1 and tmp2:
                try:
                    xs = crystal.symmetry(tmp1, tmp2)
                except:
                    raise Exception("Invalid space_group and/or unit_cell (space_group=%s unit_cell=%s)"%(tmp2, tmp1))
                
                ret[sample_name][1]["reference_sym"]  = xs
            
    return ret
# read_sample_info()

def choose_best_result(summarydat, log_out):
    wdir = os.path.dirname(summarydat)

    lines = [x for x in open(summarydat) if not x.startswith("#")]
    
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
    results = [x for x in results if x[2]==max(cls_runs[x[1]])]

    results.sort(key=lambda x:x[5], reverse=True)
    if len(results) > 2:
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

def filter_datasets(cell_and_files, params, log_out): # Not tested!!
    cell_and_files = copy.copy(cell_and_files)
    deleted = {}
    if "cell" in params.choice and len(cell_and_files) > 2:
        iqrc = params.cell_iqr_scale
        log_out.write("\nCell-based filtering (iqr_coeffs=%.3f):\n"%iqrc)   #Filtering: %d files were excluded.\n" % len(deleted))
        cells = numpy.array([x[0] for x in list(cell_and_files.values())])
        q25, q50, q75 = numpy.percentile(cells, [25, 50, 75], axis=0)
        iqr = q75-q25
        lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
        log_out.write("  median: %6.2f %6.2f %6.2f %5.1f %5.1f %5.1f\n"%tuple(q50))
        log_out.write("     IQR: %6.2f %6.2f %6.2f %5.1f %5.1f %5.1f\n"%tuple(iqr))
        log_out.write("  lowlim: %6.2f %6.2f %6.2f %5.1f %5.1f %5.1f\n"%tuple(lowlim))
        log_out.write(" highlim: %6.2f %6.2f %6.2f %5.1f %5.1f %5.1f\n"%tuple(highlim))
        bad_idxes = numpy.where(numpy.any(cells>highlim, axis=1) | numpy.any(cells<lowlim, axis=1))[0]
        log_out.write(" %d files were excluded.\n"%len(bad_idxes))

        keys = list(cell_and_files.keys())
        for i in bad_idxes:
            deleted[keys[i]] = cell_and_files[keys[i]]
            del cell_and_files[keys[i]]
        
    return cell_and_files, deleted
# filter_datasets()

def auto_merge(workdir, topdirs, cell_method, ref_array, ref_sym, merge_params, rescut_params, filter_params, log_out_all=None):
    log_out = multi_out()
    log_out.register("log", open(os.path.join(workdir, "multi_merge.log"), "w"))
    if log_out_all: log_out.register("original", log_out_all)

    log_out.write("Inspecting data for %s\n" % workdir)

    xdsdirs = []
    for topdir in topdirs:
        for root, dirnames, filenames in os.walk(topdir, followlinks=True):
            if any(y.startswith("XDS_ASCII.HKL") or y == "DIALS.HKL" for y in filenames): xdsdirs.append(root)

    log_out.write("NOTE: %6d xds directories detected.\n" % len(xdsdirs))
    log_out.flush()
    
    cm = CellGraph()
    for i, xd in enumerate(xdsdirs):
        cm.add_proc_result(i, xd)

    topdir = os.path.dirname(os.path.commonprefix(topdirs))
    lstname = os.path.join(workdir, "formerge.lst")
    
    pm = PrepMerging(cm)
    pm.find_groups()
    
    if len(cm.groups) == 0:
        log_out.write("Error: No data!\n")
        return

    if ref_array or ref_sym:
        group = cm.get_group_symmetry_reference_matched(ref_array if ref_array else ref_sym)
        reference_symm = cm.get_symmetry_reference_matched(group-1, ref_array if ref_array else ref_sym)
    else:
        group = 1
        reference_symm = cm.get_most_frequent_symmetry(group-1)
        if reference_symm is None:
            raise "Failed to get reference symmetry"

    # TODO This writes formerge.lst but `lstname` should be passed to this function.
    msg, reidx_ops = pm.prep_merging(workdir=workdir, group=group, reference_symm=reference_symm,
                                     topdir=topdir, cell_method=cell_method,
                                     nproc=1, prep_dials_files=False, into_workdir=True)
    log_out.write(msg+"\n")
    cell_and_files = pm.cell_and_files

    if filter_params.choice:
        cell_and_files, deleted = filter_datasets(cell_and_files, filter_params, log_out)
        with open(os.path.join(workdir, "excluded.lst"), "w") as ofs:
            deleted_files = [deleted[x][1] for x in deleted]
            deleted_files.sort()
            ofs.write("\n".join(deleted_files)+"\n")

        with open(lstname, "w") as ofs: # Overwrite lstname
            for wd in sorted(cell_and_files):
                _, xas = cell_and_files[wd]
                ofs.write(xas+"\n")


    if len(reidx_ops) > 1:
        if ref_array and ref_array.size() > 0:
            rb = resolve_reindex.ReferenceBased([x[1] for x in list(cell_and_files.values())],
                                                ref_array, max_delta=5,
                                                d_min=3, min_ios=None,
                                                nproc=1, log_out=log_out)
            rb.assign_operators()
        elif len(cell_and_files) > 1: # no need if only 1 data
            rb = resolve_reindex.KabschSelectiveBreeding([x[1] for x in list(cell_and_files.values())], max_delta=5,
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


    if (params.space_group, params.unit_cell).count(None) == 1:
        log_out.write("Error: Specify both space_group and unit_cell!")
        return

    ref_sym_global = None
    if params.space_group is not None:
        try:
            ref_sym_global = crystal.symmetry(params.unit_cell, params.space_group)
            if not ref_sym_global.change_of_basis_op_to_reference_setting().is_identity_op():
                xs_refset = ref_sym_global.as_reference_setting()
                log_out.write('Sorry. Currently space group in non-reference setting is not supported. In this case please give space_group=%s unit_cell="%s" instead.' % (str(xs_refset.space_group_info()).replace(" ",""), format_unit_cell(xs_refset.unit_cell())))
                return
        except:
            log_out.write("Invalid crystal symmetry. Check space_group= and unit_cell=.")
            return

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
        batchjobs = batchjob.ExecLocal(max_parallel=params.batch.sh_max_jobs)
    else:
        batchjobs = None

    jobs = []

    for k in samples:
        params2 = copy.deepcopy(params)
        ref_array = ref_arrays.get(params.reference, None)
        ref_sym = ref_sym_global
        params2.merge.reference.data = params.reference
            

        # Reflect custom parameters
        if "anomalous" in samples[k][1]:
            params2.merge.anomalous = samples[k][1]["anomalous"]
        if "reference" in samples[k][1]:
            ref_array = ref_arrays[samples[k][1]["reference"]]
            params2.merge.reference.data = samples[k][1]["reference"]
        if "reference_sym" in samples[k][1]:
            ref_sym = samples[k][1]["reference_sym"]

        if params2.merge.reference.data:
            params2.merge.reference.data = os.path.abspath(params2.merge.reference.data)

        log_out.write("\n\n")
        workdir = os.path.join(params2.workdir, "%s%s" % (params2.prefix,
                                                         replace_forbidden_chars(k).replace(" ","_")))
        os.mkdir(workdir)
        if batchjobs:
            shname = "multimerge.sh"
            pickle.dump(dict(workdir=workdir, topdirs=samples[k][0],
                             cell_method=params2.cell_method, ref_array=ref_array, ref_sym=ref_sym,
                             merge_params=params2.merge, rescut_params=params2.rescut, filter_params=params2.filtering),
                        open(os.path.join(workdir, "kwargs.pkl"), "wb"), -1)
            job = batchjob.Job(workdir, shname, nproc=params2.batch.nproc_each)
            job.write_script("""\
cd "%s" || exit 1
"%s" -c '\
import pickle; \
from yamtbx.dataproc.auto.command_line.auto_multi_merge import auto_merge; \
kwargs = pickle.load(open("kwargs.pkl", "rb")); \
auto_merge(**kwargs); \
'
""" % (os.path.abspath(workdir), sys.executable))
            batchjobs.submit(job)
            jobs.append(job)
        else:
            try:
                auto_merge(workdir=workdir, topdirs=samples[k][0],
                           cell_method=params2.cell_method, ref_array=ref_array, ref_sym=ref_sym,
                           merge_params=params2.merge, rescut_params=params2.rescut, filter_params=params2.filtering,
                           log_out_all=log_out)
            except:
                log_out.write("Error occurred in %s\n%s\n"%(workdir, traceback.format_exc()))

        log_out.flush()

    if batchjobs:
        batchjobs.wait_all(jobs)

# run()
def run_from_args(argv):
    if "-h" in argv or "--help" in argv or not argv:
        print("All parameters:\n")
        iotbx.phil.parse(gui_phil_str).show(prefix="  ", attributes_level=1)
        print()
        print()
        print(r"""Typical usage:

kamo.auto_multi_merge.dev \
  csv=automerge.csv \
  workdir=$PWD \
  prefix=merge_ \
  cell_method=reindex \
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
""")

        return

    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=gui_phil_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    if params.workdir is None:
        params.workdir = os.getcwd()

    for arg in args:
        if not os.path.exists(arg):
            print("Error: Given path does not exist: %s" % arg)
            return
        if params.csv is None and arg.endswith(".csv"):
            params.csv = arg

    run(params)
# run_from_args()

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
